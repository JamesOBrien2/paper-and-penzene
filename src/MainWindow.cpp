#include "MainWindow.h"
#include "Canvas.h"
#include "Chem.h"

#include <QActionGroup>
#include <QApplication>
#include <QClipboard>
#include <QCloseEvent>
#include <QColorDialog>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QSettings>
#include <QStyle>
#include <QStyleHints>
#include <QPainter>
#include <QtMath>
#include <functional>
#include <QLabel>
#include <QRegularExpression>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QStatusBar>
#include <QToolBar>
#include <QUndoStack>

static const char* kMolMime = "chemical/x-mdl-molfile";
static const char* kPenzMime = "application/x-penzene";  // full fidelity: arrows and text too

MainWindow::MainWindow() : undo_(new QUndoStack(this)), canvas_(new Canvas(undo_, this)) {
    setCentralWidget(canvas_);
    setWindowTitle("Penzene");
    resize(1100, 750);
    buildTools();
    buildMenus();
    connect(undo_, &QUndoStack::cleanChanged, this, &MainWindow::updateTitle);
    updateTitle();
    info_ = new QLabel;
    info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    statusBar()->addPermanentWidget(info_);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateInfo);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateInfo);
}

// Formula and masses of the selection, or of everything.
void MainWindow::updateInfo() {
    auto p = chem::properties(canvas_->selectedSubset());
    if (!p) return info_->clear();
    QString f = QString::fromStdString(p->formula).toHtmlEscaped();
    f.replace(QRegularExpression("(\\d+)"), "<sub>\\1</sub>").replace(QRegularExpression("([+-])$"), "<sup>\\1</sup>");
    info_->setText(tr("%1 &nbsp;·&nbsp; MW %2 &nbsp;·&nbsp; exact mass %3")
                       .arg(f)
                       .arg(p->mw, 0, 'f', 2)
                       .arg(p->exactMass, 0, 'f', 4));
}

void MainWindow::updateTitle() {
    QString name = path_.isEmpty() ? tr("Untitled") : QFileInfo(path_).fileName();
    setWindowTitle(name + "[*] — Penzene");
    setWindowModified(!undo_->isClean());
}

bool MainWindow::openFile(const QString& path) {
    const QString ext = QFileInfo(path).suffix().toLower();
    std::optional<Document> doc = chem::readFile(path);
    if (!doc) {
        QMessageBox::warning(this, tr("Open"), tr("%1 is not a structure file I can read.").arg(path));
        return false;
    }
    undo_->clear();
    canvas_->setDocumentSilently(*doc);
    canvas_->fitToDocument();
    path_ = ext == "cdxml" || ext == "cdx" ? QString() : path;  // never save over a ChemDraw file
    updateTitle();
    return true;
}

bool MainWindow::saveTo(const QString& path) {
    const auto& doc = canvas_->document();
    QByteArray data = path.endsWith(".penz", Qt::CaseInsensitive)
                          ? doc.toJson()
                          : QByteArray::fromStdString(chem::toMolBlock(doc));
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly) || f.write(data) != data.size()) {
        QMessageBox::warning(this, tr("Save"), tr("Cannot write %1").arg(path));
        return false;
    }
    path_ = path;
    undo_->setClean();
    updateTitle();
    return true;
}

bool MainWindow::save() {
    // MOL can't hold everything .penz will (text, arrows), so only .penz saves silently.
    return path_.endsWith(".penz", Qt::CaseInsensitive) ? saveTo(path_) : saveAs();
}

bool MainWindow::saveAs() {
    QString path = QFileDialog::getSaveFileName(this, tr("Save As"), path_,
                                                tr("Penzene document (*.penz);;MDL Molfile (*.mol)"));
    return !path.isEmpty() && saveTo(path);
}

bool MainWindow::maybeSave() {
    if (undo_->isClean()) return true;
    auto r = QMessageBox::question(this, tr("Unsaved changes"), tr("Save changes to this document?"),
                                   QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    return r == QMessageBox::Discard || (r == QMessageBox::Save && save());
}

void MainWindow::closeEvent(QCloseEvent* e) {
    if (maybeSave()) e->accept();
    else e->ignore();
}

void MainWindow::exportImage() {
    QString base = path_.isEmpty() ? QString("structure") : QFileInfo(path_).completeBaseName();
    QString path = QFileDialog::getSaveFileName(this, tr("Export"), base + ".svg",
                                                tr("SVG (*.svg);;PNG image (*.png);;PDF (*.pdf)"));
    if (path.isEmpty()) return;
    if (!exportDocument(canvas_->selectedSubset(), path))
        QMessageBox::warning(this, tr("Export"), tr("Nothing to export, or cannot write %1").arg(path));
}

void MainWindow::importSmiles() {
    bool ok = false;
    QString s = QInputDialog::getText(this, tr("Import SMILES"), tr("SMILES:"), QLineEdit::Normal, {}, &ok);
    if (!ok || s.trimmed().isEmpty()) return;
    if (auto doc = chem::fromSmiles(s.trimmed().toStdString())) canvas_->insert(*doc, tr("Import SMILES"));
    else QMessageBox::warning(this, tr("Import SMILES"), tr("Not a valid SMILES string."));
}

void MainWindow::copy() {
    Document doc = canvas_->selectedSubset();
    if (doc.empty()) return;
    auto* mime = new QMimeData;
    mime->setImageData(renderImage(doc));
    mime->setData("image/svg+xml", renderSvg(doc));
    mime->setData(kPenzMime, doc.toJson());
    if (!doc.atoms.empty()) {
        std::string mol = chem::toMolBlock(doc), smi = chem::toSmiles(doc);
        mime->setData(kMolMime, QByteArray::fromStdString(mol));
        mime->setText(QString::fromStdString(smi.empty() ? mol : smi));
    }
    QApplication::clipboard()->setMimeData(mime);
}

void MainWindow::paste() {
    const QMimeData* mime = QApplication::clipboard()->mimeData();
    if (auto doc = Document::fromJson(mime->data(kPenzMime)); doc && !doc->empty())
        return canvas_->insert(*doc, tr("Paste"));
    std::string text = mime->hasFormat(kMolMime) ? mime->data(kMolMime).toStdString()
                                                 : mime->text().trimmed().toStdString();
    if (text.empty()) return;
    auto doc = text.find("M  END") != std::string::npos ? chem::fromMolBlock(text) : chem::fromSmiles(text);
    if (doc) canvas_->insert(*doc, tr("Paste"));
    else statusBar()->showMessage(tr("Clipboard has no structure or SMILES"), 4000);
}

// Tool icons are drawn with the same renderer as the canvas, in the palette's ink.
// Returned as makers so they can be repainted when the theme changes.
using IconMaker = std::function<QIcon()>;
static IconMaker paintedIcon(std::function<void(QPainter&, QColor)> paint) {
    return [paint] {
    QPixmap pm(48, 48);
    pm.setDevicePixelRatio(2);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setRenderHint(QPainter::Antialiasing);
    paint(p, QApplication::palette().color(QPalette::WindowText));
    return QIcon(pm);
    };
}

static IconMaker docIcon(const Document& d) {
    return paintedIcon([d](QPainter& p, QColor ink) {
        QRectF r;
        for (const auto& a : d.atoms) r |= QRectF(a.pos, QSizeF(0.01, 0.01));
        for (const auto& a : d.arrows) r |= arrowPath(a).boundingRect().adjusted(-3, -3, 3, 3);
        for (const auto& t : d.texts) r |= textPath(t).boundingRect();
        const double s = 18 / std::max(r.width(), r.height());
        p.translate(12, 12);
        p.scale(s, s);
        p.translate(-r.center());
        paintDocument(p, d, {ink, ink, 1.3 / s});  // constant stroke whatever the scale
    });
}

static Document chainDoc(std::vector<QPointF> pts, int order = 1, BondStereo stereo = BondStereo::None) {
    Document d;
    for (QPointF q : pts) d.addAtom(q * kBondLength);
    for (int i = 1; i < int(pts.size()); ++i) d.bonds.push_back({i - 1, i, i == 1 ? order : 1, i == 1 ? stereo : BondStereo::None});
    return d;
}

static Document ringDoc(int n, bool aromatic) {
    Document d;
    const double r = kBondLength / (2 * std::sin(M_PI / n));
    for (int k = 0; k < n; ++k) d.addAtom(r * QPointF(std::sin(2 * M_PI * k / n), -std::cos(2 * M_PI * k / n)));
    for (int k = 0; k < n; ++k) d.bonds.push_back({k, (k + 1) % n, aromatic && k % 2 == 0 ? 2 : 1});
    return d;
}

static Document arrowDoc(ArrowKind kind, double bend = 0) {
    Document d;
    d.arrows.push_back({{0, 0}, {16, 0}, kind, bend});
    return d;
}

static Document textDoc(const QString& s) {
    Document d;
    d.texts.push_back({{0, 0}, s});
    return d;
}

void MainWindow::buildTools() {
    auto* bar = addToolBar(tr("Tools"));
    bar->setObjectName("tools");
    addToolBar(Qt::LeftToolBarArea, bar);
    bar->setIconSize({22, 22});
    auto* group = new QActionGroup(this);

    auto add = [&](const IconMaker& icon, const QString& tip, auto setup) {
        auto* a = bar->addAction(icon(), {});
        icons_.push_back({a, icon});
        a->setToolTip(tip);
        a->setCheckable(true);
        group->addAction(a);
        connect(a, &QAction::triggered, this, setup);
        return a;
    };
    using T = Canvas::Tool;
    auto tool = [this](T t) { return [this, t] { canvas_->setTool(t); }; };
    auto bond = [this](int order) {
        return [this, order] { canvas_->setTool(T::Bond), canvas_->setBondOrder(order); };
    };
    // Keys that pick a tool when no atom or bond is the hotspot (ChemDraw).
    QHash<QString, QAction*> keys;
    const IconMaker select = paintedIcon([](QPainter& p, QColor ink) {
        p.setPen(QPen(ink, 1.2, Qt::DashLine));
        p.drawRect(QRectF(4.5, 5.5, 15, 13));
    });
    keys[" "] = add(select, tr("Select (drag to move, Alt+drag to rotate, double-click for fragment) — Space"),
                    tool(T::Select));
    const QPointF bondPts[] = {{0, 0}, {0.87, -0.5}};
    auto bondIcon = [&](int order, BondStereo st = BondStereo::None) {
        return docIcon(chainDoc({std::begin(bondPts), std::end(bondPts)}, order, st));
    };
    keys["x"] = add(bondIcon(1), tr("Single bond / chain start — x"), bond(1));
    keys["x"]->setChecked(true);
    add(bondIcon(2), tr("Double bond"), bond(2));
    add(bondIcon(3), tr("Triple bond"), bond(3));
    add(bondIcon(1, BondStereo::Wedge), tr("Wedge bond"), tool(T::Wedge));
    add(bondIcon(1, BondStereo::Hash), tr("Hashed bond"), tool(T::Hash));
    keys["X"] = add(docIcon(chainDoc({{0, 0}, {0.87, -0.5}, {1.73, 0}, {2.6, -0.5}})), tr("Chain — X"), tool(T::Chain));
    bar->addSeparator();

    auto ring = [this](int n, bool arom) {
        return [this, n, arom] { canvas_->setTool(T::Ring), canvas_->setRing(n, arom); };
    };
    keys["j"] = add(docIcon(ringDoc(6, true)), tr("Benzene — j"), ring(6, true));
    for (int n = 3; n <= 8; ++n) add(docIcon(ringDoc(n, false)), tr("%1-membered ring").arg(n), ring(n, false));
    Document filled = ringDoc(6, false);
    filled.fills.push_back({{0, 1, 2, 3, 4, 5}, QColor(120, 170, 255)});
    add(docIcon(filled), tr("Ring fill (click inside a ring; again to clear) — colour in Structure menu"), tool(T::Fill));
    bar->addSeparator();

    auto* elements = new QComboBox;
    for (auto s : {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H", "B", "Si"}) elements->addItem(s);
    elements->setToolTip(tr("Element for the atom tool (or hover an atom and press C, N, O…)"));
    auto* atom = add(docIcon(textDoc("N")), tr("Atom"), tool(T::Atom));
    bar->addWidget(elements);
    connect(elements, &QComboBox::currentTextChanged, this, [this, atom](const QString& s) {
        canvas_->setElement(chem::atomicNumber(s.toStdString()));
        canvas_->setTool(T::Atom);
        atom->setChecked(true);
    });
    bar->addSeparator();
    auto arrow = [this](ArrowKind k, bool curved) {
        return [this, k, curved] { canvas_->setTool(T::Arrow), canvas_->setArrow(k, curved); };
    };
    const QString drag = tr(" (drag to draw; click an arrow to restyle it)");
    keys["e"] = add(docIcon(arrowDoc(ArrowKind::Reaction)), tr("Reaction arrow — e") + drag,
                    arrow(ArrowKind::Reaction, false));
    add(docIcon(arrowDoc(ArrowKind::Equilibrium)), tr("Equilibrium arrow") + drag, arrow(ArrowKind::Equilibrium, false));
    add(docIcon(arrowDoc(ArrowKind::Resonance)), tr("Resonance arrow") + drag, arrow(ArrowKind::Resonance, false));
    add(docIcon(arrowDoc(ArrowKind::Retro)), tr("Retrosynthesis arrow") + drag, arrow(ArrowKind::Retro, false));
    add(docIcon(arrowDoc(ArrowKind::Reaction, 10)), tr("Curved arrow, electron pair (click it again to flip the curve)"),
        arrow(ArrowKind::Reaction, true));
    add(docIcon(arrowDoc(ArrowKind::Fishhook, 10)), tr("Fishhook arrow, single electron (click it again to flip)"),
        arrow(ArrowKind::Fishhook, true));
    keys["t"] = add(docIcon(textDoc("T")), tr("Text (click to add or edit; H2O is set as H₂O) — t"), tool(T::Text));
    connect(canvas_, &Canvas::toolKey, this, [keys](const QString& k) {
        if (auto* a = keys.value(k)) a->trigger();
    });
    bar->addSeparator();
    auto charge = [](bool plus) {
        return paintedIcon([plus](QPainter& p, QColor ink) {
            p.setPen(QPen(ink, 1.3));
            p.drawEllipse(QPointF(12, 12), 7, 7);
            p.drawLine(QPointF(8.5, 12), QPointF(15.5, 12));
            if (plus) p.drawLine(QPointF(12, 8.5), QPointF(12, 15.5));
        });
    };
    add(charge(true), tr("Positive charge"), tool(T::ChargePlus));
    add(charge(false), tr("Negative charge"), tool(T::ChargeMinus));
    const IconMaker eraser = paintedIcon([](QPainter& p, QColor ink) {
        p.translate(12, 12);
        p.rotate(-40);
        p.setPen(QPen(ink, 1.3));
        p.drawRoundedRect(QRectF(-8, -4, 16, 8), 1.5, 1.5);
        p.drawLine(QPointF(-2, -4), QPointF(-2, 4));
    });
    add(eraser, tr("Eraser"), tool(T::Erase));
}

// "System" follows the OS; the others force light or dark, and Catppuccin
// also recolours the UI. Canvas colours always come from the theme.
void MainWindow::applyTheme(const QString& name) {
    const Theme& chosen = theme(name);
    auto* hints = QGuiApplication::styleHints();
    hints->setColorScheme(chosen.name == "System" ? Qt::ColorScheme::Unknown
                          : chosen.dark           ? Qt::ColorScheme::Dark
                                                  : Qt::ColorScheme::Light);
    QPalette pal = QApplication::style()->standardPalette();
    if (chosen.window.isValid()) {
        for (auto role : {QPalette::Window, QPalette::Button}) pal.setColor(role, chosen.window);
        for (auto role : {QPalette::WindowText, QPalette::Text, QPalette::ButtonText, QPalette::ToolTipText})
            pal.setColor(role, chosen.text);
        pal.setColor(QPalette::Base, chosen.paper);
        pal.setColor(QPalette::AlternateBase, chosen.surface);
        pal.setColor(QPalette::ToolTipBase, chosen.surface);
        pal.setColor(QPalette::Highlight, chosen.accent);
        pal.setColor(QPalette::HighlightedText, chosen.paper);
        pal.setColor(QPalette::Mid, chosen.surface);
        QApplication::setPalette(pal);
    } else {
        QApplication::setPalette(QPalette());  // back to the platform's own
    }
    const bool dark = chosen.name == "System" ? hints->colorScheme() == Qt::ColorScheme::Dark : chosen.dark;
    canvas_->setTheme(chosen.name == "System" ? theme(dark ? "Dark" : "Light") : chosen);
    for (auto& [action, make] : icons_) action->setIcon(make());
    QSettings().setValue("theme", chosen.name);
}

void MainWindow::buildMenus() {
    auto* file = menuBar()->addMenu(tr("&File"));
    file->addAction(tr("&New"), QKeySequence::New, this, [this] {
        if (!maybeSave()) return;
        undo_->clear();
        canvas_->setDocumentSilently({});
        path_.clear();
        updateTitle();
    });
    file->addAction(tr("&Open…"), QKeySequence::Open, this, [this] {
        if (!maybeSave()) return;
        QString p = QFileDialog::getOpenFileName(this, tr("Open"), {},
                                                 tr("Structures (*.penz *.mol *.sdf *.cdxml *.cdx);;All files (*)"));
        if (!p.isEmpty()) openFile(p);
    });
    file->addAction(tr("&Save"), QKeySequence::Save, this, &MainWindow::save);
    file->addAction(tr("Save &As…"), QKeySequence::SaveAs, this, &MainWindow::saveAs);
    file->addSeparator();
    file->addAction(tr("Import &SMILES…"), QKeySequence(tr("Ctrl+Shift+I")), this, &MainWindow::importSmiles);
    file->addAction(tr("&Export…"), QKeySequence(tr("Ctrl+E")), this, &MainWindow::exportImage);
    file->addSeparator();
    file->addAction(tr("&Quit"), QKeySequence::Quit, this, &QWidget::close);

    auto* edit = menuBar()->addMenu(tr("&Edit"));
    auto* u = undo_->createUndoAction(this);
    u->setShortcut(QKeySequence::Undo);
    auto* r = undo_->createRedoAction(this);
    r->setShortcut(QKeySequence::Redo);
    edit->addAction(u);
    edit->addAction(r);
    edit->addSeparator();
    edit->addAction(tr("Cu&t"), QKeySequence::Cut, this, [this] {
        copy();
        canvas_->deleteSelection();
    });
    edit->addAction(tr("&Copy"), QKeySequence::Copy, this, &MainWindow::copy);
    edit->addAction(tr("Copy as S&MILES"), QKeySequence(tr("Ctrl+Alt+C")), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toSmiles(canvas_->selectedSubset())));
    });
    edit->addAction(tr("Copy as &InChI"), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toInchi(canvas_->selectedSubset())));
    });
    edit->addAction(tr("Copy as InChI&Key"), this, [this] {
        QApplication::clipboard()->setText(QString::fromStdString(chem::toInchiKey(canvas_->selectedSubset())));
    });
    edit->addAction(tr("&Paste"), QKeySequence::Paste, this, &MainWindow::paste);
    edit->addAction(tr("&Delete"), canvas_, &Canvas::deleteSelection);
    edit->addSeparator();
    edit->addAction(tr("Select &All"), QKeySequence::SelectAll, canvas_, &Canvas::selectAll);

    auto* structure = menuBar()->addMenu(tr("&Structure"));
    // ponytail: cleans the whole document; clean just the selection when someone asks.
    structure->addAction(tr("&Clean Structure"), QKeySequence(tr("Ctrl+Shift+K")), this, [this] {
        const auto& sel = canvas_->selection();  // selected molecules only, else everything
        canvas_->commit(chem::clean2D(canvas_->document(), {sel.begin(), sel.end()}), tr("Clean"));
    });
    // Drawing style presets, like ChemDraw's document settings; stored in the .penz.
    auto* styles = structure->addMenu(tr("Drawing &Style"));
    auto* styleGroup = new QActionGroup(this);
    for (const auto& st : drawingStyles()) {
        auto* act = styles->addAction(st.name);
        act->setCheckable(true);
        styleGroup->addAction(act);
        connect(act, &QAction::triggered, this, [this, name = st.name] {
            Document next = canvas_->document();
            next.style = name == drawingStyles()[0].name ? QString() : name;
            if (!(next == canvas_->document())) canvas_->commit(next, tr("Drawing style"));
        });
    }
    auto syncStyle = [this, styleGroup] {
        const QString current = drawingStyle(canvas_->document().style).name;
        for (auto* act : styleGroup->actions()) act->setChecked(act->text() == current);
    };
    connect(canvas_, &Canvas::documentChanged, this, syncStyle);
    syncStyle();
    structure->addAction(tr("Ring &Fill Colour…"), this, [this] {
        QColor c = QColorDialog::getColor(canvas_->fillColor(), this, tr("Ring fill colour"));
        if (c.isValid()) canvas_->setFillColor(c);
    });
    structure->addAction(tr("&Expand Abbreviations"), QKeySequence(tr("Ctrl+Shift+E")), canvas_,
                         &Canvas::expandAbbreviations);

    auto* view = menuBar()->addMenu(tr("&View"));
    view->addAction(tr("Zoom &In"), QKeySequence::ZoomIn, this, [this] { canvas_->zoomBy(1.25); });
    view->addAction(tr("Zoom &Out"), QKeySequence::ZoomOut, this, [this] { canvas_->zoomBy(0.8); });
    view->addAction(tr("&Fit to Window"), QKeySequence(tr("Ctrl+0")), canvas_, &Canvas::fitToDocument);
    auto* themeMenu = view->addMenu(tr("&Theme"));
    auto* themeGroup = new QActionGroup(this);
    const QString current = QSettings().value("theme", "System").toString();
    for (const auto& t : themes()) {
        auto* a = themeMenu->addAction(t.name, this, [this, n = t.name] { applyTheme(n); });
        a->setCheckable(true);
        a->setChecked(t.name == theme(current).name);
        themeGroup->addAction(a);
        if (t.name == "Dark") themeMenu->addSeparator();
    }
    applyTheme(current);
    // Following the OS: repaint the canvas and icons when it switches.
    connect(QGuiApplication::styleHints(), &QStyleHints::colorSchemeChanged, this, [this] {
        if (QSettings().value("theme", "System").toString() == "System") applyTheme("System");
    });

    auto* help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&Keyboard Shortcuts"), QKeySequence(tr("F1")), this, [this] {
        QMessageBox box(this);
        box.setWindowTitle(tr("Keyboard Shortcuts"));
        box.setTextFormat(Qt::RichText);
        box.setText(tr(R"(<p>Point at an atom or bond to make it the <b>hotspot</b>. It stays put when the mouse
moves off, so you can keep typing. Follows ChemDraw's hotkeys.</p>
<table cellspacing="5">
<tr><th colspan="2" align="left">Moving the hotspot</th></tr>
<tr><td><b>←↑→↓</b></td><td>atom → bond → atom; with <b>Shift</b>: atom → atom, bond → bond</td></tr>
<tr><td><b>Esc</b></td><td>clear hotspot and selection</td></tr>
<tr><th colspan="2" align="left">Atom: sprout</th></tr>
<tr><td><b>1</b> / <b>0</b></td><td>single bond, linear / cyclic mode (0 is longer on 2°/3° carbons)</td></tr>
<tr><td><b>2</b></td><td>acetyl (1°), C=O (2°), CH<sub>2</sub>-acetyl (3°/aromatic)</td></tr>
<tr><td><b>3</b> or <b>a</b></td><td>phenyl</td></tr>
<tr><td><b>4</b> / <b>5</b></td><td>wedged / hashed methyl</td></tr>
<tr><td><b>6 7 u v</b></td><td>cyclohexane, cyclopentane, cyclobutane, cyclopropane (spiro on 2°)</td></tr>
<tr><td><b>8 9 z</b></td><td>methylidene, dimethyl / gem-dimethyl / isopropyl, alkyne</td></tr>
<tr><td><b>k K</b></td><td>sulfonyl, t-Bu</td></tr>
<tr><th colspan="2" align="left">Atom: label</th></tr>
<tr><td><b>c n/w o/q s p f l b i h</b></td><td>C N O S P F Cl Br I H</td></tr>
<tr><td><b>B S L</b></td><td>B, Si, Li</td></tr>
<tr><td><b>m e P A</b></td><td>Me, Et, Ph, Ac</td></tr>
<tr><td><b>O N F E Z</b></td><td>OMe, NO<sub>2</sub>, CF<sub>3</sub>, CO<sub>2</sub>Me, N<sub>3</sub></td></tr>
<tr><td><b>Y H Q M</b></td><td>Boc, Cbz, Fmoc, MgBr</td></tr>
<tr><td><b>+ −</b></td><td>charge</td></tr>
<tr><td><b>Enter</b> or <b>=</b></td><td>type a label: element, abbreviation (OMe, Boc, TBS…) or SMILES</td></tr>
<tr><td><b>Delete</b></td><td>remove label (C stays), or delete a carbon</td></tr>
<tr><th colspan="2" align="left">Bond</th></tr>
<tr><td><b>1 2 3</b></td><td>single, double, triple; <b>2</b> on a double bond swaps the side of its second line</td></tr>
<tr><td><b>w</b> / <b>h</b></td><td>wedged / hashed (press again to flip)</td></tr>
<tr><td><b>a z</b></td><td>fuse benzene / cyclopentadiene</td></tr>
<tr><td><b>v 4–8</b></td><td>fuse ring of that size (v = 3)</td></tr>
<tr><td><b>9</b> / <b>0</b></td><td>fuse chair cyclohexane (two orientations)</td></tr>
<tr><td><b>d b y</b></td><td>dashed, bold, wavy</td></tr>
<tr><td><b>D</b> / <b>B</b></td><td>dashed double / bold double</td></tr>
<tr><td><b>l c r</b></td><td>double bond's second line left / centred / right</td></tr>
<tr><th colspan="2" align="left">No hotspot (Esc)</th></tr>
<tr><td><b>x X j e t Space</b></td><td>bond, chain, benzene, arrow, text, select tool</td></tr>
<tr><th colspan="2" align="left">Selection</th></tr>
<tr><td><b>Ctrl+←↑→↓</b></td><td>duplicate across the next arrow that way (or alongside)</td></tr>
<tr><td><b>Alt+← →</b></td><td>rotate 15° &nbsp;•&nbsp; <b>Alt+drag</b> rotate freely • <b>double-click</b> select fragment, or edit text</td></tr>
</table>)"));
        box.exec();
    });
    help->addAction(tr("&About Penzene"), this, [this] {
        QMessageBox::about(this, tr("About Penzene"),
                           tr("<h3>Penzene %1</h3><p>An open-source chemical structure editor.</p>"
                              "<p>GPL-3.0 • <a href='https://github.com/JamesOBrien2/penzene'>GitHub</a></p>"
                              "<p>Chemistry by RDKit. GUI by Qt.</p>").arg(PENZENE_VERSION));
    });
}
