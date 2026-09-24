#include "MainWindow.h"
#include "Canvas.h"
#include "Chem.h"
#include "PubChem.h"

#include <QActionGroup>
#include <QApplication>
#include <QClipboard>
#include <QCloseEvent>
#include <QColorDialog>
#include <QFile>
#include <QComboBox>
#include <QFileDialog>
#include <QSpinBox>
#include <QFormLayout>
#include <QVBoxLayout>
#include <QListWidget>
#include <QDialogButtonBox>
#include <QDialog>
#include <QVBoxLayout>
#include <QPushButton>
#include <QDockWidget>
#include <QFileInfo>
#include <QInputDialog>
#include <memory>
#include <QSet>
#include <QWidgetAction>
#include <QToolButton>
#include <QMenu>
#include <QGridLayout>
#include <QFrame>
#include <QDir>
#include <QSettings>
#include <QSignalBlocker>
#include <algorithm>
#include <QStandardPaths>
#include <QTimer>
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
    setWindowTitle("Penzene " PENZENE_BUILD);
    resize(1100, 750);
    // Properties panel (built before the menus, which offer its toggle): descriptors for the selection or everything.
    profileDock_ = new QDockWidget(tr("Properties"), this);
    profileDock_->setObjectName("properties");
    auto* panel = new QWidget;
    auto* panelLayout = new QVBoxLayout(panel);
    profile_ = new QLabel;
    profile_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    profile_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    auto* copyProfile = new QPushButton(tr("Copy as Text"));
    connect(copyProfile, &QPushButton::clicked, this, [this] { QApplication::clipboard()->setText(profileText_); });
    panelLayout->addWidget(profile_);
    panelLayout->addWidget(copyProfile);
    panelLayout->addStretch();
    panel->setMinimumWidth(300);  // room for the values beside their names
    profileDock_->setWidget(panel);
    addDockWidget(Qt::RightDockWidgetArea, profileDock_);
    profileDock_->hide();
    connect(profileDock_, &QDockWidget::visibilityChanged, this, &MainWindow::updateProfile);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateProfile);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateProfile);
    buildTools();
    buildMenus();
    connect(undo_, &QUndoStack::cleanChanged, this, &MainWindow::updateTitle);
    updateTitle();
    info_ = new QLabel;
    info_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    statusBar()->addPermanentWidget(info_);
    connect(canvas_, &Canvas::documentChanged, this, &MainWindow::updateInfo);
    connect(canvas_, &Canvas::selectionChanged, this, &MainWindow::updateInfo);

    auto* autosaver = new QTimer(this);
    connect(autosaver, &QTimer::timeout, this, &MainWindow::autosave);
    autosaver->start(60 * 1000);
}

void MainWindow::updateProfile() {
    if (!profileDock_->isVisible()) return;
    const auto p = chem::profile(canvas_->selectedSubset());
    if (!p) {
        profile_->setText(tr("Draw or select a valid structure."));
        profileText_.clear();
        return;
    }
    QString formula = QString::fromStdString(p->basic.formula).toHtmlEscaped();
    formula.replace(QRegularExpression("(\\d+)"), "<sub>\\1</sub>");
    QStringList analysis, plain;
    for (const auto& [el, pct] : p->elemental) {
        analysis << QString("%1 %2").arg(QString::fromStdString(el)).arg(pct, 0, 'f', 2);
    }
    const QList<std::pair<QString, QString>> rows{
        {tr("Formula"), formula},
        {tr("MW"), QString::number(p->basic.mw, 'f', 2)},
        {tr("Exact mass"), QString::number(p->basic.exactMass, 'f', 4)},
        {tr("Elemental (%)"), analysis.join(", ")},
        {tr("cLogP"), QString::number(p->logP, 'f', 2)},
        {tr("TPSA (Å²)"), QString::number(p->tpsa, 'f', 1)},
        {tr("H-bond donors"), QString::number(p->hbd)},
        {tr("H-bond acceptors"), QString::number(p->hba)},
        {tr("Rotatable bonds"), QString::number(p->rotatable)},
        {tr("Heavy atoms"), QString::number(p->heavyAtoms)},
        {tr("Lipinski (Ro5)"), p->lipinskiViolations ? tr("%n violation(s)", "", p->lipinskiViolations) : tr("passes")},
        {tr("Veber"), p->veber ? tr("passes") : tr("fails")},
    };
    QString html = "<table cellspacing='4'>";
    for (const auto& [k, v] : rows) {
        html += QString("<tr><td><b>%1</b></td><td>%2</td></tr>").arg(k, v);
        plain << k + "\t" + QString(v).remove(QRegularExpression("<[^>]*>"));
    }
    profile_->setText(html + "</table>");
    profileText_ = plain.join("\n");
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

// Children are deleted after this destructor has run, and some signal on the way
// out (the undo stack's cleanChanged, a dock's visibilityChanged as it hides):
// none of that may reach a half-destroyed window.
MainWindow::~MainWindow() {
    for (QObject* child : findChildren<QObject*>()) child->disconnect(this);
}

void MainWindow::updateTitle() {
    QString name = path_.isEmpty() ? tr("Untitled") : QFileInfo(path_).fileName();
    setWindowTitle(name + "[*] — Penzene " PENZENE_BUILD);
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
    remember(path);
    updateTitle();
    return true;
}

QStringList MainWindow::recentFiles() const { return QSettings().value("recentFiles").toStringList(); }

void MainWindow::remember(const QString& path) {
    QStringList files = recentFiles();
    files.removeAll(QFileInfo(path).absoluteFilePath());
    files.prepend(QFileInfo(path).absoluteFilePath());
    QSettings().setValue("recentFiles", files.mid(0, 10));
}

QString MainWindow::autosavePath() {
    return QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/autosave.penz";
}

void MainWindow::autosave() {
    if (undo_->isClean()) return QFile::remove(autosavePath()), void();
    QDir().mkpath(QFileInfo(autosavePath()).path());
    QFile f(autosavePath());
    if (f.open(QIODevice::WriteOnly)) f.write(canvas_->document().toJson());
}

void MainWindow::offerRecovery() {
    QFile f(autosavePath());
    if (!f.open(QIODevice::ReadOnly)) return;
    auto doc = Document::fromJson(f.readAll());
    f.close();
    if (doc && !doc->empty() &&
        QMessageBox::question(this, tr("Recover"),
                              tr("Penzene closed without saving your last drawing. Recover it?")) == QMessageBox::Yes)
        canvas_->commit(*doc, tr("Recover"));  // unsaved, so Save asks where to put it
    QFile::remove(autosavePath());
}

bool MainWindow::saveTo(const QString& path, bool v3000) {
    const auto& doc = canvas_->document();
    QByteArray data = path.endsWith(".penz", Qt::CaseInsensitive)
                          ? doc.toJson()
                          : QByteArray::fromStdString(chem::toMolBlock(doc, v3000));
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly) || f.write(data) != data.size()) {
        QMessageBox::warning(this, tr("Save"), tr("Cannot write %1").arg(path));
        return false;
    }
    path_ = path;
    undo_->setClean();
    QFile::remove(autosavePath());
    remember(path);
    updateTitle();
    return true;
}

bool MainWindow::save() {
    // MOL can't hold everything .penz will (text, arrows), so only .penz saves silently.
    return path_.endsWith(".penz", Qt::CaseInsensitive) ? saveTo(path_) : saveAs();
}

bool MainWindow::saveAs() {
    const QString v3000 = tr("MDL Molfile V3000 (*.mol)");
    QString filter;
    QString path = QFileDialog::getSaveFileName(this, tr("Save As"), path_,
                                                tr("Penzene document (*.penz);;MDL Molfile (*.mol);;") + v3000, &filter);
    return !path.isEmpty() && saveTo(path, filter == v3000);
}

bool MainWindow::maybeSave() {
    if (undo_->isClean()) return true;
    auto r = QMessageBox::question(this, tr("Unsaved changes"), tr("Save changes to this document?"),
                                   QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    return r == QMessageBox::Discard || (r == QMessageBox::Save && save());
}

void MainWindow::closeEvent(QCloseEvent* e) {
    if (!maybeSave()) return e->ignore();
    QFile::remove(autosavePath());  // a deliberate quit: nothing to recover
    e->accept();
}

// Export preferences (Edit > Preferences), used by Export and Copy.
static double exportDpi() { return QSettings().value("exportDpi", 300).toDouble(); }
static QColor exportBackground() {
    return QSettings().value("exportBackground").toString() == "white" ? QColor(Qt::white) : QColor(Qt::transparent);
}

void MainWindow::showPreferences() {
    QDialog dialog(this);
    dialog.setWindowTitle(tr("Preferences"));
    auto* form = new QFormLayout(&dialog);
    auto* themeBox = new QComboBox;
    for (const auto& t : themes()) themeBox->addItem(t.name);
    themeBox->setCurrentText(theme(QSettings().value("theme", "System").toString()).name);
    auto* styleBox = new QComboBox;
    for (const auto& s : drawingStyles()) styleBox->addItem(s.name);
    styleBox->setCurrentText(drawingStyle(QSettings().value("defaultStyle").toString()).name);
    auto* dpiBox = new QSpinBox;
    dpiBox->setRange(72, 1200);
    dpiBox->setSingleStep(50);
    dpiBox->setSuffix(tr(" dpi"));
    dpiBox->setValue(int(exportDpi()));
    auto* backgroundBox = new QComboBox;
    backgroundBox->addItems({tr("Clear"), tr("White")});
    backgroundBox->setCurrentIndex(exportBackground().alpha() ? 1 : 0);
    form->addRow(tr("Theme:"), themeBox);
    form->addRow(tr("Drawing style for new documents:"), styleBox);
    form->addRow(tr("PNG resolution:"), dpiBox);
    form->addRow(tr("Export and copy background:"), backgroundBox);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted) return;
    QSettings settings;
    settings.setValue("defaultStyle", styleBox->currentIndex() ? styleBox->currentText() : QString());
    settings.setValue("exportDpi", dpiBox->value());
    settings.setValue("exportBackground", backgroundBox->currentIndex() ? "white" : "clear");
    applyTheme(themeBox->currentText());
    for (auto* a : themeGroup_->actions()) a->setChecked(a->text() == themeBox->currentText());
}

QWidget* MainWindow::checkStructure() {
    auto* dialog = new QDialog(this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowTitle(tr("Check Structure"));
    auto* layout = new QVBoxLayout(dialog);
    auto* list = new QListWidget;
    layout->addWidget(list);
    const auto problems = chem::checkStructure(canvas_->document());
    for (const auto& p : problems) {
        auto* item = new QListWidgetItem(p.message, list);
        QVariantList atoms;
        for (int i : p.atoms) atoms << i;
        item->setData(Qt::UserRole, atoms);
    }
    if (problems.empty()) list->addItem(tr("No problems found."));
    connect(list, &QListWidget::currentItemChanged, this, [this](QListWidgetItem* item) {
        if (!item) return;
        QSet<int> atoms;
        for (const QVariant& v : item->data(Qt::UserRole).toList()) atoms.insert(v.toInt());
        canvas_->setSelection(atoms);  // show where the problem is
    });
    auto* close = new QDialogButtonBox(QDialogButtonBox::Close);
    connect(close, &QDialogButtonBox::rejected, dialog, &QDialog::close);
    layout->addWidget(close);
    dialog->resize(460, 260);
    dialog->show();
    return dialog;
}

void MainWindow::exportImage() {
    QString base = path_.isEmpty() ? QString("structure") : QFileInfo(path_).completeBaseName();
    QString path = QFileDialog::getSaveFileName(this, tr("Export"), base + ".svg",
                                                tr("SVG (*.svg);;PNG image (*.png);;PDF (*.pdf)"));
    if (path.isEmpty()) return;
    if (!exportDocument(canvas_->selectedSubset(), path, exportDpi(), exportBackground()))
        QMessageBox::warning(this, tr("Export"), tr("Nothing to export, or cannot write %1").arg(path));
}

void MainWindow::importSmiles() {
    bool ok = false;
    QString s = QInputDialog::getText(this, tr("Import SMILES"), tr("SMILES:"), QLineEdit::Normal, {}, &ok);
    if (!ok || s.trimmed().isEmpty()) return;
    if (auto doc = chem::fromSmiles(s.trimmed().toStdString())) canvas_->insert(*doc, tr("Import SMILES"));
    else QMessageBox::warning(this, tr("Import SMILES"), tr("Not a valid SMILES string."));
}

void MainWindow::importName() {
    bool ok = false;
    const QString name = QInputDialog::getText(this, tr("Import Name"),
                                               tr("Compound name (looked up on PubChem, online):"),
                                               QLineEdit::Normal, {}, &ok);
    if (!ok || name.trimmed().isEmpty()) return;
    QString error;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    const QString smiles = pubchem::fetch(pubchem::nameToSmilesUrl(name), "SMILES", &error);
    QApplication::restoreOverrideCursor();
    if (auto doc = chem::fromSmiles(smiles.toStdString()); doc && !smiles.isEmpty())
        canvas_->insert(*doc, tr("Import %1").arg(name.trimmed()));
    else
        QMessageBox::warning(this, tr("Import Name"), tr("Could not look up “%1”: %2").arg(name.trimmed(), error));
}

void MainWindow::copy() {
    Document doc = canvas_->selectedSubset();
    if (doc.empty()) return;
    auto* mime = new QMimeData;
    mime->setImageData(renderImage(doc, exportDpi(), exportBackground()));
    mime->setData("image/svg+xml", renderSvg(doc, exportBackground()));
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
    // A figure Penzene exported, copied from another app or as a file: the drawing inside it.
    for (const char* type : {"image/svg+xml", "image/png"})
        if (auto doc = Document::fromEmbedded(mime->data(type)); doc && !doc->empty())
            return canvas_->insert(*doc, tr("Paste"));
    for (const QUrl& url : mime->urls())
        if (auto doc = url.isLocalFile() ? chem::readFile(url.toLocalFile()) : std::nullopt; doc && !doc->empty())
            return canvas_->insert(*doc, tr("Paste"));
    std::string text = mime->hasFormat(kMolMime) ? mime->data(kMolMime).toStdString()
                                                 : mime->text().trimmed().toStdString();
    if (text.empty()) return;
    auto doc = text.find("M  END") != std::string::npos ? chem::fromMolBlock(text)
               : text.starts_with("InChI=")         ? chem::fromInchi(text)
                                                    : chem::fromSmiles(text);
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
    const double turn = n % 4 == 0 ? M_PI / n : 0;  // squares sit flat, not as diamonds
    for (int k = 0; k < n; ++k)
        d.addAtom(r * QPointF(std::sin(2 * M_PI * k / n + turn), -std::cos(2 * M_PI * k / n + turn)));
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

// Periodic table: main block by group and period, lanthanides and actinides
// underneath. Organic elements are bold, since they're the ones drawn most.
// A drop-down of colour swatches plus "Custom…", for the colour and ring fill tools.
static QMenu* colourMenu(QWidget* parent, const QList<QColor>& presets, std::function<QColor()> current,
                         std::function<void(QColor)> picked) {
    auto* menu = new QMenu(parent);
    auto* w = new QWidget;
    auto* grid = new QGridLayout(w);
    grid->setSpacing(3);
    grid->setContentsMargins(6, 6, 6, 6);
    for (int i = 0; i < presets.size(); ++i) {
        auto* b = new QToolButton;
        b->setFixedSize(24, 24);
        b->setAutoRaise(true);
        b->setToolTip(presets[i].name());
        QPixmap swatch(16, 16);
        swatch.fill(presets[i]);
        b->setIcon(QIcon(swatch));
        QObject::connect(b, &QToolButton::clicked, menu, [=] { picked(presets[i]), menu->close(); });
        grid->addWidget(b, i / 4, i % 4);
    }
    auto* custom = new QToolButton;
    custom->setText(QObject::tr("Custom…"));
    custom->setAutoRaise(true);
    QObject::connect(custom, &QToolButton::clicked, menu, [=] {
        menu->close();
        QColor c = QColorDialog::getColor(current(), parent);
        if (c.isValid()) picked(c);
    });
    grid->addWidget(custom, (presets.size() + 3) / 4, 0, 1, 4);
    auto* action = new QWidgetAction(menu);
    action->setDefaultWidget(w);
    menu->addAction(action);
    return menu;
}

static QWidget* periodicTable(const std::function<void(int)>& picked) {
    auto* w = new QWidget;
    auto* grid = new QGridLayout(w);
    grid->setSpacing(2);
    grid->setContentsMargins(6, 6, 6, 6);
    auto place = [&](int z, int row, int col) {
        const QString sym = QString::fromStdString(chem::symbol(z));
        auto* b = new QToolButton;
        b->setText(sym);
        b->setToolTip(QString("%1 (%2)").arg(sym).arg(z));
        b->setFixedSize(30, 26);
        b->setAutoRaise(true);
        static const QSet<int> organic{1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53};
        if (organic.contains(z)) {
            QFont f = b->font();
            f.setBold(true);
            b->setFont(f);
        }
        QObject::connect(b, &QToolButton::clicked, w, [picked, z] { picked(z); });
        grid->addWidget(b, row, col);
    };
    place(1, 0, 0), place(2, 0, 17);
    for (int p = 1, z = 3; p <= 2; ++p) {  // periods 2-3: s block, then p block
        place(z++, p, 0), place(z++, p, 1);
        for (int c = 12; c < 18; ++c) place(z++, p, c);
    }
    for (int p = 3, z = 19; p <= 4; ++p)  // periods 4-5 are full
        for (int c = 0; c < 18; ++c) place(z++, p, c);
    for (int p = 5, z = 55; p <= 6; ++p, z += 32) {  // periods 6-7: Cs/Fr, Ba/Ra, then Hf/Rf onwards
        place(z, p, 0), place(z + 1, p, 1);
        for (int c = 3; c < 18; ++c) place(z + 14 + c, p, c);  // 72 (Hf) at column 3
        for (int k = 0; k < 15; ++k) place(z + 2 + k, p + 3, 2 + k);  // La-Lu, Ac-Lr below
    }
    grid->setRowMinimumHeight(7, 8);  // gap above the f block
    return w;
}

void MainWindow::buildTools() {
    auto* bar = addToolBar(tr("Tools"));
    bar->setObjectName("tools");
    addToolBar(Qt::LeftToolBarArea, bar);
    bar->setMovable(false);
    // A two-column palette, like ChemDraw's, so related tools sit together.
    auto* palette = new QWidget;
    auto* grid = new QGridLayout(palette);
    grid->setSpacing(2);
    grid->setContentsMargins(4, 4, 4, 4);
    bar->addWidget(palette);
    int slot = 0;  // next free cell, counted left to right
    auto section = [&] {
        if (slot % 2) ++slot;
        auto* line = new QFrame;
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        grid->addWidget(line, slot / 2, 0, 1, 2);
        slot += 2;
    };
    auto* group = new QActionGroup(this);
    auto button = [&](QAction* a) {
        auto* b = new QToolButton;
        b->setDefaultAction(a);
        b->setIconSize({26, 26});
        b->setAutoRaise(true);
        grid->addWidget(b, slot / 2, slot % 2);
        ++slot;
        return b;
    };
    auto add = [&](const IconMaker& icon, const QString& tip, auto setup) {
        auto* a = new QAction(icon(), {}, this);
        icons_.push_back({a, icon});
        a->setToolTip(tip);
        a->setStatusTip(tip);  // the status bar explains the tool while it's chosen
        a->setCheckable(true);
        group->addAction(a);
        connect(a, &QAction::triggered, this, setup);
        connect(a, &QAction::triggered, this, [this, tip] { statusBar()->showMessage(tip); });
        button(a);
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
    const IconMaker eraser = paintedIcon([](QPainter& p, QColor ink) {
        p.translate(12, 12);
        p.rotate(-40);
        p.setPen(QPen(ink, 1.3));
        p.drawRoundedRect(QRectF(-8, -4, 16, 8), 1.5, 1.5);
        p.drawLine(QPointF(-2, -4), QPointF(-2, 4));
    });
    add(eraser, tr("Eraser (click an atom, bond, arrow or text)"), tool(T::Erase));
    section();
    const QPointF bondPts[] = {{0, 0}, {0.87, -0.5}};
    auto bondIcon = [&](int order, BondStereo st = BondStereo::None) {
        return docIcon(chainDoc({std::begin(bondPts), std::end(bondPts)}, order, st));
    };
    keys["x"] = add(bondIcon(1), tr("Single bond — x: click empty space or an atom to add a bond; drag to aim it; click a bond to change it"), bond(1));
    keys["x"]->setChecked(true);
    add(bondIcon(2), tr("Double bond: click an atom to add one, or a bond to make it double"), bond(2));
    add(bondIcon(3), tr("Triple bond: click an atom to add one, or a bond to make it triple"), bond(3));
    add(bondIcon(1, BondStereo::Wedge), tr("Wedge bond: points from the atom you start at; click a wedge again to flip it"), tool(T::Wedge));
    add(bondIcon(1, BondStereo::Hash), tr("Hashed bond: points from the atom you start at; click a hash again to flip it"), tool(T::Hash));
    keys["X"] = add(docIcon(chainDoc({{0, 0}, {0.87, -0.5}, {1.73, 0}, {2.6, -0.5}})), tr("Chain — X: drag to draw a zig-zag chain; it grows with the drag"), tool(T::Chain));
    section();

    auto ring = [this](int n, bool arom) {
        return [this, n, arom] { canvas_->setTool(T::Ring), canvas_->setRing(n, arom); };
    };
    keys["j"] = add(docIcon(ringDoc(6, true)), tr("Benzene — j: click empty space for a ring, an atom to attach one, or a bond to fuse one"), ring(6, true));
    for (int n = 3; n <= 8; ++n) add(docIcon(ringDoc(n, false)), tr("%1-membered ring: click empty space, an atom (spiro/attached) or a bond (fused)").arg(n), ring(n, false));
    // Ring fill: the icon shows the current fill colour; the arrow picks it.
    auto fill = std::make_shared<QColor>(canvas_->fillColor());
    const IconMaker fillIcon = [fill] {
        Document filled = ringDoc(6, false);
        filled.fills.push_back({{0, 1, 2, 3, 4, 5}, *fill});
        return docIcon(filled)();
    };
    auto* fillTool = add(fillIcon, tr("Ring fill: click inside a ring to shade it (again to clear); pick the colour from the arrow"),
                         tool(T::Fill));
    for (auto* b : palette->findChildren<QToolButton*>())
        if (b->defaultAction() == fillTool) {
            b->setPopupMode(QToolButton::MenuButtonPopup);
            // Light tints, so bonds and labels stay readable on top.
            b->setMenu(colourMenu(b,
                                  {QColor(207, 227, 255), QColor(255, 214, 214), QColor(212, 240, 210), QColor(255, 236, 196),
                                   QColor(232, 218, 250), QColor(255, 222, 240), QColor(220, 220, 220), QColor(255, 250, 200)},
                                  [this] { return canvas_->fillColor(); },
                                  [=, this](QColor c) {
                                      *fill = c;
                                      canvas_->setFillColor(c);
                                      canvas_->setTool(T::Fill);
                                      fillTool->setChecked(true);
                                      fillTool->setIcon(fillIcon());
                                  }));
        }
    section();

    // Element: the button shows the current element and draws it; its arrow
    // opens the periodic table, and picking one switches to the atom tool.
    // Until an element has been picked (ever: it's remembered), the button looks
    // like a small periodic table, so what it opens is obvious.
    auto element = std::make_shared<QString>(QSettings().value("element").toString());
    auto* atom = new QAction(this);
    const IconMaker tableIcon = paintedIcon([](QPainter& p, QColor ink) {
        p.setPen(Qt::NoPen);
        p.setBrush(ink);
        auto cell = [&](int col, int row) { p.drawRect(QRectF(3 + col * 2.6, 6 + row * 2.6, 2, 2)); };
        for (int row = 0; row < 4; ++row) cell(0, row), cell(6, row);  // groups 1 and 18
        for (int row = 1; row < 4; ++row) cell(1, row), cell(4, row), cell(5, row);
        for (int row = 2; row < 4; ++row) cell(2, row), cell(3, row);  // the d block
        for (int col = 1; col < 6; ++col) cell(col, 5);                 // f block underneath
    });
    const IconMaker atomIcon = [element, tableIcon] {
        return element->isEmpty() ? tableIcon() : docIcon(textDoc(*element))();
    };
    if (!element->isEmpty()) canvas_->setElement(chem::atomicNumber(element->toStdString()));
    atom->setIcon(atomIcon());
    icons_.push_back({atom, atomIcon});
    atom->setToolTip(tr("Atom: click to place or relabel (element from the arrow's periodic table; "
                        "or point at an atom and type N, O, S…)"));
    atom->setCheckable(true);
    group->addAction(atom);
    connect(atom, &QAction::toggled, this, [this, atom](bool on) {
        if (on) statusBar()->showMessage(atom->toolTip());
    });
    connect(atom, &QAction::triggered, this, [this] { canvas_->setTool(T::Atom); });
    if (slot % 2) ++slot;
    auto* atomButton = button(atom);
    grid->addWidget(atomButton, (slot - 1) / 2, 0, 1, 2);  // full width
    ++slot;
    atomButton->setPopupMode(QToolButton::MenuButtonPopup);
    atomButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
    auto* menu = new QMenu(atomButton);
    auto* table = new QWidgetAction(menu);
    table->setDefaultWidget(periodicTable([=, this](int z) {
        *element = QString::fromStdString(chem::symbol(z));
        QSettings().setValue("element", *element);
        canvas_->setElement(z);
        canvas_->setTool(T::Atom);
        atom->setChecked(true);
        atom->setIcon(atomIcon());
        menu->close();
    }));
    menu->addAction(table);
    atomButton->setMenu(menu);
    auto charge = [](bool plus) {
        return paintedIcon([plus](QPainter& p, QColor ink) {
            p.setPen(QPen(ink, 1.3));
            p.drawEllipse(QPointF(12, 12), 7, 7);
            p.drawLine(QPointF(8.5, 12), QPointF(15.5, 12));
            if (plus) p.drawLine(QPointF(12, 8.5), QPointF(12, 15.5));
        });
    };
    add(charge(true), tr("Positive charge: click an atom to add +1"), tool(T::ChargePlus));
    add(charge(false), tr("Negative charge: click an atom to add −1"), tool(T::ChargeMinus));
    // Colour tool: paints atoms, bonds, arrows and text; the arrow picks the colour.
    auto colour = std::make_shared<QColor>(canvas_->colour());
    const IconMaker colourIcon = paintedIcon([colour](QPainter& p, QColor ink) {
        p.setPen(QPen(ink, 1));
        p.setBrush(*colour);
        p.drawRoundedRect(QRectF(5, 5, 14, 14), 3, 3);
    });
    auto* colourTool = add(colourIcon, tr("Colour: click an atom, bond, arrow or text to paint it (again to clear); "
                                          "pick the colour from the arrow"),
                           tool(T::Colour));
    for (auto* b : palette->findChildren<QToolButton*>())
        if (b->defaultAction() == colourTool) {
            b->setPopupMode(QToolButton::MenuButtonPopup);
            b->setMenu(colourMenu(b,
                                  {QColor(214, 39, 40), QColor(255, 127, 14), QColor(44, 160, 44), QColor(31, 119, 180),
                                   QColor(148, 103, 189), QColor(227, 119, 194), QColor(127, 127, 127), Qt::black},
                                  [this] { return canvas_->colour(); },
                                  [=, this](QColor c) {
                                      *colour = c;
                                      canvas_->setColour(c);
                                      canvas_->setTool(T::Colour);
                                      colourTool->setChecked(true);
                                      colourTool->setIcon(colourIcon());
                                  }));
        }
    section();
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
        Document blank;
        blank.style = QSettings().value("defaultStyle").toString();  // Edit > Preferences
        canvas_->setDocumentSilently(blank);
        path_.clear();
        updateTitle();
    });
    file->addAction(tr("&Open…"), QKeySequence::Open, this, [this] {
        if (!maybeSave()) return;
        QString p = QFileDialog::getOpenFileName(this, tr("Open"), {},
                                                 tr("Structures (*.penz *.mol *.sdf *.smi *.inchi *.cdxml *.cdx);;Penzene figures (*.svg *.png);;All files (*)"));
        if (!p.isEmpty()) openFile(p);
    });
    auto* recent = file->addMenu(tr("Open &Recent"));
    connect(recent, &QMenu::aboutToShow, this, [this, recent] {
        recent->clear();
        for (const QString& p : recentFiles())
            recent->addAction(QFileInfo(p).fileName(), this, [this, p] {
                if (maybeSave()) openFile(p);
            })->setToolTip(p);
        if (recent->isEmpty()) recent->addAction(tr("No recent files"))->setEnabled(false);
        recent->addSeparator();
        recent->addAction(tr("Clear Menu"), this, [] { QSettings().remove("recentFiles"); });
    });
    file->addAction(tr("&Save"), QKeySequence::Save, this, &MainWindow::save);
    file->addAction(tr("Save &As…"), QKeySequence::SaveAs, this, &MainWindow::saveAs);
    file->addSeparator();
    file->addAction(tr("Import &SMILES…"), QKeySequence(tr("Ctrl+Shift+I")), this, &MainWindow::importSmiles);
    file->addAction(tr("Import &Name from PubChem…"), this, &MainWindow::importName);
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
    edit->addSeparator();
    auto* prefs = edit->addAction(tr("&Preferences…"), QKeySequence::Preferences, this, &MainWindow::showPreferences);
    prefs->setMenuRole(QAction::PreferencesRole);  // the app menu on macOS

    auto* structure = menuBar()->addMenu(tr("&Structure"));
    structure->addAction(tr("Flip &Horizontal"), QKeySequence(tr("Ctrl+Shift+H")), this,
                         [this] { canvas_->flipSelection(true); });
    structure->addAction(tr("Flip &Vertical"), QKeySequence(tr("Ctrl+Shift+V")), this,
                         [this] { canvas_->flipSelection(false); });
    auto* arrange = structure->addMenu(tr("&Align and Distribute"));
    using A = Canvas::Align;
    for (auto [label, edge] : {std::pair{tr("Align &Left"), A::Left}, {tr("Align &Centres"), A::HCentre},
                               {tr("Align &Right"), A::Right}, {tr("Align &Top"), A::Top},
                               {tr("Align &Middles"), A::VCentre}, {tr("Align &Bottom"), A::Bottom}})
        arrange->addAction(label, this, [this, edge] { canvas_->alignSelection(edge); });
    arrange->addSeparator();
    arrange->addAction(tr("Distribute &Horizontally"), this, [this] { canvas_->distributeSelection(true); });
    arrange->addAction(tr("Distribute &Vertically"), this, [this] { canvas_->distributeSelection(false); });
    structure->addSeparator();
    structure->addAction(tr("Add Explicit &Hydrogens"), this, [this] {
        canvas_->commit(chem::addHydrogens(canvas_->document()), tr("Add hydrogens"));
    });
    structure->addAction(tr("Remove Explicit Hydro&gens"), this, [this] {
        canvas_->commit(chem::removeHydrogens(canvas_->document()), tr("Remove hydrogens"));
    });
    structure->addAction(tr("&Name from PubChem"), this, [this] {
        const Document doc = canvas_->selectedSubset();
        if (doc.atoms.empty()) return;
        QString error;
        QApplication::setOverrideCursor(Qt::WaitCursor);
        const QString smiles = QString::fromStdString(chem::toSmiles(doc));
        const QString name = pubchem::fetch(pubchem::smilesToNameUrl(), "IUPACName", &error,
                                            pubchem::smilesToNameForm(smiles));
        QApplication::restoreOverrideCursor();
        if (name.isEmpty()) {
            QMessageBox::warning(this, tr("Name from PubChem"),
                                 tr("No name for %1: %2").arg(smiles, error.isEmpty() ? tr("PubChem has no match.") : error));
            return;
        }
        QApplication::clipboard()->setText(name);
        QMessageBox::information(this, tr("Name from PubChem"), tr("%1\n\n(copied to the clipboard)").arg(name));
    });
    structure->addAction(tr("Chec&k Structure…"), QKeySequence(tr("Ctrl+Alt+K")), this, [this] { checkStructure(); });
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
    structure->addAction(tr("C&olour Selection"), this, [this] { canvas_->colourSelection(); });
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
    // Display options belong to the document (saved, and in exports), so changing one is an edit.
    auto setDisplay = [this](auto change, const QString& what) {
        Document next = canvas_->document();
        change(next);
        if (!(next == canvas_->document())) canvas_->commit(next, what);
    };
    auto* carbons = view->addMenu(tr("&Carbon Labels"));
    auto* carbonGroup = new QActionGroup(carbons);
    using CL = Document::CarbonLabels;
    for (auto [text, mode] : {std::pair{tr("&None (skeletal)"), CL::None}, {tr("&Terminal CH₃"), CL::Terminal},
                              {tr("&All Carbons"), CL::All}}) {
        auto* a = carbons->addAction(text, this, [=] {
            setDisplay([mode](Document& d) { d.carbonLabels = mode; }, tr("Carbon labels"));
        });
        a->setCheckable(true);
        a->setData(int(mode));
        carbonGroup->addAction(a);
    }
    auto* implicitH = view->addAction(tr("Show &Implicit Hydrogens"));
    implicitH->setCheckable(true);
    connect(implicitH, &QAction::triggered, this, [=](bool on) {
        setDisplay([on](Document& d) { d.hideImplicitH = !on; }, tr("Implicit hydrogens"));
    });
    connect(canvas_, &Canvas::documentChanged, this, [this, carbonGroup, implicitH] {
        for (auto* a : carbonGroup->actions()) a->setChecked(a->data().toInt() == int(canvas_->document().carbonLabels));
        QSignalBlocker quiet(implicitH);
        implicitH->setChecked(!canvas_->document().hideImplicitH);
    });
    implicitH->setChecked(true);
    carbonGroup->actions().first()->setChecked(true);
    // Stereo labels belong to the document (saved, and in exports), so toggling is an edit.
    auto* stereo = view->addAction(tr("Show &Stereo Labels"));
    stereo->setCheckable(true);
    connect(stereo, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().showStereo == on) return;
        Document next = canvas_->document();
        next.showStereo = on;
        canvas_->commit(next, on ? tr("Show Stereo Labels") : tr("Hide Stereo Labels"));
    });
    connect(canvas_, &Canvas::documentChanged, stereo, [this, stereo] {
        QSignalBlocker quiet(stereo);
        stereo->setChecked(canvas_->document().showStereo);
    });
    auto* numbers = view->addAction(tr("Atom &Numbers"));
    numbers->setCheckable(true);
    numbers->setStatusTip(tr("Number every atom; ' on an atom sets its reaction map number"));
    connect(numbers, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().showAtomNumbers == on) return;
        Document next = canvas_->document();
        next.showAtomNumbers = on;
        canvas_->commit(next, on ? tr("Show atom numbers") : tr("Hide atom numbers"));
    });
    connect(canvas_, &Canvas::documentChanged, numbers, [this, numbers] {
        QSignalBlocker quiet(numbers);
        numbers->setChecked(canvas_->document().showAtomNumbers);
    });
    auto* circles = view->addAction(tr("&Aromatic Circles"));
    circles->setCheckable(true);
    connect(circles, &QAction::toggled, this, [this](bool on) {
        if (canvas_->document().aromaticCircles == on) return;
        Document next = canvas_->document();
        next.aromaticCircles = on;
        next.aromaticCircleOverrides.clear();
        canvas_->commit(next, on ? tr("Aromatic circles") : tr("Kekulé rings"));
    });
    connect(canvas_, &Canvas::documentChanged, circles, [this, circles] {
        QSignalBlocker quiet(circles);
        circles->setChecked(canvas_->document().aromaticCircles);
    });
    auto* selectedCircles = view->addAction(tr("Circles for Selected &Rings"));
    selectedCircles->setStatusTip(tr("Select every atom in an aromatic ring"));
    connect(selectedCircles, &QAction::triggered, this, [this] {
        Document next = canvas_->document();
        const auto& selected = canvas_->selection();
        for (auto ring : chem::aromaticRings(next)) {
            if (!std::all_of(ring.begin(), ring.end(), [&](int i) { return selected.contains(i); })) continue;
            std::sort(ring.begin(), ring.end());
            auto it = std::find(next.aromaticCircleOverrides.begin(), next.aromaticCircleOverrides.end(), ring);
            if (it == next.aromaticCircleOverrides.end()) next.aromaticCircleOverrides.push_back(ring);
            else next.aromaticCircleOverrides.erase(it);
        }
        if (!(next == canvas_->document())) canvas_->commit(next, tr("Toggle aromatic circles"));
    });
    view->addSeparator();
    auto* panelToggle = profileDock_->toggleViewAction();
    panelToggle->setText(tr("&Properties Panel"));
    panelToggle->setShortcut(QKeySequence(tr("Ctrl+I")));
    view->addAction(panelToggle);
    auto* themeMenu = view->addMenu(tr("&Theme"));
    auto* themeGroup = themeGroup_ = new QActionGroup(this);
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
<tr><td><b>Drag onto an atom</b></td><td>merge (Select tool) &nbsp;•&nbsp; <b>Shift+drag</b> move straight; draw a bond at any angle</td></tr>
<tr><td><b>Ctrl+←↑→↓</b></td><td>duplicate across the next arrow that way (or alongside)</td></tr>
<tr><td><b>Alt+← →</b></td><td>rotate 15° &nbsp;•&nbsp; <b>Alt+drag</b> rotate freely • <b>double-click</b> select fragment, or edit text</td></tr>
</table>)"));
        box.exec();
    });
    help->addAction(tr("&About Penzene"), this, [this] {
        QMessageBox::about(this, tr("About Penzene"),
                           tr("<h3>Penzene %1</h3><p>An open-source chemical structure editor.</p>"
                              "<p>GPL-3.0 • <a href='https://github.com/JamesOBrien2/penzene'>GitHub</a></p>"
                              "<p>Chemistry by RDKit. GUI by Qt.</p>").arg(PENZENE_BUILD));
    });
}
