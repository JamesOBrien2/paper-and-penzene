#include "MainWindow.h"
#include "Canvas.h"
#include "Chem.h"

#include <QActionGroup>
#include <QApplication>
#include <QClipboard>
#include <QCloseEvent>
#include <QComboBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QStatusBar>
#include <QToolBar>
#include <QUndoStack>

static const char* kMolMime = "chemical/x-mdl-molfile";

MainWindow::MainWindow() : undo_(new QUndoStack(this)), canvas_(new Canvas(undo_, this)) {
    setCentralWidget(canvas_);
    setWindowTitle("Paper & Penzene");
    resize(1100, 750);
    buildTools();
    buildMenus();
    connect(undo_, &QUndoStack::cleanChanged, this, &MainWindow::updateTitle);
    updateTitle();
}

void MainWindow::updateTitle() {
    QString name = path_.isEmpty() ? tr("Untitled") : QFileInfo(path_).fileName();
    setWindowTitle(name + "[*] — Paper & Penzene");
    setWindowModified(!undo_->isClean());
}

bool MainWindow::openFile(const QString& path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        QMessageBox::warning(this, tr("Open"), tr("Cannot read %1").arg(path));
        return false;
    }
    QByteArray data = f.readAll();
    std::optional<Document> doc = path.endsWith(".penz", Qt::CaseInsensitive)
                                      ? Document::fromJson(data)
                                      : chem::fromMolBlock(data.toStdString());
    if (!doc) {
        QMessageBox::warning(this, tr("Open"), tr("%1 is not a structure file I can read.").arg(path));
        return false;
    }
    undo_->clear();
    canvas_->setDocumentSilently(*doc);
    canvas_->fitToDocument();
    path_ = path;
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
    if (doc.atoms.empty()) return;
    auto* mime = new QMimeData;
    mime->setImageData(renderImage(doc));
    mime->setData("image/svg+xml", renderSvg(doc));
    std::string mol = chem::toMolBlock(doc), smi = chem::toSmiles(doc);
    mime->setData(kMolMime, QByteArray::fromStdString(mol));
    mime->setText(QString::fromStdString(smi.empty() ? mol : smi));
    QApplication::clipboard()->setMimeData(mime);
}

void MainWindow::paste() {
    const QMimeData* mime = QApplication::clipboard()->mimeData();
    std::string text = mime->hasFormat(kMolMime) ? mime->data(kMolMime).toStdString()
                                                 : mime->text().trimmed().toStdString();
    if (text.empty()) return;
    auto doc = text.find("M  END") != std::string::npos ? chem::fromMolBlock(text) : chem::fromSmiles(text);
    if (doc) canvas_->insert(*doc, tr("Paste"));
    else statusBar()->showMessage(tr("Clipboard has no structure or SMILES"), 4000);
}

void MainWindow::buildTools() {
    auto* bar = addToolBar(tr("Tools"));
    bar->setObjectName("tools");
    addToolBar(Qt::LeftToolBarArea, bar);
    auto* group = new QActionGroup(this);

    // ponytail: text glyphs for tool icons; draw real icons when someone cares.
    auto add = [&](const QString& label, const QString& tip, auto setup) {
        auto* a = bar->addAction(label);
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
    add("⬚", tr("Select (drag to move, Alt+drag to rotate, double-click for fragment)"), tool(T::Select));
    add("╱", tr("Single bond / chain start"), bond(1))->setChecked(true);
    add("═", tr("Double bond"), bond(2));
    add("≡", tr("Triple bond"), bond(3));
    add("▶", tr("Wedge bond"), tool(T::Wedge));
    add("┇", tr("Hashed bond"), tool(T::Hash));
    add("⦚", tr("Chain"), tool(T::Chain));
    bar->addSeparator();

    auto ring = [this](int n, bool arom) {
        return [this, n, arom] { canvas_->setTool(T::Ring), canvas_->setRing(n, arom); };
    };
    add("⌬", tr("Benzene"), ring(6, true));
    const char* shapes[] = {"△", "□", "⬠", "⬡", "7", "8"};
    for (int n = 3; n <= 8; ++n) add(shapes[n - 3], tr("%1-membered ring").arg(n), ring(n, false));
    bar->addSeparator();

    auto* elements = new QComboBox;
    for (auto s : {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H", "B", "Si"}) elements->addItem(s);
    elements->setToolTip(tr("Element for the atom tool (or hover an atom and press C, N, O…)"));
    auto* atom = add("A", tr("Atom"), tool(T::Atom));
    bar->addWidget(elements);
    connect(elements, &QComboBox::currentTextChanged, this, [this, atom](const QString& s) {
        canvas_->setElement(chem::atomicNumber(s.toStdString()));
        canvas_->setTool(T::Atom);
        atom->setChecked(true);
    });
    add("⊕", tr("Positive charge"), tool(T::ChargePlus));
    add("⊖", tr("Negative charge"), tool(T::ChargeMinus));
    add("⌫", tr("Eraser"), tool(T::Erase));
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
                                                 tr("Structures (*.penz *.mol *.sdf);;All files (*)"));
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
    edit->addAction(tr("&Paste"), QKeySequence::Paste, this, &MainWindow::paste);
    edit->addAction(tr("&Delete"), canvas_, &Canvas::deleteSelection);
    edit->addSeparator();
    edit->addAction(tr("Select &All"), QKeySequence::SelectAll, canvas_, &Canvas::selectAll);

    auto* structure = menuBar()->addMenu(tr("&Structure"));
    // ponytail: cleans the whole document; clean just the selection when someone asks.
    structure->addAction(tr("&Clean Structure"), QKeySequence(tr("Ctrl+Shift+K")), this, [this] {
        canvas_->commit(chem::clean2D(canvas_->document()), tr("Clean"));
    });

    auto* view = menuBar()->addMenu(tr("&View"));
    view->addAction(tr("Zoom &In"), QKeySequence::ZoomIn, this, [this] { canvas_->zoomBy(1.25); });
    view->addAction(tr("Zoom &Out"), QKeySequence::ZoomOut, this, [this] { canvas_->zoomBy(0.8); });
    view->addAction(tr("&Fit to Window"), QKeySequence(tr("Ctrl+0")), canvas_, &Canvas::fitToDocument);

    auto* help = menuBar()->addMenu(tr("&Help"));
    help->addAction(tr("&Keyboard Shortcuts"), QKeySequence(tr("F1")), this, [this] {
        QMessageBox box(this);
        box.setWindowTitle(tr("Keyboard Shortcuts"));
        box.setTextFormat(Qt::RichText);
        box.setText(tr(R"(<p>Hover an atom or bond and press a key.</p>
<table cellspacing="6">
<tr><th colspan="2" align="left">Atom</th></tr>
<tr><td><b>1 2 3</b></td><td>add single / double / triple bond (hotspot follows, so 1111 draws a chain)</td></tr>
<tr><td><b>← ↑ → ↓</b></td><td>move the hotspot to a neighbouring atom</td></tr>
<tr><td><b>4–8</b></td><td>spiro ring of that size</td></tr>
<tr><td><b>c n o s P f l b i h B</b></td><td>C N O S P F Cl Br I H B &nbsp;(w = N, q = O)</td></tr>
<tr><td><b>O N F S</b></td><td>OMe, NO<sub>2</sub>, CF<sub>3</sub>, SiH<sub>3</sub></td></tr>
<tr><td><b>m e p/a K v u</b></td><td>Me, Et, Ph, tBu, cyclopropyl, cyclobutyl</td></tr>
<tr><td><b>E x y t Z</b></td><td>CO<sub>2</sub>Me, Ac, CN, Boc, N<sub>3</sub></td></tr>
<tr><td><b>+ &minus;</b></td><td>charge</td></tr>
<tr><td><b>Enter</b></td><td>type a label: element, group (OMe, Boc…) or SMILES</td></tr>
<tr><td><b>Delete</b></td><td>delete atom (or selection)</td></tr>
<tr><th colspan="2" align="left">Bond</th></tr>
<tr><td><b>1 2 3</b></td><td>bond order</td></tr>
<tr><td><b>w h</b></td><td>wedge / hashed (press again to flip)</td></tr>
<tr><td><b>a v 4–8</b></td><td>fuse benzene / cyclopropane / ring of that size</td></tr>
<tr><th colspan="2" align="left">Selection</th></tr>
<tr><td><b>Alt+← →</b></td><td>rotate 15°</td></tr>
<tr><td><b>Alt+drag</b></td><td>rotate freely &nbsp; <b>double-click</b>: select fragment</td></tr>
</table>)"));
        box.exec();
    });
    help->addAction(tr("&About Paper && Penzene"), this, [this] {
        QMessageBox::about(this, tr("About Paper & Penzene"),
                           tr("<h3>Paper &amp; Penzene %1</h3><p>An open-source chemical structure editor.</p>"
                              "<p>GPL-3.0 • <a href='https://github.com/JamesOBrien2/paper-and-penzene'>GitHub</a></p>"
                              "<p>Chemistry by RDKit. GUI by Qt.</p>").arg(PENZENE_VERSION));
    });
}
