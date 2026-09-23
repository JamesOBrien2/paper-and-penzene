#include "MainWindow.h"
#include "Canvas.h"
#include "Chem.h"

#include <QActionGroup>
#include <QComboBox>
#include <QMenuBar>
#include <QToolBar>
#include <QUndoStack>

MainWindow::MainWindow() : undo_(new QUndoStack(this)), canvas_(new Canvas(undo_, this)) {
    setCentralWidget(canvas_);
    setWindowTitle("Paper & Penzene");
    resize(1100, 750);
    buildTools();
    buildMenus();
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
    auto* edit = menuBar()->addMenu(tr("&Edit"));
    auto* u = undo_->createUndoAction(this);
    u->setShortcut(QKeySequence::Undo);
    auto* r = undo_->createRedoAction(this);
    r->setShortcut(QKeySequence::Redo);
    edit->addAction(u);
    edit->addAction(r);
    edit->addSeparator();
    edit->addAction(tr("Select &All"), QKeySequence::SelectAll, canvas_, &Canvas::selectAll);

    auto* view = menuBar()->addMenu(tr("&View"));
    view->addAction(tr("Zoom &In"), QKeySequence::ZoomIn, this, [this] { canvas_->zoomBy(1.25); });
    view->addAction(tr("Zoom &Out"), QKeySequence::ZoomOut, this, [this] { canvas_->zoomBy(0.8); });
    view->addAction(tr("&Fit to Window"), QKeySequence(tr("Ctrl+0")), canvas_, &Canvas::fitToDocument);
}
