#include "Canvas.h"
#include "Chem.h"
#include "MainWindow.h"

#include <QApplication>
#include <QTest>
#include <QUndoStack>
#include <catch2/catch_test_macros.hpp>

// Drives the real canvas with synthetic mouse events (QT_QPA_PLATFORM=offscreen).
struct App {  // base class so the QApplication exists before any widget member
    App() {
        static int argc = 1;
        static char name[] = "tests";
        static char* argv[] = {name};
        static QApplication app(argc, argv);
    }
};

struct Fixture : App {
    Fixture() {
        canvas.resize(800, 600);
        canvas.show();
    }
    QPoint at(QPointF scene) { return canvas.mapFromScene(scene); }
    void click(QPointF scene) { QTest::mouseClick(canvas.viewport(), Qt::LeftButton, {}, at(scene)); }
    void drag(QPointF from, QPointF to) {
        QTest::mousePress(canvas.viewport(), Qt::LeftButton, {}, at(from));
        QMouseEvent move(QEvent::MouseMove, at(to), canvas.viewport()->mapToGlobal(at(to)), Qt::NoButton,
                         Qt::LeftButton, {});
        QApplication::sendEvent(canvas.viewport(), &move);
        QTest::mouseRelease(canvas.viewport(), Qt::LeftButton, {}, at(to));
    }
    void hover(QPointF scene) {
        QMouseEvent move(QEvent::MouseMove, at(scene), canvas.viewport()->mapToGlobal(at(scene)), Qt::NoButton,
                         Qt::NoButton, {});
        QApplication::sendEvent(canvas.viewport(), &move);
    }
    void key(const QString& k) { QTest::keyClicks(canvas.viewport(), k); }
    const Document& doc() const { return canvas.document(); }
    int doubles() const {
        int n = 0;
        for (auto& b : canvas.document().bonds) n += b.order == 2;
        return n;
    }
    QUndoStack undo;
    Canvas canvas{&undo};
};

TEST_CASE("benzene then fused ring gives naphthalene") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Ring);
    f.canvas.setRing(6, true);
    f.click({0, 0});
    REQUIRE(f.canvas.document().atoms.size() == 6);
    CHECK(f.canvas.document().bonds.size() == 6);
    CHECK(f.doubles() == 3);

    // Click the midpoint of a bond to fuse.
    const auto& d = f.canvas.document();
    QPointF mid = (d.atoms[d.bonds[1].a].pos + d.atoms[d.bonds[1].b].pos) / 2;
    f.click(mid);
    CHECK(f.canvas.document().atoms.size() == 10);
    CHECK(f.canvas.document().bonds.size() == 11);
    CHECK(f.doubles() == 5);

    f.undo.undo();
    CHECK(f.canvas.document().atoms.size() == 6);
    f.undo.redo();
    CHECK(f.canvas.document().atoms.size() == 10);
}

TEST_CASE("bond click, drag, cycle, chain and erase") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});  // empty: new bond
    REQUIRE(f.canvas.document().bonds.size() == 1);
    f.click(f.canvas.document().atoms[1].pos);  // atom: grow zig-zag
    CHECK(f.canvas.document().bonds.size() == 2);

    auto& d = f.canvas.document();
    f.click((d.atoms[0].pos + d.atoms[1].pos) / 2);  // bond: cycle order
    CHECK(f.canvas.document().bonds[0].order == 2);

    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 100}, {60, 100});
    CHECK(f.canvas.document().bonds.size() == 2 + 5);

    f.canvas.setTool(Canvas::Tool::Erase);
    f.click(f.canvas.document().atoms[0].pos);
    CHECK(f.canvas.document().bonds.size() == 6);
}

TEST_CASE("select and delete, charges") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.canvas.setTool(Canvas::Tool::ChargePlus);
    f.click(f.canvas.document().atoms[0].pos);
    CHECK(f.canvas.document().atoms[0].charge == 1);
    f.canvas.selectAll();
    f.canvas.deleteSelection();
    CHECK(f.canvas.document().atoms.empty());
    f.undo.undo();
    CHECK(f.canvas.document().atoms.size() == 2);
}

TEST_CASE("insert centres the fragment and selects it") {
    Fixture f;
    auto frag = chem::fromSmiles("CCO");
    REQUIRE(frag);
    f.canvas.insert(*frag, "Paste");
    f.canvas.insert(*frag, "Paste");
    CHECK(f.canvas.document().atoms.size() == 6);
    CHECK(f.canvas.document().bonds[2].a == 3);
    CHECK(f.canvas.selection() == QSet<int>{3, 4, 5});
}

TEST_CASE("main window screenshot") {
    App app;
    MainWindow w;
    w.resize(1000, 650);
    w.show();
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));
    if (auto out = qgetenv("PENZENE_SCREENSHOT"); !out.isEmpty()) w.grab().save(out);
}

// #42: a ring on a terminal atom must continue straight on, so the substituent
// bond bisects the ring's outside angle.
TEST_CASE("ring on a terminal atom bisects the outside angle") {
    for (int n : {3, 6}) {
        Fixture f;
        f.canvas.setTool(Canvas::Tool::Bond);
        f.click({0, 0});
        const QPointF stem = f.doc().atoms[0].pos;
        f.canvas.setTool(Canvas::Tool::Ring);
        f.canvas.setRing(n, n == 6);
        f.click(f.doc().atoms[1].pos);
        REQUIRE(f.doc().atoms.size() == size_t(n + 1));
        QPointF p = f.doc().atoms[1].pos, sum;
        for (int nb : f.doc().neighbors(1))
            if (f.doc().atoms[nb].pos != stem) {
                QPointF v = f.doc().atoms[nb].pos - p;
                sum += v / std::hypot(v.x(), v.y());
            }
        QPointF s = stem - p;
        double cosang = (sum.x() * s.x() + sum.y() * s.y()) / (std::hypot(sum.x(), sum.y()) * std::hypot(s.x(), s.y()));
        CHECK(cosang < -0.999);
    }
}

TEST_CASE("hotkeys: chain, groups, bonds, fused rings, arrows") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[1].pos);
    f.key("111");  // hotspot follows each new atom
    CHECK(f.doc().atoms.size() == 5);
    CHECK(f.doc().bonds.size() == 4);

    f.key("O");  // last atom -> OMe
    REQUIRE(f.doc().atoms.size() == 6);
    CHECK(f.doc().atoms[4].z == 8);
    CHECK(f.doc().atoms[5].z == 6);

    f.key("F");  // hotspot still on the O: now CF3 (C + 3 F), replacing O
    CHECK(f.doc().atoms[4].z == 6);
    CHECK(f.doc().atoms.size() == 9);

    // Arrow keys walk the hotspot back along the chain.
    f.hover(f.doc().atoms[1].pos);
    QPointF dir = f.doc().atoms[0].pos - f.doc().atoms[1].pos;
    QTest::keyClick(f.canvas.viewport(), dir.x() < 0 ? Qt::Key_Left : Qt::Key_Right);
    f.key("n");
    CHECK(f.doc().atoms[0].z == 7);

    // Bond hotkeys.
    const auto& d = f.doc();
    QPointF mid = (d.atoms[d.bonds[1].a].pos + d.atoms[d.bonds[1].b].pos) / 2;
    f.hover(mid);
    f.key("2");
    CHECK(f.doc().bonds[1].order == 2);
    f.key("w");
    CHECK(f.doc().bonds[1].stereo == BondStereo::Wedge);
    int a = f.doc().bonds[1].a;
    f.key("w");
    CHECK(f.doc().bonds[1].b == a);  // flipped
}

TEST_CASE("hotkeys: fused ring on a bond") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover((f.doc().atoms[0].pos + f.doc().atoms[1].pos) / 2);
    f.key("a");
    CHECK(f.doc().atoms.size() == 6);  // benzene shares the bond's 2 atoms
    CHECK(f.doubles() == 3);
}

TEST_CASE("applyLabel understands elements, groups and SMILES") {
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    CHECK(Canvas::applyLabel(d, 1, "Br"));
    CHECK(d.atoms[1].z == 35);
    CHECK(Canvas::applyLabel(d, 1, "NO2"));
    CHECK(d.atoms[1].charge == 1);
    CHECK(d.atoms.size() == 4);
    CHECK(chem::toSmiles(d) == "C[N+](=O)[O-]");
    CHECK_FALSE(Canvas::applyLabel(d, 0, "notachem!!"));
}
