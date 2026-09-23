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
    QApplication::processEvents();
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

TEST_CASE("hotkeys: chain, labels, groups, bonds") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[1].pos);
    f.key("111");  // hotspot follows each new atom
    CHECK(f.doc().atoms.size() == 5);
    CHECK(f.doc().bonds.size() == 4);

    f.key("O");  // Shift+o: OMe, as an abbreviation
    REQUIRE(f.doc().atoms.size() == 5);
    CHECK(f.doc().atoms[4].z == 8);
    CHECK(f.doc().atoms[4].label == "OMe");
    f.key("F");  // Shift+f: CF3 replaces it
    CHECK(f.doc().atoms[4].z == 6);
    CHECK(f.doc().atoms[4].label == "CF3");

    // Shift+arrow jumps atom to atom back along the chain.
    f.hover(f.doc().atoms[1].pos);
    QPointF dir = f.doc().atoms[0].pos - f.doc().atoms[1].pos;
    QTest::keyClick(f.canvas.viewport(), dir.x() < 0 ? Qt::Key_Left : Qt::Key_Right, Qt::ShiftModifier);
    f.key("n");
    CHECK(f.doc().atoms[0].z == 7);

    // Bond hotkeys: 2 on a bond makes it double (on an atom it sprouts a carbonyl).
    const auto& d = f.doc();
    f.hover((d.atoms[d.bonds[1].a].pos + d.atoms[d.bonds[1].b].pos) / 2);
    f.key("2");
    CHECK(f.doc().bonds[1].order == 2);
    f.key("w");
    CHECK(f.doc().bonds[1].stereo == BondStereo::Wedge);
    int a = f.doc().bonds[1].a;
    f.key("w");
    CHECK(f.doc().bonds[1].b == a);  // flipped
}

static std::string noStereo(std::string smi) {
    std::erase(smi, '@');
    return chem::toSmiles(*chem::fromSmiles(smi));
}

// The worked example from ChemDraw's cheat sheet: from H2N-CH3, "42n152o" builds Ala-Ala.
TEST_CASE("hotkeys: ChemDraw dipeptide example") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[0].pos);
    f.key("n");
    f.hover(f.doc().atoms[1].pos);
    f.key("42n152o");
    CHECK(noStereo(chem::toSmiles(f.doc())) == noStereo("CC(N)C(=O)NC(C)C(=O)O"));
    int wedges = 0, hashes = 0;
    for (auto& b : f.doc().bonds) wedges += b.stereo == BondStereo::Wedge, hashes += b.stereo == BondStereo::Hash;
    CHECK(wedges == 1);
    CHECK(hashes == 1);
}

TEST_CASE("hotkeys: context-dependent sprouts") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 0}, {40, 0});  // propane-ish chain
    REQUIRE(f.doc().atoms.size() >= 3);
    f.hover(f.doc().atoms[1].pos);  // secondary carbon
    f.key("2");
    CHECK(f.doc().atoms.back().z == 8);  // ketone on the hotspot
    CHECK(f.doc().bonds.back().order == 2);

    Fixture g;
    g.canvas.setTool(Canvas::Tool::Ring);
    g.canvas.setRing(6, false);
    g.click({0, 0});
    g.hover(g.doc().atoms[0].pos);  // secondary ring carbon
    g.key("9");                     // gem-dimethyl
    CHECK(g.doc().atoms.size() == 8);
    CHECK(g.doc().neighbors(0).size() == 4);

    Fixture t;  // tertiary ring carbon: "6" adds a C-C bond, then a cyclohexane
    t.canvas.setTool(Canvas::Tool::Ring);
    t.canvas.setRing(6, false);
    t.click({0, 0});
    t.hover(t.doc().atoms[0].pos);
    t.key("1");
    t.hover(t.doc().atoms[0].pos);
    t.key("6");
    CHECK(t.doc().atoms.size() == 6 + 1 + 1 + 5);
    CHECK(t.doc().neighbors(0).size() == 4);

    Fixture h;  // phenyl on an aromatic carbon goes via a C-C bond (biphenyl)
    h.canvas.setTool(Canvas::Tool::Ring);
    h.canvas.setRing(6, true);
    h.click({0, 0});
    h.hover(h.doc().atoms[0].pos);
    h.key("a");
    CHECK(noStereo(chem::toSmiles(h.doc())) == noStereo("c1ccc(-c2ccccc2)cc1"));
}

TEST_CASE("hotspot is sticky and arrows walk atom -> bond -> atom") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Bond);
    f.click({0, 0});
    f.hover(f.doc().atoms[1].pos);
    f.hover({300, 300});  // drift off into empty space
    CHECK(f.canvas.hotspotAtom() == 1);
    f.key("1");
    CHECK(f.doc().atoms.size() == 3);
    CHECK(f.canvas.hotspotAtom() == 2);

    QPointF back = f.doc().atoms[1].pos - f.doc().atoms[2].pos;
    auto toward = [&](QPointF v) {
        return std::abs(v.x()) > std::abs(v.y()) ? (v.x() < 0 ? Qt::Key_Left : Qt::Key_Right)
                                                 : (v.y() < 0 ? Qt::Key_Up : Qt::Key_Down);
    };
    QTest::keyClick(f.canvas.viewport(), toward(back));
    CHECK(f.canvas.hotspotBond() == 1);
    QTest::keyClick(f.canvas.viewport(), toward(back));
    CHECK(f.canvas.hotspotAtom() == 1);
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);
    CHECK(f.canvas.hotspotAtom() == -1);
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
    CHECK(Canvas::applyLabel(d, 1, "NO2"));  // abbreviation: one labelled atom
    CHECK(d.atoms[1].label == "NO2");
    CHECK(d.atoms[1].charge == 1);
    CHECK(d.atoms.size() == 2);
    CHECK(chem::toSmiles(d) == "C[N+](=O)[O-]");  // chemistry sees the full group
    CHECK(Canvas::applyLabel(d, 1, "OH"));
    CHECK(d.atoms[1].z == 8);
    CHECK(d.atoms[1].label.isEmpty());
    CHECK(Canvas::applyLabel(d, 1, "C(=O)Cl"));  // SMILES: drawn out
    CHECK(d.atoms.size() == 4);
    CHECK_FALSE(Canvas::applyLabel(d, 0, "notachem!!"));
}

TEST_CASE("abbreviations: valence, clean, expand") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 5}}};
    d.bonds = {{0, 1}, {1, 2}};
    REQUIRE(Canvas::applyLabel(d, 2, "Boc"));
    CHECK_FALSE(chem::atomInfo(d)[2].valenceError);
    CHECK(chem::toSmiles(d) == "CCC(=O)OC(C)(C)C");
    Document clean = chem::clean2D(d);
    CHECK(clean.atoms.size() == 3);  // still abbreviated
    CHECK(clean.bonds.size() == 2);
    CHECK(clean.atoms[2].label == "Boc");
    f.canvas.setDocumentSilently(d);
    f.canvas.expandAbbreviations();
    CHECK(f.doc().atoms.size() == 9);
    CHECK(f.doc().atoms[2].label.isEmpty());
    CHECK(chem::toSmiles(f.doc()) == "CCC(=O)OC(C)(C)C");
    auto back = Document::fromJson(d.toJson());
    REQUIRE(back);
    CHECK(back->atoms[2].label == "Boc");
}

static double angleAt(const Document& d, int centre, int x, int y) {
    QPointF u = d.atoms[x].pos - d.atoms[centre].pos, v = d.atoms[y].pos - d.atoms[centre].pos;
    return std::acos((u.x() * v.x() + u.y() * v.y()) / (std::hypot(u.x(), u.y()) * std::hypot(v.x(), v.y()))) * 180 / M_PI;
}

// #43: sp centres (allenes, alkynes) are linear.
TEST_CASE("allene and alkyne centres are linear") {
    Fixture f;  // drawing double bonds onto a double bond
    f.canvas.setTool(Canvas::Tool::Bond);
    f.canvas.setBondOrder(2);
    f.click({0, 0});
    f.click(f.doc().atoms[1].pos);
    REQUIRE(f.doc().atoms.size() == 3);
    CHECK(angleAt(f.doc(), 1, 0, 2) > 179);

    Fixture g;  // zig-zag chain, then make both bonds double: terminal atom swings into line
    g.canvas.setTool(Canvas::Tool::Chain);
    g.drag({0, 0}, {40, 0});
    REQUIRE(g.doc().atoms.size() == 4);
    for (int bi : {1, 2}) {
        const auto& d = g.doc();
        g.hover((d.atoms[d.bonds[bi].a].pos + d.atoms[d.bonds[bi].b].pos) / 2);
        g.key("2");
    }
    CHECK(angleAt(g.doc(), 2, 1, 3) > 179);

    Fixture h;  // growing from an alkyne carbon continues straight
    h.canvas.setTool(Canvas::Tool::Bond);
    h.canvas.setBondOrder(3);
    h.click({0, 0});
    h.hover(h.doc().atoms[1].pos);
    h.key("1");
    CHECK(angleAt(h.doc(), 1, 0, 2) > 179);
}

TEST_CASE("arrows: draw, restyle, select, move, delete; text subscripts") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Arrow);
    f.canvas.setArrow(ArrowKind::Reaction, false);
    f.drag({0, 0}, {50, 3});  // snaps to horizontal
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(std::abs(f.doc().arrows[0].to.y()) < 0.5);

    f.canvas.setArrow(ArrowKind::Equilibrium, false);
    f.click({25, 0});  // restyles instead of adding
    REQUIRE(f.doc().arrows.size() == 1);
    CHECK(f.doc().arrows[0].kind == ArrowKind::Equilibrium);

    f.canvas.setArrow(ArrowKind::Reaction, true);
    f.drag({0, 40}, {40, 40});
    REQUIRE(f.doc().arrows.size() == 2);
    const double bend = f.doc().arrows[1].bend;
    CHECK(bend != 0);
    QPointF mid = arrowPath(f.doc().arrows[1]).pointAtPercent(0.5);
    f.click(mid);  // same tool again: flips the curve
    CHECK(f.doc().arrows[1].bend == -bend);

    Document withText = f.doc();
    withText.texts.push_back({{0, -30}, "CH2Cl2"});
    f.canvas.setDocumentSilently(withText);
    f.canvas.setTool(Canvas::Tool::Select);
    f.drag({-10, -50}, {60, 10});  // rubber band: straight arrow and text, not the curve
    CHECK(f.canvas.selectedArrows() == QSet<int>{0});
    CHECK(f.canvas.selectedTexts() == QSet<int>{0});
    f.drag({25, 0}, {25, 20});  // drag the arrow: text moves with it
    CHECK(std::abs(f.doc().arrows[0].from.y() - 20) < 1);
    CHECK(std::abs(f.doc().texts[0].pos.y() + 10) < 1);
    f.canvas.deleteSelection();
    CHECK(f.doc().arrows.size() == 1);
    CHECK(f.doc().texts.empty());

    // Formula subscripts sit below the baseline; digits after a space don't.
    auto bottom = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().bottom(); };
    CHECK(bottom("H2") > bottom("H") + 1);
    CHECK(bottom("80 C") <= bottom("H") + 0.5);
    CHECK(bottom("(2 equiv)") <= bottom("(") + 0.5);

    // Tabs jump to stops, so columns line up; leading spaces indent.
    auto left = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().left(); };
    auto right = [](const QString& s) { return textPath({{0, 0}, s}).boundingRect().right(); };
    CHECK(std::abs(right("\tA") - right("ab\tA")) < 0.01);
    CHECK(left("\tA") > left("A") + 5);
    CHECK(left("  A") > left("A") + 2);
    CHECK(bottom("CH2Cl2") > bottom("CHCl") + 1);  // run layout keeps subscripts
}

TEST_CASE("hotkeys: bond styles, positions, chair; duplicate across an arrow; tool keys") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Chain);
    f.drag({0, 0}, {40, 0});
    REQUIRE(f.doc().bonds.size() >= 2);
    auto hoverBond = [&](int i) {
        const auto& d = f.doc();
        f.hover((d.atoms[d.bonds[i].a].pos + d.atoms[d.bonds[i].b].pos) / 2);
    };
    hoverBond(0);
    f.key("y");
    CHECK(f.doc().bonds[0].stereo == BondStereo::Wavy);
    f.key("B");
    CHECK(f.doc().bonds[0].order == 2);
    CHECK(f.doc().bonds[0].stereo == BondStereo::Bold);
    f.key("r");
    CHECK(f.doc().bonds[0].position == BondPosition::Right);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());

    const size_t before = f.doc().atoms.size();
    hoverBond(1);
    f.key("9");  // chair: four new atoms, bonds all about one bond long
    CHECK(f.doc().atoms.size() == before + 4);
    for (const auto& b : f.doc().bonds) {
        QPointF v = f.doc().atoms[b.a].pos - f.doc().atoms[b.b].pos;
        CHECK(std::abs(std::hypot(v.x(), v.y()) - kBondLength) < 0.15 * kBondLength);
    }

    // Duplicate across an arrow on the right.
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    d.arrows = {{{40, 0}, {80, 0}}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setSelection({0, 1});
    f.canvas.duplicateSelection({1, 0});
    REQUIRE(f.doc().atoms.size() == 4);
    CHECK(f.doc().atoms[2].pos.x() > 80);  // beyond the arrow head
    CHECK(f.canvas.selection() == QSet<int>{2, 3});

    QString picked;
    QObject::connect(&f.canvas, &Canvas::toolKey, [&](const QString& k) { picked = k; });
    QTest::keyClick(f.canvas.viewport(), Qt::Key_Escape);
    f.key("e");
    CHECK(picked == "e");
}
