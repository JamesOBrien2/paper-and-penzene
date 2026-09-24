#include "Canvas.h"
#include "Chem.h"
#include "Edit.h"
#include "MainWindow.h"
#include "PubChem.h"
#include "Render.h"

#include <QApplication>
#include <QSettings>
#include <QStatusBar>
#include <QMessageBox>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QListWidget>

#include <QLabel>
#include <QTest>

#include <QComboBox>
#include <QDialog>
#include <QSpinBox>
#include <QTest>
#include <QTemporaryDir>
#include <QLineEdit>
#include <QTcpServer>
#include <QTcpSocket>
#include <QTimer>
#include <QClipboard>
#include <QMimeData>
#include <qpa/qwindowsysteminterface.h>
#include <QMenu>
#include <QToolButton>
#include <QWidgetAction>
#include <QUndoStack>
#include <catch2/catch_test_macros.hpp>

// Drives the real canvas with synthetic mouse events (QT_QPA_PLATFORM=offscreen).
struct App {  // base class so the QApplication exists before any widget member
    App() {
        static int argc = 1;
        static char name[] = "tests";
        static char* argv[] = {name};
        static QApplication app(argc, argv);
        QApplication::setOrganizationName("penzene-tests");  // keep the user's settings out of it
        QStandardPaths::setTestModeEnabled(true);  // and their app data (autosave)
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
    QSettings().setValue("theme", qEnvironmentVariable("PENZENE_THEME", "Light"));
    QSettings().remove("element");  // the periodic-table icon, as on first run
    MainWindow w;
    w.resize(1000, 800);
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
    CHECK(edit::applyLabel(d, 1, "Br"));
    CHECK(d.atoms[1].z == 35);
    CHECK(edit::applyLabel(d, 1, "NO2"));  // abbreviation: one labelled atom
    CHECK(d.atoms[1].label == "NO2");
    CHECK(d.atoms[1].charge == 1);
    CHECK(d.atoms.size() == 2);
    CHECK(chem::toSmiles(d) == "C[N+](=O)[O-]");  // chemistry sees the full group
    CHECK(edit::applyLabel(d, 1, "OH"));
    CHECK(d.atoms[1].z == 8);
    CHECK(d.atoms[1].label.isEmpty());
    CHECK(edit::applyLabel(d, 1, "C(=O)Cl"));  // SMILES: drawn out
    CHECK(d.atoms.size() == 4);
    CHECK_FALSE(edit::applyLabel(d, 0, "notachem!!"));
}

TEST_CASE("abbreviations: valence, clean, expand") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 5}}};
    d.bonds = {{0, 1}, {1, 2}};
    REQUIRE(edit::applyLabel(d, 2, "Boc"));
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

TEST_CASE("drawing style presets: JDP from its ChemDraw stationery") {
    App app;
    const auto& jdp = drawingStyle("JDP");
    CHECK(jdp.name == "JDP");
    CHECK(std::abs(jdp.lineWidth * jdp.bondLength / kBondLength - 0.879) < 1e-9);  // native pt in exports
    CHECK(drawingStyle("").name == "ACS 1996");
    CHECK(drawingStyle("no such style").name == "ACS 1996");

    // The style is part of the document: it changes the rendering.
    Document d = *chem::fromSmiles("CC(=O)O");
    Document j = d;
    j.style = "JDP";
    CHECK(renderSvg(d) != renderSvg(j));
}

TEST_CASE("themes: Catppuccin palettes; exports stay black") {
    CHECK(theme("Catppuccin Mocha").paper == QColor("#1e1e2e"));
    CHECK(theme("Catppuccin Latte").ink == QColor("#4c4f69"));
    CHECK(theme("no such theme").name == "System");
    App app;
    Fixture f;
    f.canvas.setTheme(theme("Catppuccin Mocha"));
    auto doc = chem::fromSmiles("CO");
    REQUIRE(doc);
    QImage img = renderImage(*doc, 72);
    bool dark = false;
    for (int y = 0; y < img.height(); ++y)
        for (int x = 0; x < img.width(); ++x)
            if (QColor c = img.pixelColor(x, y); c.alpha() > 200 && c.lightness() < 60) dark = true;
    CHECK(dark);  // black ink, not the theme's pale text
}

TEST_CASE("ring fill: click inside toggles; survives delete, copy and save") {
    Fixture f;
    auto nap = chem::fromSmiles("c1ccc2ccccc2c1");
    REQUIRE(nap);
    f.canvas.setDocumentSilently(*nap);
    // Centre of one ring: mean of the atoms in its fill.
    auto rings = chem::rings(f.doc());
    REQUIRE(rings.size() == 2);
    QPointF c;
    for (int i : rings[0]) c += f.doc().atoms[i].pos / 6;
    f.canvas.setTool(Canvas::Tool::Fill);
    f.click(c);
    REQUIRE(f.doc().fills.size() == 1);
    CHECK(f.doc().fills[0].atoms.size() == 6);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());
    f.click(c);  // same colour again: cleared
    CHECK(f.doc().fills.empty());
    f.click(c);
    REQUIRE(f.doc().fills.size() == 1);
    Document d = f.doc();
    d.removeAtoms({d.fills[0].atoms[0]});
    CHECK(d.fills.empty());  // a ring missing an atom loses its fill
}

TEST_CASE("2 on a double bond swaps the side of its second line") {
    Fixture f;
    auto d = chem::fromSmiles("CC=CC");  // trans-2-butene: an offset (not centred) double bond
    REQUIRE(d);
    f.canvas.setDocumentSilently(*d);
    int db = -1;
    for (int i = 0; i < int(f.doc().bonds.size()); ++i)
        if (f.doc().bonds[i].order == 2) db = i;
    REQUIRE(db >= 0);
    const auto& doc = f.doc();
    f.hover((doc.atoms[doc.bonds[db].a].pos + doc.atoms[doc.bonds[db].b].pos) / 2);
    f.key("2");
    CHECK(f.doc().bonds[db].order == 2);
    const BondPosition first = f.doc().bonds[db].position;
    CHECK(first != BondPosition::Auto);
    f.key("2");
    CHECK(f.doc().bonds[db].position != first);  // and back again
    CHECK(f.doc().bonds[db].position != BondPosition::Auto);
}

TEST_CASE("picking from the periodic table switches to the atom tool") {
    App app;
    QSettings().remove("element");  // first run: nothing picked yet
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(canvas);
    QToolButton* nitrogen = nullptr;  // in the periodic table popup
    for (auto* wa : w.findChildren<QWidgetAction*>())
        if (wa->defaultWidget())
            for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
                if (b->text() == "N") nitrogen = b;
    REQUIRE(nitrogen);
    nitrogen->click();  // no need to pick the atom tool first
    QTest::mouseClick(canvas->viewport(), Qt::LeftButton, {}, canvas->viewport()->rect().center());
    REQUIRE(canvas->document().atoms.size() == 1);
    CHECK(canvas->document().atoms[0].z == 7);
    CHECK(QSettings().value("element").toString() == "N");  // remembered for next time (#128)
    MainWindow again;  // a new window starts with the picked element
    auto* c2 = again.findChild<Canvas*>();
    c2->setTool(Canvas::Tool::Atom);
    again.show();
    QTest::mouseClick(c2->viewport(), Qt::LeftButton, {}, c2->viewport()->rect().center());
    REQUIRE(c2->document().atoms.size() == 1);
    CHECK(c2->document().atoms[0].z == 7);
}

TEST_CASE("Ac, Pr and Ts are groups, not actinium, praseodymium and tennessine (#125)") {
    for (auto [label, formula] : {std::pair{"Ac", "C3H6O"}, {"Pr", "C4H10"}, {"Ts", "C8H10O2S"}}) {
        Document d;
        d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
        d.bonds = {{0, 1}};
        REQUIRE(edit::applyLabel(d, 1, label));
        CHECK(d.atoms[1].label == label);
        CHECK(chem::properties(d)->formula == formula);
    }
    Fixture f;  // and the Shift+A hotkey
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);
    f.hover(f.doc().atoms[1].pos);
    f.key("A");
    CHECK(f.doc().atoms[1].label == "Ac");
}

TEST_CASE("curved arrows are circular arcs, exact past 180 degrees (#86)") {
    // Chord 40, bend 30: more than a semicircle, radius (20² + 30²) / 60 = 21.7.
    const QRectF r = arrowPath({{0, 0}, {40, 0}, ArrowKind::Reaction, 30}).boundingRect();
    CHECK(std::abs(r.top() + 30) < 0.2);          // the arc's midpoint sits 30 above the chord
    CHECK(r.left() < -1.4);                          // and it bulges past both ends
    CHECK(r.right() > 41.4);
    // A small bend still matches the old midpoint.
    QPointF mid = arrowPath({{0, 0}, {40, 0}, ArrowKind::Reaction, 8}).pointAtPercent(0.5);
    CHECK(std::abs(mid.y() + 8) < 0.2);
}

TEST_CASE("flip mirrors (enantiomer with wedges kept), align and distribute (#87)") {
    Fixture f;
    auto ala = chem::fromSmiles("C[C@H](N)C(=O)O");
    REQUIRE(ala);
    Document d = *ala;
    d.arrows.push_back({{60, 0}, {100, 0}, ArrowKind::Reaction, 10});
    int dbl = -1;
    for (int i = 0; i < int(d.bonds.size()); ++i)
        if (d.bonds[i].order == 2) dbl = i;
    d.bonds[dbl].position = BondPosition::Left;
    f.canvas.setDocumentSilently(d);
    const std::string before = chem::toSmiles(f.doc());
    f.canvas.flipSelection(true);
    CHECK(chem::toSmiles(f.doc()) == chem::toSmiles(*chem::fromSmiles("C[C@@H](N)C(=O)O")));  // mirror image
    CHECK(f.doc().bonds[dbl].position == BondPosition::Right);
    CHECK(f.doc().arrows[0].bend == -10);
    f.canvas.flipSelection(true);  // flipping back restores everything
    CHECK(chem::toSmiles(f.doc()) == before);
    for (size_t i = 0; i < d.atoms.size(); ++i)
        CHECK(std::hypot(f.doc().atoms[i].pos.x() - d.atoms[i].pos.x(), f.doc().atoms[i].pos.y() - d.atoms[i].pos.y()) < 1e-6);

    // Three methanols at uneven spacing and heights.
    Document three;
    for (double x : {0.0, 30.0, 100.0}) {
        int c = three.addAtom({x, x / 5});
        int o = three.addAtom({x + kBondLength, x / 5}, 8);
        three.bonds.push_back({c, o});
    }
    f.canvas.setDocumentSilently(three);
    f.canvas.alignSelection(Canvas::Align::Top);
    for (int i = 0; i < 6; ++i) CHECK(std::abs(f.doc().atoms[i].pos.y() - f.doc().atoms[0].pos.y()) < 1e-6);
    f.canvas.distributeSelection(true);
    const double gap1 = f.doc().atoms[2].pos.x() - f.doc().atoms[1].pos.x();
    const double gap2 = f.doc().atoms[4].pos.x() - f.doc().atoms[3].pos.x();
    CHECK(std::abs(gap1 - gap2) < 1e-6);
    CHECK(std::abs(f.doc().atoms[0].pos.x()) < 1e-9);  // outermost objects stay put
    CHECK(std::abs(f.doc().atoms[4].pos.x() - 100) < 1e-9);
}

TEST_CASE("exports come out at the drawing style's own bond length (#92)") {
    App app;
    auto d = chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    REQUIRE(d);
    const QImage acs = renderImage(*d, 72);
    d->style = "RSC";
    const QImage rsc = renderImage(*d, 72);
    CHECK(std::abs(exportScale(*d) - 12.2 / 14.4) < 1e-12);
    // Same drawing, same model bounds (RSC's labels are smaller, so compare against its own bounds).
    const QRectF model = documentBounds(*d);
    CHECK(std::abs(rsc.width() - model.width() * 12.2 / 14.4) <= 1);
    CHECK(rsc.width() < acs.width());
}

// Finds a menu entry by its text, looking inside submenus.
static QAction* findAction(QMenu* menu, const QString& text) {
    for (QAction* a : menu->actions()) {
        if (a->text() == text) return a;
        if (a->menu())
            if (QAction* sub = findAction(a->menu(), text)) return sub;
    }
    return nullptr;
}

TEST_CASE("terminal bond deletion drops the end atom through each UI path") {
    for (int path = 0; path < 3; ++path) {
        Fixture f;
        Document d;
        d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{2 * kBondLength, 0}}};
        d.bonds = {{0, 1}, {1, 2}};
        f.canvas.setDocumentSilently(d);
        const QPointF mid(1.5 * kBondLength, 0);
        if (path == 0) {
            f.hover(mid);
            QTest::keyClick(f.canvas.viewport(), Qt::Key_Delete);
        } else if (path == 1) {
            f.canvas.setTool(Canvas::Tool::Erase);
            f.click(mid);
        } else {
            QMenu* menu = f.canvas.contextMenuAt(mid);
            QAction* action = findAction(menu, "Delete Bond");
            REQUIRE(action);
            action->trigger();
        }
        REQUIRE(f.doc().atoms.size() == 2);
        CHECK(f.doc().bonds.size() == 1);
        CHECK(f.doc().atoms[1].pos == QPointF(kBondLength, 0));
        if (path == 0)
            if (auto out = qgetenv("PENZENE_DELETE_SHOT"); !out.isEmpty()) f.canvas.grab().save(out);
        f.undo.undo();
        CHECK(f.doc() == d);
    }
}

TEST_CASE("right-click menus for atoms, bonds, selection and canvas (#89)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    d.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(d);

    QMenu* atomMenu = f.canvas.contextMenuAt(f.doc().atoms[1].pos);
    REQUIRE(findAction(atomMenu, "N"));
    findAction(atomMenu, "N")->trigger();
    CHECK(f.doc().atoms[1].z == 7);
    findAction(f.canvas.contextMenuAt(f.doc().atoms[1].pos), "Boc")->trigger();
    CHECK(f.doc().atoms[1].label == "Boc");

    QMenu* bondMenu = f.canvas.contextMenuAt({kBondLength / 2, 0});
    REQUIRE(findAction(bondMenu, "Double"));
    findAction(bondMenu, "Double")->trigger();
    CHECK(f.doc().bonds[0].order == 2);
    CHECK(findAction(f.canvas.contextMenuAt({kBondLength / 2, 0}), "Right"));  // double: position submenu

    f.canvas.selectAll();
    CHECK(findAction(f.canvas.contextMenuAt(f.doc().atoms[0].pos), "Flip Horizontal"));
    f.canvas.setSelection({});
    CHECK(findAction(f.canvas.contextMenuAt({200, 200}), "Select All"));
}

TEST_CASE("colour atoms, bonds, arrows and text; exports keep the colour (#82)") {
    Fixture f;
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}, 8}};
    d.bonds = {{0, 1}};
    d.arrows = {{{40, 0}, {80, 0}}};
    d.texts = {{{40, -10}, "heat"}};
    f.canvas.setDocumentSilently(d);
    const QColor red(214, 39, 40);
    f.canvas.setColour(red);
    f.canvas.setTool(Canvas::Tool::Colour);
    f.click(f.doc().atoms[1].pos);
    CHECK(f.doc().atoms[1].color == red);
    f.click(f.doc().atoms[1].pos);  // same colour again: cleared
    CHECK_FALSE(f.doc().atoms[1].color.isValid());

    f.canvas.selectAll();
    f.canvas.colourSelection();
    CHECK(f.doc().bonds[0].color == red);
    CHECK(f.doc().arrows[0].color == red);
    CHECK(f.doc().texts[0].color == red);
    auto back = Document::fromJson(f.doc().toJson());
    REQUIRE(back);
    CHECK(*back == f.doc());

    Document bond;  // a red bond exports red (colours are the user's, unlike theme ink)
    bond.atoms = {{{0, 0}}, {{kBondLength * 3, 0}}};
    bond.bonds = {{0, 1}};
    bond.bonds[0].color = red;
    const QImage img = renderImage(bond, 150);
    bool sawRed = false;
    for (int y = 0; y < img.height(); ++y)
        for (int x = 0; x < img.width(); ++x)
            if (QColor c = img.pixelColor(x, y); c.alpha() > 200 && c.red() > 150 && c.green() < 90) sawRed = true;
    CHECK(sawRed);
}

TEST_CASE("choosing a tool explains it in the status bar (#93)") {
    App app;
    MainWindow w;
    w.show();
    QAction* chain = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->toolTip().startsWith("Chain")) chain = a;
    REQUIRE(chain);
    chain->trigger();
    CHECK(w.statusBar()->currentMessage().startsWith("Chain"));
    CHECK(w.statusBar()->currentMessage().contains("drag"));
}

TEST_CASE("drop an atom on another to merge; Shift for free angles and straight moves (#88)") {
    Fixture f;
    auto drag = [&](QPointF from, QPointF to, Qt::KeyboardModifiers mods) {
        QTest::mousePress(f.canvas.viewport(), Qt::LeftButton, {}, f.at(from));
        QMouseEvent move(QEvent::MouseMove, f.at(to), f.canvas.viewport()->mapToGlobal(f.at(to)), Qt::NoButton,
                         Qt::LeftButton, mods);
        QApplication::sendEvent(f.canvas.viewport(), &move);
        QTest::mouseRelease(f.canvas.viewport(), Qt::LeftButton, mods, f.at(to));
    };
    // Two separate bonds; drag the second's end atom onto the first's.
    Document d;
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{60, 30}}, {{60 + kBondLength, 30}}};
    d.bonds = {{0, 1}, {2, 3}};
    f.canvas.setDocumentSilently(d);
    f.canvas.setTool(Canvas::Tool::Select);
    drag({60, 30}, {kBondLength, 0}, {});
    CHECK(f.doc().atoms.size() == 3);
    CHECK(f.doc().bonds.size() == 2);
    CHECK(f.doc().neighbors(1).size() == 2);  // one connected chain now
    f.undo.undo();
    CHECK(f.doc().atoms.size() == 4);

    // Shift-drag: moves along one axis only.
    f.canvas.setSelection({2, 3});
    drag({60, 30}, {90, 36}, Qt::ShiftModifier);
    CHECK(std::abs(f.doc().atoms[2].pos.y() - 30) < 1e-6);
    CHECK(std::abs(f.doc().atoms[2].pos.x() - 90) < 0.5);

    // Shift while drawing a bond: any angle, not snapped to 30 degrees.
    f.canvas.setDocumentSilently({});
    f.canvas.setTool(Canvas::Tool::Bond);
    drag({0, 0}, {20, 20}, Qt::ShiftModifier);
    REQUIRE(f.doc().bonds.size() == 1);
    QPointF v = f.doc().atoms[1].pos - f.doc().atoms[0].pos;
    CHECK(std::abs(std::atan2(v.y(), v.x()) * 180 / M_PI - 45) < 2);
}

TEST_CASE("drag the ring tool to size a ring (#127)") {
    Fixture f;
    f.canvas.setTool(Canvas::Tool::Ring);
    f.canvas.setRing(6, false);
    // One atom per half bond length of drag, from 3.
    Document none;
    f.canvas.setDocumentSilently(none);
    f.drag({0, 0}, {5.2 * kBondLength / 2, 0});  // 5 steps: an 8-membered ring
    CHECK(f.doc().atoms.size() == 8);
    CHECK(f.doc().bonds.size() == 8);
    f.canvas.setDocumentSilently(none);
    f.drag({0, 0}, {0.5, 0});  // barely moved: a click, so the chosen 6-ring
    CHECK(f.doc().atoms.size() == 6);
    // Dragging from a bond fuses the sized ring onto it.
    Document bond;
    bond.atoms = {{{0, 0}}, {{kBondLength, 0}}};
    bond.bonds = {{0, 1}};
    f.canvas.setDocumentSilently(bond);
    f.drag({kBondLength / 2, 0}, {kBondLength / 2, 2.2 * kBondLength / 2});  // 2 steps: a 5-ring
    CHECK(f.doc().atoms.size() == 5);  // 2 shared + 3 new
    CHECK(f.doc().bonds.size() == 5);
}

TEST_CASE("ring fill colour is picked from the fill button (#126)") {
    App app;
    MainWindow w;
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    QToolButton* fillButton = nullptr;
    for (auto* b : w.findChildren<QToolButton*>())
        if (b->defaultAction() && b->defaultAction()->toolTip().startsWith("Ring fill")) fillButton = b;
    REQUIRE(fillButton);
    REQUIRE(fillButton->menu());
    QToolButton* pink = nullptr;  // a swatch in the drop-down
    for (auto* wa : fillButton->menu()->findChildren<QWidgetAction*>())
        for (auto* b : wa->defaultWidget()->findChildren<QToolButton*>())
            if (b->toolTip() == QColor(255, 214, 214).name()) pink = b;
    REQUIRE(pink);
    pink->click();
    CHECK(canvas->fillColor() == QColor(255, 214, 214));
    CHECK(fillButton->defaultAction()->isChecked());  // and the fill tool is chosen
}

TEST_CASE("recent files, autosave and crash recovery (#91)") {
    App app;
    QSettings().remove("recentFiles");
    QFile::remove(MainWindow::autosavePath());
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    REQUIRE(w.openFile(QString(PENZENE_TEST_DATA) + "/aspirin.mol"));
    CHECK(w.recentFiles().value(0).endsWith("aspirin.mol"));

    canvas->commit(*chem::fromSmiles("CCO"), "edit");  // unsaved changes
    w.autosave();
    REQUIRE(QFile::exists(MainWindow::autosavePath()));

    // A fresh window after a "crash" offers the autosave back; answer Yes.
    MainWindow after;
    QTimer::singleShot(0, [] {
        if (auto* box = qobject_cast<QMessageBox*>(QApplication::activeModalWidget()))
            box->button(QMessageBox::Yes)->click();
    });
    after.offerRecovery();
    auto* c2 = after.findChild<Canvas*>();
    CHECK(chem::toSmiles(c2->document()) == "CCO");
    CHECK_FALSE(QFile::exists(MainWindow::autosavePath()));  // offered once only

    after.autosave();  // once changes are saved (the stack is clean), autosave removes its copy
    REQUIRE(QFile::exists(MainWindow::autosavePath()));
    after.findChild<QUndoStack*>()->setClean();
    after.autosave();
    CHECK_FALSE(QFile::exists(MainWindow::autosavePath()));
}

TEST_CASE("preferences: default style for new documents, export resolution and background (#90)") {
    App app;
    QSettings().remove("defaultStyle");
    QSettings().remove("exportBackground");
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    QTimer::singleShot(0, [] {
        auto* dialog = qobject_cast<QDialog*>(QApplication::activeModalWidget());
        REQUIRE(dialog);
        auto boxes = dialog->findChildren<QComboBox*>();  // theme, style, background
        REQUIRE(boxes.size() == 3);
        boxes[1]->setCurrentText("RSC");
        boxes[2]->setCurrentIndex(1);  // white
        dialog->findChild<QSpinBox*>()->setValue(150);
        dialog->accept();
    });
    w.showPreferences();
    CHECK(QSettings().value("defaultStyle").toString() == "RSC");
    CHECK(QSettings().value("exportDpi").toInt() == 150);
    CHECK(QSettings().value("exportBackground").toString() == "white");

    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&New") a->trigger();
    CHECK(canvas->document().style == "RSC");

    // White background: the corner pixel of an export is opaque white, not clear.
    QTemporaryDir dir;
    auto doc = chem::fromSmiles("CCO");
    REQUIRE(exportDocument(*doc, dir.filePath("x.png"), 150, Qt::white));
    QImage img(dir.filePath("x.png"));
    CHECK(img.pixelColor(0, 0) == QColor(Qt::white));
    REQUIRE(exportDocument(*doc, dir.filePath("y.png")));
    CHECK(QImage(dir.filePath("y.png")).pixelColor(0, 0).alpha() == 0);
    QSettings().remove("defaultStyle");
    QSettings().remove("exportBackground");
    QSettings().remove("exportDpi");
}

TEST_CASE("Check Structure dialog selects the problem's atoms (#95)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CCC(C)O"));
    QWidget* dialog = w.checkStructure();
    auto* list = dialog->findChild<QListWidget*>();
    REQUIRE(list);
    REQUIRE(list->count() == 1);
    list->setCurrentRow(0);
    CHECK(canvas->selection() == QSet<int>{2});  // the unassigned stereocentre
    if (auto out = qgetenv("PENZENE_CHECK_SHOT"); !out.isEmpty()) dialog->grab().save(out);
    dialog->close();
}

TEST_CASE("selected aromatic rings can use circles independently (#98)") {
    App app;
    MainWindow w;
    auto* canvas = w.findChild<Canvas*>();
    auto doc = chem::fromSmiles("c1ccc2ccccc2c1");
    REQUIRE(doc);
    const auto before = chem::toSmiles(*doc);
    auto rings = chem::aromaticRings(*doc);
    REQUIRE(rings.size() == 2);
    canvas->setDocumentSilently(*doc);
    QSet<int> selected(rings[0].begin(), rings[0].end());
    canvas->setSelection(selected);
    QAction* oneRing = nullptr;
    QAction* allRings = nullptr;
    for (auto* action : w.findChildren<QAction*>()) {
        if (action->text() == "Circles for Selected &Rings") oneRing = action;
        if (action->text() == "&Aromatic Circles") allRings = action;
    }
    REQUIRE(oneRing);
    REQUIRE(allRings);
    oneRing->trigger();
    CHECK_FALSE(canvas->document().aromaticCircles);
    CHECK(canvas->document().aromaticCircleOverrides.size() == 1);
    CHECK(chem::toSmiles(canvas->document()) == before);
    auto back = Document::fromJson(canvas->document().toJson());
    REQUIRE(back);
    CHECK(back->aromaticCircleOverrides == canvas->document().aromaticCircleOverrides);
    if (auto out = qgetenv("PENZENE_ONE_CIRCLE_SHOT"); !out.isEmpty())
        CHECK(exportDocument(canvas->document(), QString::fromUtf8(out), 150, Qt::white));
    canvas->setSelection(selected);
    oneRing->trigger();
    CHECK(canvas->document().aromaticCircleOverrides.empty());
    allRings->trigger();
    CHECK(canvas->document().aromaticCircles);
    CHECK(canvas->document().aromaticCircleOverrides.empty());
}

TEST_CASE("properties panel shows descriptors for the selection (#96)") {
    App app;
    MainWindow w;
    w.resize(1100, 700);
    w.show();
    auto* canvas = w.findChild<Canvas*>();
    canvas->setDocumentSilently(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
    QAction* toggle = nullptr;
    for (auto* a : w.findChildren<QAction*>())
        if (a->text() == "&Properties Panel") toggle = a;
    REQUIRE(toggle);
    toggle->trigger();
    QApplication::processEvents();
    QLabel* panel = nullptr;
    for (auto* l : w.findChildren<QLabel*>())
        if (l->text().contains("cLogP")) panel = l;
    REQUIRE(panel);
    CHECK(panel->text().contains("C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>"));
    CHECK(panel->text().contains("63.6"));
    if (auto out = qgetenv("PENZENE_PANEL_SHOT"); !out.isEmpty()) w.grab().save(out);
}

TEST_CASE("PubChem name lookup: URL and response parsing (#28)") {
    CHECK(pubchem::nameToSmilesUrl(" acetylsalicylic acid ").toString(QUrl::FullyEncoded) ==
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/acetylsalicylic%20acid/property/SMILES/JSON");
    const QByteArray found = R"({"PropertyTable": {"Properties": [{"CID": 2244, "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O"}]}})";
    CHECK(pubchem::property(found, "SMILES") == "CC(=O)OC1=CC=CC=C1C(=O)O");
    CHECK(pubchem::property(found, "IUPACName").isEmpty());
    CHECK(pubchem::property(R"({"Fault": {"Code": "PUGREST.NotFound"}})", "SMILES").isEmpty());
    CHECK(pubchem::property("not json", "SMILES").isEmpty());
}

TEST_CASE("PubChem structure lookup: POSTed SMILES and the IUPAC name (#27)") {
    CHECK(pubchem::smilesToNameUrl().toString() ==
          "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/IUPACName/JSON");
    CHECK(pubchem::smilesToNameForm("C#N/C=C/O") == "smiles=C%23N%2FC%3DC%2FO");
    CHECK(pubchem::property(R"({"PropertyTable": {"Properties": [{"CID": 2244, "IUPACName": "2-acetyloxybenzoic acid"}]}})",
                            "IUPACName") == "2-acetyloxybenzoic acid");
}

TEST_CASE("x and r label an atom X and R; free-text labels keep unspecified chemistry (#156)") {
    Document d = *chem::fromSmiles("CCO");
    CHECK(edit::hotkey(d, {0, -1}, "x").atom == 0);
    CHECK(edit::hotkey(d, {2, -1}, "r").atom == 2);
    CHECK(d.atoms[0].label == "X");
    CHECK(d.atoms[2].label == "R");
    CHECK(d.atoms[0].z == 0);
    CHECK(d.atoms.size() == 3);  // nothing invented
    CHECK(chem::toSmiles(d) == "*C*");

    CHECK_FALSE(edit::applyLabel(d, 1, "MgEt"));  // strict unless asked (the Python API)
    REQUIRE(edit::applyLabel(d, 1, "MgEt", true));
    CHECK(d.atoms[1].label == "MgEt");
    auto back = Document::fromJson(d.toJson());
    REQUIRE(back);
    CHECK(back->atoms[1].label == "MgEt");
    CHECK(edit::applyLabel(d, 1, "OMe", true));  // real groups still win
    CHECK(d.atoms[1].z == 8);

    // With no hotspot, x is still the bond tool; on a bond, r still places the double bond.
    Document e = *chem::fromSmiles("C=C");
    CHECK(edit::hotkey(e, {-1, 0}, "r").valid());
    CHECK(e.atoms[0].label.isEmpty());
}

TEST_CASE("the window title names the build (#165)") {
    App app;
    MainWindow w;
    CHECK(w.windowTitle().endsWith("Penzene " PENZENE_BUILD));
    CHECK(QString(PENZENE_BUILD).startsWith(PENZENE_VERSION));
}

TEST_CASE("a PubChem lookup holds user input until it returns (#171)") {
    App app;
    QTcpServer server;  // accepts, stays silent, then hangs up
    REQUIRE(server.listen(QHostAddress::LocalHost));
    QObject::connect(&server, &QTcpServer::newConnection, [&] {
        QTcpSocket* s = server.nextPendingConnection();
        QTimer::singleShot(300, s, [s] { s->close(); });
    });
    QLineEdit typing;
    typing.show();
    QTimer::singleShot(50, [&] {
        // As the window system delivers it (a posted QKeyEvent would bypass the filter).
        QWindowSystemInterface::handleKeyEvent(typing.windowHandle(), QEvent::KeyPress, Qt::Key_A, Qt::NoModifier, "a");
    });
    QString error;
    const QUrl url(QString("http://127.0.0.1:%1/").arg(server.serverPort()));
    CHECK(pubchem::fetch(url, "SMILES", &error).isEmpty());
    CHECK(typing.text().isEmpty());  // not handled mid-request
    QApplication::processEvents();
    CHECK(typing.text() == "a");  // delivered afterwards
}

TEST_CASE("exported SVG and PNG reopen as the editable drawing (#100)") {
    App app;
    Document doc = *chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    doc.atoms[0].color = Qt::red;
    doc.atoms[3].map = 4;
    doc.texts.push_back({{0, 60}, "aspirin"});
    QTemporaryDir dir;
    for (const char* ext : {"svg", "png"}) {
        const QString path = dir.filePath(QString("aspirin.") + ext);
        REQUIRE(exportDocument(doc, path));
        auto back = chem::readFile(path);
        REQUIRE(back);
        CHECK(*back == doc);
    }
    CHECK_FALSE(Document::fromEmbedded(QByteArray("<svg xmlns='http://www.w3.org/2000/svg'/>")));
    CHECK_FALSE(Document::fromEmbedded(QByteArray()));

    // A copied figure pastes back as structure, not as a picture.
    MainWindow w;
    auto* mime = new QMimeData;
    mime->setData("image/svg+xml", renderSvg(doc));
    QApplication::clipboard()->setMimeData(mime);
    w.findChild<Canvas*>()->setDocumentSilently(Document{});
    for (auto* a : w.findChildren<QAction*>())
        if (a->shortcut() == QKeySequence::Paste) a->trigger();
    CHECK(w.findChild<Canvas*>()->document().atoms.size() == doc.atoms.size());
}

TEST_CASE("paste from ChemDraw: CDX/CDXML clipboard formats (#104)") {
    App app;
    const QString path = QString(PENZENE_TEST_DATA) + "/scheme.cdxml";
    QFile f(path);
    REQUIRE(f.open(QIODevice::ReadOnly));
    const QByteArray cdxml = f.readAll();
    const auto expected = chem::readFile(path);
    REQUIRE(expected);
    // macOS (through ChemDrawPasteboard) and Windows name the format differently.
    for (const char* type : {"chemical/x-cdx", "application/x-qt-windows-mime;value=\"ChemDraw Interchange Format\""}) {
        INFO(type);
        MainWindow w;
        auto* mime = new QMimeData;
        mime->setData(type, cdxml);
        mime->setText("not a structure");  // ChemDraw also offers text; the CDX wins
        QApplication::clipboard()->setMimeData(mime);
        for (auto* a : w.findChildren<QAction*>())
            if (a->shortcut() == QKeySequence::Paste) a->trigger();
        CHECK(w.findChild<Canvas*>()->document().atoms.size() == expected->atoms.size());
    }
#ifdef Q_OS_MACOS
    ChemDrawPasteboard uti;
    CHECK(uti.mimeForUti("com.perkinelmer.chemdraw.cdx-clipboard") == "chemical/x-cdx");
    CHECK(uti.mimeForUti("public.utf8-plain-text").isEmpty());
    CHECK(uti.convertToMime("chemical/x-cdx", {cdxml}, {}).toByteArray() == cdxml);
#endif
}
