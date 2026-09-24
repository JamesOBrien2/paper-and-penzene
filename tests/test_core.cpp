#include "Chem.h"
#include "Edit.h"

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("SMILES gives 2D coordinates") {
    auto doc = chem::fromSmiles("c1ccccc1O");
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 7);
    REQUIRE(doc->bonds.size() == 7);
    int doubles = 0;
    for (auto& b : doc->bonds) doubles += b.order == 2;
    CHECK(doubles == 3);  // kekulized
    auto& b = doc->bonds[0];
    auto d = doc->atoms[b.a].pos - doc->atoms[b.b].pos;
    CHECK(std::abs(std::hypot(d.x(), d.y()) - kBondLength) < 0.5);
}

TEST_CASE("bad SMILES is rejected") {
    CHECK_FALSE(chem::fromSmiles("C1CC"));
}

TEST_CASE(".penz round-trips") {
    auto doc = chem::fromSmiles("C[C@H](N)C(=O)[O-]");
    REQUIRE(doc);
    doc->bonds[0].stereo = BondStereo::Wedge;
    doc->arrows.push_back({{0, 0}, {40, 0}, ArrowKind::Equilibrium});
    doc->arrows.push_back({{0, 10}, {20, 10}, ArrowKind::Fishhook, -6});
    doc->texts.push_back({{5, -8}, "Pd(PPh3)4\n80 °C"});
    doc->style = "JDP";
    auto back = Document::fromJson(doc->toJson());
    REQUIRE(back);
    CHECK(*back == *doc);
}

TEST_CASE(".penz rejects junk and dangling bonds") {
    CHECK_FALSE(Document::fromJson("not json"));
    CHECK_FALSE(Document::fromJson(R"({"format":"penzene","version":1,
        "atoms":[{"x":0,"y":0,"z":6}],"bonds":[{"a":0,"b":5}]})"));
}

TEST_CASE("MOL round-trips atoms, bonds, coordinates and wedges") {
    auto doc = chem::fromSmiles("C[C@H](N)C(=O)O");
    REQUIRE(doc);
    int wedged = 0;
    for (auto& b : doc->bonds) wedged += b.stereo != BondStereo::None;
    CHECK(wedged == 1);
    auto back = chem::fromMolBlock(chem::toMolBlock(*doc));
    REQUIRE(back);
    REQUIRE(back->atoms.size() == doc->atoms.size());
    REQUIRE(back->bonds.size() == doc->bonds.size());
    for (size_t i = 0; i < doc->atoms.size(); ++i) {
        CHECK(back->atoms[i].z == doc->atoms[i].z);
        auto d = back->atoms[i].pos - doc->atoms[i].pos;
        CHECK(std::hypot(d.x(), d.y()) < 0.01);
    }
    for (size_t i = 0; i < doc->bonds.size(); ++i) CHECK(back->bonds[i] == doc->bonds[i]);
    CHECK(chem::toSmiles(*back) == chem::toSmiles(*doc));
    CHECK(chem::toSmiles(*doc) == "C[C@H](N)C(=O)O");
}

TEST_CASE("implicit hydrogens and valence errors") {
    Document d;
    d.atoms = {{{0, 0}, 8}, {{kBondLength, 0}, 6}};
    d.bonds = {{0, 1}};
    auto info = chem::atomInfo(d);
    CHECK(info[0].hydrogens == 1);  // OH
    CHECK(info[1].hydrogens == 3);  // CH3
    d.bonds[0].order = 3;           // O#C: oxygen over valence
    CHECK(chem::atomInfo(d)[0].valenceError);
}

TEST_CASE("clean keeps atom order and centroid") {
    Document d;
    for (int i = 0; i < 6; ++i) d.atoms.push_back({QPointF(100 + i * 3, 50 + i * i)});
    for (int i = 0; i < 6; ++i) d.bonds.push_back({i, (i + 1) % 6, i % 2 ? 2 : 1});
    auto c = chem::clean2D(d);
    REQUIRE(c.atoms.size() == 6);
    auto e = c.atoms[0].pos - c.atoms[1].pos;
    CHECK(std::abs(std::hypot(e.x(), e.y()) - kBondLength) < 0.5);
    QPointF c0, c1;
    for (int i = 0; i < 6; ++i) c0 += d.atoms[i].pos, c1 += c.atoms[i].pos;
    CHECK(std::hypot((c0 - c1).x(), (c0 - c1).y()) < 0.01);
}

TEST_CASE("element symbols") {
    CHECK(chem::symbol(17) == "Cl");
    CHECK(chem::atomicNumber("Br") == 35);
    CHECK(chem::atomicNumber("Xx") == 0);
}

TEST_CASE(".penz rejects bad arrows; v0.1 files still load") {
    CHECK_FALSE(Document::fromJson(R"({"format":"penzene","version":1,"arrows":[{"kind":"wiggly"}]})"));
    auto old = Document::fromJson(R"({"format":"penzene","version":1,"atoms":[{"x":0,"y":0,"z":6}],"bonds":[]})");
    REQUIRE(old);
    CHECK(old->arrows.empty());
}

TEST_CASE("clean lays out each fragment in place and keeps arrows and text") {
    auto left = chem::fromSmiles("CCO"), right = chem::fromSmiles("CC=O");
    REQUIRE(left);
    REQUIRE(right);
    Document scheme;
    scheme.append(*left, {-100, 0});
    scheme.append(*right, {100, 0});
    scheme.arrows.push_back({{-30, 0}, {30, 0}});
    scheme.texts.push_back({{-20, -5}, "PCC"});
    for (auto& a : scheme.atoms) a.pos += QPointF(0, a.pos.x() * 0.1);  // skew it
    Document clean = chem::clean2D(scheme);
    CHECK(clean.arrows == scheme.arrows);
    CHECK(clean.texts == scheme.texts);
    CHECK(clean.bonds.size() == scheme.bonds.size());
    for (int i : {0, 1, 2}) CHECK(clean.atoms[i].pos.x() < -50);  // reactant stays left
    for (int i : {3, 4, 5}) CHECK(clean.atoms[i].pos.x() > 50);

    // Only the molecule with a selected atom moves; the other is untouched.
    Document partial = chem::clean2D(scheme, {4});
    for (int i : {0, 1, 2}) CHECK(partial.atoms[i] == scheme.atoms[i]);
    CHECK(partial.atoms[3] == clean.atoms[3]);  // the product is laid out as in a full clean
    CHECK(partial.atoms[5] == clean.atoms[5]);
    CHECK(partial.bonds.size() == scheme.bonds.size());
}

TEST_CASE("formula, weights and InChI") {
    auto aspirin = chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
    REQUIRE(aspirin);
    auto p = chem::properties(*aspirin);
    REQUIRE(p);
    CHECK(p->formula == "C9H8O4");
    CHECK(std::abs(p->mw - 180.159) < 0.01);
    CHECK(std::abs(p->exactMass - 180.0423) < 0.001);
    CHECK(chem::toInchiKey(*aspirin) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N");
    CHECK(chem::toInchi(*aspirin).rfind("InChI=1S/C9H8O4/", 0) == 0);
    CHECK_FALSE(chem::properties(Document{}));
}

TEST_CASE("CDXML import: molecules, arrows and text in place") {
    auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/scheme.cdxml");
    REQUIRE(doc);
    CHECK(chem::toSmiles(*doc) == "CCO");
    REQUIRE(doc->atoms.size() == 3);
    CHECK(std::abs(doc->atoms[0].pos.x() - 100) < 0.5);  // 14.4 pt bonds: points map 1:1
    CHECK(std::abs(doc->atoms[0].pos.y() - 100) < 0.5);
    REQUIRE(doc->arrows.size() == 3);
    CHECK(doc->arrows[0].kind == ArrowKind::Reaction);
    CHECK(doc->arrows[0].to == QPointF(190, 104));
    CHECK(doc->arrows[1].kind == ArrowKind::Equilibrium);
    // Quarter circle of radius 25 about (165,175): its top is 7.3 pt above the chord.
    CHECK(std::abs(doc->arrows[2].bend - 7.32) < 0.1);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].text == "PCC");
}

TEST_CASE("hotkeys without a canvas: ChemDraw's dipeptide example") {
    Document doc;
    int n = doc.addAtom({0, 0}, 7);
    edit::link(doc, n, doc.addAtom({kBondLength, 0}));
    edit::Hotspot h{1, -1};
    for (QChar k : QString("42n152o")) {
        h = edit::hotkey(doc, h, k);
        REQUIRE(h.valid());
    }
    std::string smi = chem::toSmiles(doc);
    std::erase(smi, '@');
    CHECK(chem::toSmiles(*chem::fromSmiles(smi)) == chem::toSmiles(*chem::fromSmiles("CC(N)C(=O)NC(C)C(=O)O")));
    CHECK_FALSE(edit::hotkey(doc, {0, -1}, "~").valid());  // not a hotkey
}

TEST_CASE("clean keeps bond display styles and double-bond positions (#84)") {
    auto d = chem::fromSmiles("CC=CC(C)C");
    REQUIRE(d);
    int dbl = -1, single = -1;
    for (int i = 0; i < int(d->bonds.size()); ++i)
        (d->bonds[i].order == 2 ? dbl : single) = i;
    d->bonds[dbl].position = BondPosition::Right;
    d->bonds[single].stereo = BondStereo::Bold;
    Document clean = chem::clean2D(*d);
    const Bond& b1 = clean.bonds[clean.bondBetween(d->bonds[dbl].a, d->bonds[dbl].b)];
    const Bond& b2 = clean.bonds[clean.bondBetween(d->bonds[single].a, d->bonds[single].b)];
    CHECK(b1.position == BondPosition::Right);
    CHECK(b2.stereo == BondStereo::Bold);
}

TEST_CASE("clean lays abbreviations out as single nodes, so bonds stay even (#85)") {
    // Found by fuzzing: "79P7" puts a Ph label on a crowded atom. Expanding the
    // ring before layout squeezed some bonds to about half length.
    auto d = std::make_optional<Document>();
    d->addAtom({0, 0});
    edit::Hotspot h{0, -1};
    for (QChar k : QString("79P7")) h = edit::hotkey(*d, h, k);
    Document clean = chem::clean2D(*d);
    REQUIRE(clean.atoms.size() == d->atoms.size());
    for (const auto& b : clean.bonds) {
        QPointF v = clean.atoms[b.a].pos - clean.atoms[b.b].pos;
        CHECK(std::abs(std::hypot(v.x(), v.y()) - kBondLength) < 0.2 * kBondLength);
    }
}

TEST_CASE("CDXML label nodes: reagent labels become text, R groups stay labelled (#86)") {
    auto doc = chem::readFile(QString(PENZENE_TEST_DATA) + "/labels.cdxml");
    REQUIRE(doc);
    REQUIRE(doc->atoms.size() == 2);  // the carbon and R; no stray atoms from the reagent label
    CHECK(doc->atoms[1].label == "R");
    CHECK(doc->bonds.size() == 1);
    REQUIRE(doc->texts.size() == 1);
    CHECK(doc->texts[0].text == "LiBr, acetone");
    CHECK(std::abs(doc->texts[0].scale - 0.7) < 0.01);  // 7 pt against the 10 pt default
    CHECK(std::abs(doc->texts[0].pos.x() - 150) < 0.5);
    auto back = Document::fromJson(doc->toJson());  // text scale survives .penz
    REQUIRE(back);
    CHECK(*back == *doc);
}

TEST_CASE("explicit hydrogens and carbon/H display options (#97)") {
    auto eth = chem::fromSmiles("CCO");
    REQUIRE(eth);
    Document withH = chem::addHydrogens(*eth);
    CHECK(withH.atoms.size() == 9);  // C2H6O: 3 heavy + 6 H
    CHECK(chem::properties(withH)->formula == "C2H6O");  // chemistry unchanged
    for (const auto& b : withH.bonds) {  // placed at a bond's length, not piled up
        QPointF v = withH.atoms[b.a].pos - withH.atoms[b.b].pos;
        CHECK(std::hypot(v.x(), v.y()) > 0.4 * kBondLength);
    }
    Document without = chem::removeHydrogens(withH);
    CHECK(without.atoms.size() == 3);

    Document wedgedH = withH;  // a wedged H carries stereo, so it stays
    for (auto& b : wedgedH.bonds)
        if (wedgedH.atoms[b.b].z == 1) { b.stereo = BondStereo::Wedge; break; }
    CHECK(chem::removeHydrogens(wedgedH).atoms.size() == 4);

    eth->carbonLabels = Document::CarbonLabels::Terminal;
    eth->hideImplicitH = true;
    auto back = Document::fromJson(eth->toJson());
    REQUIRE(back);
    CHECK(back->carbonLabels == Document::CarbonLabels::Terminal);
    CHECK(back->hideImplicitH);
}

TEST_CASE("CIP stereo labels, and E/Z read from the drawing (#94)") {
    auto ala = chem::fromSmiles("C[C@H](N)C(=O)O");  // L-alanine
    REQUIRE(ala);
    auto labels = chem::stereoLabels(*ala);
    REQUIRE(labels.size() == 1);
    CHECK(labels[0].atom == 1);
    CHECK(labels[0].text == "S");

    auto ene = chem::fromSmiles("C/C=C/C");
    REQUIRE(ene);
    labels = chem::stereoLabels(*ene);
    REQUIRE(labels.size() == 1);
    CHECK(labels[0].bond >= 0);
    CHECK(labels[0].text == "E");

    // A plain zig-zag drawing of 2-butene is trans: SMILES now says so.
    Document drawn;
    drawn.atoms = {{{0, 0}}, {{12.47, -7.2}}, {{24.94, 0}}, {{37.41, -7.2}}};
    drawn.bonds = {{0, 1}, {1, 2, 2}, {2, 3}};
    CHECK(chem::toSmiles(drawn) == "C/C=C/C");

    ala->showStereo = true;
    auto back = Document::fromJson(ala->toJson());
    REQUIRE(back);
    CHECK(back->showStereo);
}

TEST_CASE("check structure finds valence, stereo, label and overlap problems (#95)") {
    auto has = [](const std::vector<chem::Problem>& ps, const QString& text) {
        return std::any_of(ps.begin(), ps.end(), [&](const auto& p) { return p.message.contains(text); });
    };
    auto butanol = chem::fromSmiles("CCC(C)O");  // a stereocentre drawn without a wedge
    REQUIRE(butanol);
    auto ps = chem::checkStructure(*butanol);
    CHECK(has(ps, "no wedge"));

    Document d;  // ethane with a wedge (not a stereocentre), a bad label, overlapping atoms
    d.atoms = {{{0, 0}}, {{kBondLength, 0}}, {{40, 0}}, {{40.5, 0}}};
    d.bonds = {{0, 1, 1, BondStereo::Wedge}};
    d.atoms[2].label = "Xyz";
    ps = chem::checkStructure(d);
    CHECK(has(ps, "not a stereocentre"));
    CHECK(has(ps, "Unknown label"));
    CHECK(has(ps, "Overlapping"));

    auto pentavalent = chem::fromSmiles("C");
    Document c5 = *pentavalent;
    for (int k = 0; k < 5; ++k) edit::link(c5, 0, c5.addAtom({10.0 * k, 10}));
    CHECK(has(chem::checkStructure(c5), "Valence error"));

    CHECK(chem::checkStructure(*chem::fromSmiles("CC(=O)Oc1ccccc1C(=O)O")).empty());  // aspirin is fine
}
