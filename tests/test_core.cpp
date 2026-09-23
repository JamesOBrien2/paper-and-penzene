#include "Chem.h"

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
