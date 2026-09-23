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
