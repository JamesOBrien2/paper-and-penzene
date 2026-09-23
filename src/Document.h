#pragma once
#include <QByteArray>
#include <QPointF>
#include <optional>
#include <vector>

// Scene units are points; a standard bond is 14.4 pt (ACS 1996).
constexpr double kBondLength = 14.4;

enum class BondStereo { None, Wedge, Hash };

struct Atom {
    QPointF pos;
    int z = 6;       // atomic number
    int charge = 0;
};

struct Bond {
    int a = 0, b = 0;  // atom indices; stereo points from a to b
    int order = 1;     // 1..3
    BondStereo stereo = BondStereo::None;
};

struct Document {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    bool operator==(const Document&) const = default;

    // .penz: {"format":"penzene","version":1,"atoms":[...],"bonds":[...]}
    QByteArray toJson() const;
    static std::optional<Document> fromJson(const QByteArray& data);
};

inline bool operator==(const Atom& x, const Atom& y) {
    return x.pos == y.pos && x.z == y.z && x.charge == y.charge;
}
inline bool operator==(const Bond& x, const Bond& y) {
    return x.a == y.a && x.b == y.b && x.order == y.order && x.stereo == y.stereo;
}
