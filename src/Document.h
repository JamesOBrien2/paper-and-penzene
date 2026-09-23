#pragma once
#include <QByteArray>
#include <QPointF>
#include <QString>
#include <optional>
#include <vector>

// Scene units are points; a standard bond is 14.4 pt (ACS 1996).
constexpr double kBondLength = 14.4;

enum class BondStereo { None, Wedge, Hash };

struct Atom {
    QPointF pos;
    int z = 6;       // atomic number
    int charge = 0;
    QString label;  // abbreviation such as "Boc"; z/charge are then its attaching atom's
};

struct Bond {
    int a = 0, b = 0;  // atom indices; stereo points from a to b
    int order = 1;     // 1..3
    BondStereo stereo = BondStereo::None;
};

enum class ArrowKind { Reaction, Equilibrium, Resonance, Retro, Fishhook };

// Straight when bend == 0; otherwise a curve whose midpoint sits `bend` points
// to the left of from->to as seen on screen (electron pushing).
struct Arrow {
    QPointF from, to;
    ArrowKind kind = ArrowKind::Reaction;
    double bend = 0;
    bool operator==(const Arrow&) const = default;
};

// Free text; `pos` is the left end of the first baseline. Digits after a
// letter or bracket render as subscripts (formula style).
struct Text {
    QPointF pos;
    QString text;
    bool operator==(const Text&) const = default;
};

struct Document {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Arrow> arrows;
    std::vector<Text> texts;
    bool operator==(const Document&) const = default;
    bool empty() const { return atoms.empty() && arrows.empty() && texts.empty(); }
    void append(const Document& other, QPointF shift = {});  // atom indices renumbered

    // .penz: {"format":"penzene","version":1,"atoms":[...],"bonds":[...],"arrows":[...],"texts":[...]}
    QByteArray toJson() const;
    static std::optional<Document> fromJson(const QByteArray& data);

    int addAtom(QPointF pos, int z = 6);
    int bondBetween(int a, int b) const;  // bond index or -1
    std::vector<int> neighbors(int atom) const;
    QPointF awayDirection(int atom) const;  // bisects the widest gap between its bonds
    void removeAtoms(const std::vector<int>& atoms);  // also drops their bonds
};

inline bool operator==(const Atom& x, const Atom& y) {
    return x.pos == y.pos && x.z == y.z && x.charge == y.charge && x.label == y.label;
}
inline bool operator==(const Bond& x, const Bond& y) {
    return x.a == y.a && x.b == y.b && x.order == y.order && x.stereo == y.stereo;
}
