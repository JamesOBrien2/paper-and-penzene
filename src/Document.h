#pragma once
#include <QByteArray>
#include <QPointF>
#include <QColor>
#include <QString>
#include <optional>
#include <vector>

// Scene units are points; a standard bond is 14.4 pt (ACS 1996).
constexpr double kBondLength = 14.4;

// Wedge/Hash/Wavy carry stereo; Bold and Dashed are drawing styles (on a
// double bond they style one line).
enum class BondStereo { None, Wedge, Hash, Bold, Dashed, Wavy };
enum class BondPosition { Auto, Left, Centre, Right };  // double bond's second line, seen from a to b

struct Atom {
    QPointF pos;
    int z = 6;       // atomic number
    int charge = 0;
    QString label;  // abbreviation such as "Boc"; z/charge are then its attaching atom's
    QColor color;   // invalid: the ink (theme on screen, black in exports); also colours its label
    int map = 0;    // reaction atom-map number (SMILES :n); 0 = none
};

struct Bond {
    int a = 0, b = 0;  // atom indices; stereo points from a to b
    int order = 1;     // 1..3
    BondStereo stereo = BondStereo::None;
    BondPosition position = BondPosition::Auto;
    QColor color;
};

enum class ArrowKind { Reaction, Equilibrium, Resonance, Retro, Fishhook };

// Straight when bend == 0; otherwise a curve whose midpoint sits `bend` points
// to the left of from->to as seen on screen (electron pushing).
struct Arrow {
    QPointF from, to;
    ArrowKind kind = ArrowKind::Reaction;
    double bend = 0;
    QColor color;
    bool operator==(const Arrow&) const = default;
};

// Free text; `pos` is the left end of the first baseline. Digits after a
// letter or bracket render as subscripts (formula style).
struct Text {
    QPointF pos;
    QString text;
    double scale = 1;  // relative to the drawing style's label size
    QColor color;
    bool operator==(const Text&) const = default;
};

// A shaded ring interior (ChemDraw ring fill); atoms in ring order.
struct Fill {
    std::vector<int> atoms;
    QColor color;
    bool operator==(const Fill&) const = default;
};

struct Document {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Arrow> arrows;
    std::vector<Text> texts;
    std::vector<Fill> fills;
    QString style;  // drawing style preset name; empty means ACS 1996
    enum class CarbonLabels { None, Terminal, All } carbonLabels = CarbonLabels::None;  // skeletal by default
    bool hideImplicitH = false;  // labels without their implicit H (NH2 drawn as N)
    bool showStereo = false;  // draw CIP (R)/(S) and (E)/(Z) labels
    bool showAtomNumbers = false;  // draw each atom's index (from 1)
    bool aromaticCircles = false;  // default for every aromatic ring
    QString page;         // a pageSizes() name: laid out at final size; "" = no page
    QPointF pageOrigin;   // the page's top-left corner
    std::vector<std::vector<int>> aromaticCircleOverrides;  // sorted ring atom IDs with the opposite display
    bool operator==(const Document&) const = default;
    bool empty() const { return atoms.empty() && arrows.empty() && texts.empty(); }
    void append(const Document& other, QPointF shift = {});  // atom indices renumbered

    // .penz: {"format":"penzene","version":1,"atoms":[...],"bonds":[...],"arrows":[...],"texts":[...]}
    QByteArray toJson() const;
    static std::optional<Document> fromJson(const QByteArray& data);
    // The drawing Penzene embedded in an exported PNG (text chunk) or SVG (<metadata>).
    static std::optional<Document> fromEmbedded(const QByteArray& file);

    int addAtom(QPointF pos, int z = 6);
    int bondBetween(int a, int b) const;  // bond index or -1
    std::vector<int> neighbors(int atom) const;
    QPointF awayDirection(int atom) const;  // bisects the widest gap between its bonds
    void removeBond(int bond);  // also drops endpoints left isolated
    void removeAtoms(const std::vector<int>& atoms);  // also drops their bonds
};

inline bool operator==(const Atom& x, const Atom& y) {
    return x.pos == y.pos && x.z == y.z && x.charge == y.charge && x.label == y.label && x.color == y.color &&
           x.map == y.map;
}
inline bool operator==(const Bond& x, const Bond& y) {
    return x.a == y.a && x.b == y.b && x.order == y.order && x.stereo == y.stereo &&
           x.position == y.position && x.color == y.color;
}
