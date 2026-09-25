#pragma once
// Small vector helpers shared by rendering, editing and the canvas.
#include "Document.h"

#include <QPointF>
#include <QtMath>
#include <cmath>

inline double len(QPointF v) { return std::hypot(v.x(), v.y()); }
inline QPointF unit(QPointF v) { double l = len(v); return l > 1e-9 ? v / l : QPointF(1, 0); }
inline QPointF perp(QPointF v) { return {-v.y(), v.x()}; }
inline double cross(QPointF a, QPointF b) { return a.x() * b.y() - a.y() * b.x(); }
inline QPointF dirAt(double deg) { return {std::cos(qDegreesToRadians(deg)), std::sin(qDegreesToRadians(deg))}; }
inline QPointF rotated(QPointF v, double deg) {
    double a = qDegreesToRadians(deg);
    return {v.x() * std::cos(a) - v.y() * std::sin(a), v.x() * std::sin(a) + v.y() * std::cos(a)};
}

using BondsAt = std::vector<std::vector<int>>;  // Document::bondsAt()

inline std::vector<int> neighbors(const Document& doc, const BondsAt& at, int atom) {
    std::vector<int> out;
    for (int b : at[atom]) out.push_back(doc.bonds[b].a == atom ? doc.bonds[b].b : doc.bonds[b].a);
    return out;
}

// sp centre: a triple bond, or two double bonds (allene). Its bonds are collinear.
inline bool isSp(const Document& doc, int atom, const BondsAt& at) {
    int doubles = 0;
    for (int b : at[atom]) {
        if (doc.bonds[b].order == 3) return true;
        doubles += doc.bonds[b].order == 2;
    }
    return doubles >= 2;
}
inline bool isSp(const Document& doc, int atom) { return isSp(doc, atom, doc.bondsAt()); }

// Where a double bond's second line goes: +1 on perp(a->b), -1 opposite, 0 centred.
// Automatic unless the bond says otherwise: toward the neighbours (inside a
// ring), centred at a terminal atom or an sp centre (so C=C=C lines meet).
inline int doubleBondSide(const Document& doc, const Bond& b, const BondsAt& at) {
    if (b.position == BondPosition::Centre) return 0;
    if (b.position != BondPosition::Auto) return b.position == BondPosition::Right ? 1 : -1;
    const auto na = neighbors(doc, at, b.a), nb = neighbors(doc, at, b.b);
    if (na.size() == 1 || nb.size() == 1 || isSp(doc, b.a, at) || isSp(doc, b.b, at)) return 0;
    const QPointF pa = doc.atoms[b.a].pos, d = unit(doc.atoms[b.b].pos - pa);
    double side = 0;
    for (const auto* list : {&na, &nb})
        for (int n : *list)
            if (n != b.a && n != b.b) side += cross(d, doc.atoms[n].pos - pa) > 0 ? 1 : -1;
    // A tie (trans chain) still offsets, or both lines would cross the single bonds.
    return side >= 0 ? 1 : -1;
}
inline int doubleBondSide(const Document& doc, const Bond& b) { return doubleBondSide(doc, b, doc.bondsAt()); }
