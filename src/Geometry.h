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

// sp centre: a triple bond, or two double bonds (allene). Its bonds are collinear.
inline bool isSp(const Document& doc, int atom) {
    int doubles = 0;
    for (const auto& b : doc.bonds)
        if (b.a == atom || b.b == atom) {
            if (b.order == 3) return true;
            doubles += b.order == 2;
        }
    return doubles >= 2;
}
