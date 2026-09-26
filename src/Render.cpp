#include "Render.h"
#include "Chem.h"
#include "Geometry.h"

#include <QBuffer>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QPainter>
#include <QSet>
#include <QPainterPath>
#include <QPdfWriter>
#include <QSvgGenerator>
#include <algorithm>
#include <limits>
#include <numbers>
#include <numeric>

// Presets. Values come from the ChemDraw stationery (.cds) of the same name;
// wedge width, hash spacing and label gap keep ACS's proportions to ours.
// Model units per point for a style's native values (bond length / 14.4 undone).
constexpr double kJdp = kBondLength / 14.17, kRsc = kBondLength / 12.2;

const std::vector<DrawingStyle>& drawingStyles() {
    // Wedge width and hash spacing keep ACS's proportion to the bold width and
    // hash setting (4.5/2.0 and 2.2/2.5); label radius is half a glyph plus the margin.
    static const std::vector<DrawingStyle> styles{
        {"ACS 1996", 14.4, 0.6, 2.0, 4.5, 2.2, 0.18, 5.5, "Arial", QFont::Normal, 10},
        // JDPReport.cds: bond 14.17 (0.5 cm), line 0.879, bold 1.814, hash 1.814, margin 1.162,
        // IBM Plex Sans Light 10 pt.
        {"JDP", 14.17, 0.879 * kJdp, 1.814 * kJdp, 4.5 * 1.814 / 2.0 * kJdp, 2.2 * 1.814 / 2.5 * kJdp, 0.18,
         (3.9 + 1.162) * kJdp, "IBM Plex Sans", QFont::Light, 10 * kJdp},
        // RSC (1 Column).cds (2 Column only differs in page size): bond 12.2, line 0.449, bold 1.602,
        // hash 1.75, margin 1.25, spacing 20%, Helvetica 7 pt.
        {"RSC", 12.2, 0.449 * kRsc, 1.602 * kRsc, 4.5 * 1.602 / 2.0 * kRsc, 2.2 * 1.75 / 2.5 * kRsc, 0.20,
         (3.9 * 0.7 + 1.25) * kRsc, "Helvetica", QFont::Normal, 7 * kRsc},
    };
    return styles;
}

const DrawingStyle& drawingStyle(const QString& name) {
    for (const auto& s : drawingStyles())
        if (s.name == name) return s;
    return drawingStyles()[0];
}
// Carbons are skeletal unless alone or the document asks for their labels.
static bool hasLabel(const Document& doc, int i, const std::vector<int>& degree) {
    if (doc.atoms[i].z == 0 && doc.atoms[i].label.isEmpty()) return false;  // attachment point, η centroid: a bare point
    if (doc.atoms[i].z != 6 || degree[i] == 0 || !doc.atoms[i].label.isEmpty()) return true;
    return doc.carbonLabels == Document::CarbonLabels::All ||
           (doc.carbonLabels == Document::CarbonLabels::Terminal && degree[i] == 1);
}

// ---------------------------------------------------------------- rendering

QFont labelFont(const DrawingStyle& s, double scale) {
    QFont f(s.font);
    f.setWeight(s.weight);
    f.setPixelSize(int(s.fontSize * scale));  // 1 px == 1 pt in our scene units
    return f;
}

// Text as outlines, so it scales identically on screen, SVG, PDF and PNG.
static void drawText(QPainter& p, const QString& s, QPointF baselineLeft, const QFont& f) {
    QPainterPath path;
    path.addText(baselineLeft, f, s);
    p.fillPath(path, p.pen().color());
}

// Abbreviation written from the right, bond side last: OMe -> MeO.
static QString reversedLabel(const QString& s) {
    static const QHash<QString, QString> r{
        {"OMe", "MeO"}, {"CO2Me", "MeO2C"}, {"CO2Et", "EtO2C"}, {"CO2H", "HO2C"}, {"NO2", "O2N"},
        {"CF3", "F3C"}, {"OAc", "AcO"},     {"CHO", "OHC"},     {"SO2Me", "MeO2S"}, {"OTf", "TfO"},
        {"OTs", "TsO"}, {"OTBS", "TBSO"},   {"CN", "NC"},       {"Bpin", "pinB"},
    };
    return r.value(s, s);
}

// Abbreviation centred on its attaching letter (the first, or the last when written from the right).
static void drawAbbreviation(QPainter& p, const Atom& a, bool fromRight, const DrawingStyle& st) {
    QFontMetricsF fm(labelFont(st));
    const QString s = fromRight ? reversedLabel(a.label) : a.label;
    const double base = a.pos.y() + fm.capHeight() / 2;
    QPainterPath path = textPath({{0, base}, s}, st);
    const double x = fromRight ? a.pos.x() + fm.horizontalAdvance(s.back()) / 2 - path.boundingRect().right()
                               : a.pos.x() - fm.horizontalAdvance(s.front()) / 2;
    p.fillPath(path.translated(x, 0), p.pen().color());
}

// Where a label's implicit H goes: after the symbol, before it, or stacked
// below/above when bonds leave neither side free (a middle CH2).
enum class HSide { Right, Left, Below, Above };

static void drawLabel(QPainter& p, const Document& doc, int i, int hydrogens, HSide side, const DrawingStyle& st) {
    const bool hLeft = side == HSide::Left;
    const auto& a = doc.atoms[i];
    if (!a.label.isEmpty()) return drawAbbreviation(p, a, hLeft, st);
    QFont f = labelFont(st), sub = labelFont(st, 0.7);
    QFontMetricsF fm(f), sm(sub);
    QString sym = QString::fromStdString(chem::symbol(a.z));
    double w = fm.horizontalAdvance(sym);
    double base = a.pos.y() + fm.capHeight() / 2;
    double x = a.pos.x() - w / 2;
    drawText(p, sym, {x, base}, f);

    double right = x + w;
    if (hydrogens > 0) {
        QString n = hydrogens > 1 ? QString::number(hydrogens) : QString();
        double hw = fm.horizontalAdvance("H"), nw = sm.horizontalAdvance(n);
        if (side == HSide::Below || side == HSide::Above) {
            const double hx = a.pos.x() - (hw + nw) / 2, dy = fm.capHeight() * 1.35 * (side == HSide::Below ? 1 : -1);
            drawText(p, "H", {hx, base + dy}, f);
            if (!n.isEmpty()) drawText(p, n, {hx + hw, base + dy + fm.capHeight() * 0.35}, sub);
        } else {
            double hx = hLeft ? x - hw - nw : right;
            drawText(p, "H", {hx, base}, f);
            if (!n.isEmpty()) drawText(p, n, {hx + hw, base + fm.capHeight() * 0.35}, sub);
            if (!hLeft) right += hw + nw;
        }
    }
    if (a.charge) {
        QString c = QString(a.charge > 0 ? "+" : "−");
        if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
        drawText(p, c, {right, base - fm.capHeight() * 0.7}, sub);
    }
}

static void drawBond(QPainter& p, const Document& doc, const Bond& b, const DrawingStyle& st, const std::vector<int>& degree,
                     const std::vector<bool>& labeled, const BondsAt& at, const std::vector<QPointF>& gaps = {}) {
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos;
    QPointF d = unit(pb - pa), n = perp(d);
    const double gap = st.bondSpacing * kBondLength;  // double-bond spacing
    // Trim at labels.
    QPointF a = labeled[b.a] ? pa + d * st.labelRadius : pa;
    QPointF e = labeled[b.b] ? pb - d * st.labelRadius : pb;

    if (b.stereo == BondStereo::Wedge) {
        QPolygonF tri{a, e + n * st.wedgeWidth / 2, e - n * st.wedgeWidth / 2};
        p.setBrush(p.pen().color());
        p.drawPolygon(tri);
        p.setBrush(Qt::NoBrush);
        return;
    }
    if (b.stereo == BondStereo::Hash) {
        double L = len(e - a);
        int count = std::max(3, int(L / st.hashSpacing));
        for (int k = 0; k <= count; ++k) {
            double t = double(k) / count;
            QPointF c = a + (e - a) * t;
            double w = st.wedgeWidth / 2 * t;
            p.drawLine(c + n * w, c - n * w);
        }
        return;
    }

    if (b.stereo == BondStereo::Interaction) {  // dotted, whatever the order
        QPen dots = p.pen();
        dots.setCapStyle(Qt::RoundCap);
        dots.setDashPattern({0.01, 2.2});
        p.setPen(dots);
        p.drawLine(a, e);
        return;
    }
    if (b.stereo == BondStereo::Wavy) {
        QPainterPath wave(a);
        const double L = len(e - a);
        const int bumps = std::max(2, int(std::round(L / 3)));
        for (int k = 0; k < bumps; ++k) {
            QPointF from = a + (e - a) * (double(k) / bumps), to = a + (e - a) * (double(k + 1) / bumps);
            wave.quadTo((from + to) / 2 + n * (k % 2 ? -2.0 : 2.0) * 1.3, to);
        }
        p.drawPath(wave);
        return;
    }
    // Bold styles the main line, Dashed and Partial the other one (or the only one).
    const QPen pen = p.pen();
    auto line = [&](QPointF x, QPointF y, bool main) {
        QPen q = pen;
        if (b.stereo == BondStereo::Bold && main) q.setWidthF(st.boldWidth), q.setCapStyle(Qt::FlatCap);
        const bool dashed = b.stereo == BondStereo::Dashed || b.stereo == BondStereo::Partial;
        if (dashed && (!main || b.order == 1)) q.setDashPattern({2.5, 2.5});
        p.setPen(q);
        // A bond drawn later (in front) crossing this one leaves a gap in it.
        std::vector<double> cuts;  // along x→y, 0..1
        const double length = len(y - x), half = 0.18 * kBondLength;
        for (QPointF g : gaps) cuts.push_back(QPointF::dotProduct(g - x, y - x) / (length * length));
        std::sort(cuts.begin(), cuts.end());
        double from = 0;
        for (double c : cuts) {
            const double to = c - half / length;
            if (to > from) p.drawLine(x + (y - x) * from, x + (y - x) * to);
            from = std::max(from, c + half / length);
        }
        if (from < 1) p.drawLine(x + (y - x) * from, y);
        p.setPen(pen);
    };

    if (b.order == 1) {
        line(a, e, true);
    } else if (b.order == 3) {
        p.drawLine(a, e);
        p.drawLine(a + n * gap, e + n * gap);
        p.drawLine(a - n * gap, e - n * gap);
    } else {
        // Offset the second line toward the side where the neighbours are
        // (inside the ring); centre it for terminal bonds like C=O.
        const int side = doubleBondSide(doc, b, at);
        if (side == 0) {
            QPointF o = n * gap / 2;
            line(a + o, e + o, true);
            line(a - o, e - o, false);
        } else {
            QPointF o = n * (side >= 0 ? gap : -gap);
            QPointF shrink = d * (0.15 * kBondLength);
            QPointF ia = labeled[b.a] ? a : a + shrink, ie = labeled[b.b] ? e : e - shrink;
            line(a, e, true);
            line(ia + o, ie + o, false);
        }
    }

    // A bold line ends flat, so at a skeletal atom fill the bevel to each other bond there, as a pen
    // joins a polyline: no notch at ring corners (#217). Labels keep their clearance.
    if (b.stereo != BondStereo::Bold || (b.order == 2 && doubleBondSide(doc, b, at) == 0) || b.order == 3) return;
    p.save();
    p.setPen(Qt::NoPen);
    p.setBrush(pen.color());
    for (int end : {b.a, b.b}) {
        if (labeled[end]) continue;
        const QPointF c = doc.atoms[end].pos, out = unit(doc.atoms[end == b.a ? b.b : b.a].pos - c);
        for (int other : at[end]) {
            const Bond& o = doc.bonds[other];
            if (std::minmax(o.a, o.b) == std::minmax(b.a, b.b)) continue;  // this bond (b may be a copy)
            const QPointF dir = unit(doc.atoms[o.a == end ? o.b : o.a].pos - c);
            const double w = st.boldWidth / 2, w2 = o.stereo == BondStereo::Bold ? w : pen.widthF() / 2;
            std::vector<QPointF> quad{c + perp(out) * w, c - perp(out) * w, c + perp(dir) * w2, c - perp(dir) * w2};
            std::sort(quad.begin(), quad.end(), [&](QPointF x, QPointF y) {
                return std::atan2(x.y() - c.y(), x.x() - c.x()) < std::atan2(y.y() - c.y(), y.x() - c.x());
            });
            p.drawPolygon(QPolygonF(QList<QPointF>(quad.begin(), quad.end())));
        }
    }
    p.restore();
}

// ---- arrows and text

constexpr double kHeadLength = 6, kHeadWidth = 2.2, kEquilibriumGap = 1.6;

// Quadratic control point: puts the curve's midpoint `bend` to the left of from->to.
// Points along the arrow. A curve is a circular arc through both ends whose
// midpoint sits `bend` off the chord, so arcs past 180° (ChemDraw's cycle
// arrows) are exact too.
static std::vector<QPointF> arrowPoints(const Arrow& a) {
    if (!a.bend || len(a.to - a.from) < 1e-6) return {a.from, a.to};
    const QPointF mid = (a.from + a.to) / 2, u = -perp(unit(a.to - a.from)) * (a.bend > 0 ? 1 : -1);
    const double h = len(a.to - a.from) / 2, sag = std::abs(a.bend), radius = (h * h + sag * sag) / (2 * sag);
    const QPointF top = mid + u * sag, c = top - u * radius;
    auto angle = [&](QPointF p) { return std::atan2(p.y() - c.y(), p.x() - c.x()); };
    const double from = angle(a.from), through = angle(top);
    double sweep = angle(a.to) - from;
    // Go the way that passes through the arc's midpoint.
    auto wrap = [](double x) { return std::remainder(x, 2 * M_PI); };
    const double half = wrap(through - from);
    if (half > 0 && sweep < 0) sweep += 2 * M_PI;
    if (half < 0 && sweep > 0) sweep -= 2 * M_PI;
    const int n = std::max(8, int(std::abs(sweep) * radius / 1.5));  // about 1.5 pt per segment
    std::vector<QPointF> pts;
    for (int k = 0; k <= n; ++k) {
        double t = from + sweep * k / n;
        pts.push_back(c + QPointF(std::cos(t), std::sin(t)) * radius);
    }
    pts.front() = a.from, pts.back() = a.to;
    return pts;
}

QPainterPath arrowPath(const Arrow& a) {
    if (isShape(a.kind) && a.kind != ArrowKind::Line) {
        const QRectF r = QRectF(a.from, a.to).normalized();
        QPainterPath path;
        if (a.kind == ArrowKind::Box) path.addRect(r);
        else if (a.kind == ArrowKind::RoundedBox) path.addRoundedRect(r, 0.35 * kBondLength, 0.35 * kBondLength);
        else path.addEllipse(r);
        return path;
    }
    const auto pts = arrowPoints(a);
    QPainterPath path(pts[0]);
    for (size_t k = 1; k < pts.size(); ++k) path.lineTo(pts[k]);
    return path;
}

// Drops `by` points' worth of length from the end of a polyline.
static void trimEnd(std::vector<QPointF>& pts, double by) {
    while (pts.size() > 2 && len(pts.back() - pts[pts.size() - 2]) <= by)
        by -= len(pts.back() - pts[pts.size() - 2]), pts.pop_back();
    pts.back() -= unit(pts.back() - pts[pts.size() - 2]) * std::min(by, len(pts.back() - pts[pts.size() - 2]) - 0.01);
}

// Where a polyline is `by` back from its end, measured along it.
static QPointF pointBack(const std::vector<QPointF>& pts, double by) {
    for (size_t i = pts.size() - 1; i > 0; --i) {
        const double seg = len(pts[i] - pts[i - 1]);
        if (seg >= by) return pts[i] - unit(pts[i] - pts[i - 1]) * by;
        by -= seg;
    }
    return pts.front();
}

// Filled head at `tip` pointing along `dir`; `sides` +1/-1 for a half head (one barb, and no
// sliver back along the shaft).
static void drawHead(QPainter& p, QPointF tip, QPointF dir, int sides = 0) {
    QPointF d = unit(dir), n = perp(d), base = tip - d * kHeadLength, notch = tip - d * (kHeadLength * 0.8);
    QPolygonF head = sides > 0   ? QPolygonF{tip, base + n * kHeadWidth, notch}
                     : sides < 0 ? QPolygonF{tip, notch, base - n * kHeadWidth}
                                 : QPolygonF{tip, base + n * kHeadWidth, notch, base - n * kHeadWidth};
    p.setBrush(p.pen().color());
    p.drawPolygon(head);
    p.setBrush(Qt::NoBrush);
}

static void drawArrow(QPainter& p, const Arrow& a) {
    if (a.dashed) {
        QPen dashes = p.pen();
        dashes.setDashPattern({3, 2.5});
        p.setPen(dashes);
    }
    if (isShape(a.kind)) {
        p.drawPath(arrowPath(a));
        return;
    }
    QPointF d = unit(a.to - a.from), n = perp(d);
    if (a.kind == ArrowKind::Equilibrium) {  // ⇌: two half-headed lines
        QPointF o = n * kEquilibriumGap;
        p.drawLine(a.from - o, a.to - o - d * kHeadLength * 0.8);
        drawHead(p, a.to - o, d, -1);
        p.drawLine(a.to + o, a.from + o + d * kHeadLength * 0.8);
        drawHead(p, a.from + o, -d, -1);
        return;
    }
    if (a.kind == ArrowKind::Retro) {  // ⇒: open double arrow
        QPointF o = n * kEquilibriumGap, back = a.to - d * kHeadLength;
        p.drawLine(a.from + o, back + o + d * kEquilibriumGap);
        p.drawLine(a.from - o, back - o + d * kEquilibriumGap);
        p.drawPolyline(QPolygonF{back + n * (kHeadWidth + kEquilibriumGap), a.to, back - n * (kHeadWidth + kEquilibriumGap)});
        return;
    }
    // Heads follow the tangent at each end; the shaft stops inside them so it
    // doesn't poke through the tip.
    // A head is aimed along the chord it covers, not the tangent at the tip: on a tight curve the
    // tangent turns the head off the shaft, which then leaves through one barb.
    std::vector<QPointF> pts = arrowPoints(a);
    const QPointF endDir = pts.back() - pointBack(pts, kHeadLength * 0.8);
    const QPointF startDir = pts.front() - pointBack(std::vector<QPointF>(pts.rbegin(), pts.rend()), kHeadLength * 0.8);
    trimEnd(pts, kHeadLength * 0.7);
    if (a.kind == ArrowKind::Resonance) {
        std::reverse(pts.begin(), pts.end());
        trimEnd(pts, kHeadLength * 0.7);
    }
    p.drawPolyline(pts.data(), int(pts.size()));
    // Fishhook: the barb sits on the outside of the curve.
    drawHead(p, a.to, endDir, a.kind == ArrowKind::Fishhook ? (a.bend >= 0 ? -1 : 1) : 0);
    if (a.kind == ArrowKind::Resonance) drawHead(p, a.from, startDir);
}

static bool subscripted(const QString& s, int i, bool prevSub) {
    if (!s[i].isDigit() || i == 0) return false;
    QChar c = s[i - 1];
    return prevSub || c.isLetter() || c == ')' || c == ']';
}

// Text as outlines, formula-style subscripts, one line per '\n'. Laid out in
// runs (not per letter) so kerning and spaces match ordinary text.
QPainterPath textPath(const Text& t, const DrawingStyle& st) {
    QFont f = labelFont(st, t.scale), sub = labelFont(st, 0.7 * t.scale);
    QFontMetricsF fm(f), sm(sub);
    const double tab = kTabSpaces * fm.horizontalAdvance(' ');
    QPainterPath path;
    const auto lines = t.text.split('\n');
    for (int li = 0; li < lines.size(); ++li) {
        const QString& s = lines[li];
        double x = 0, y = t.pos.y() + li * fm.lineSpacing();
        bool sub_ = false;
        for (int i = 0; i < s.size();) {
            if (s[i] == '\t') {
                x = (std::floor(x / tab + 1e-6) + 1) * tab;
                ++i, sub_ = false;
                continue;
            }
            if (s[i] == ' ') {  // by hand: some platforms drop leading spaces from a shaped run
                x += fm.horizontalAdvance(' ');
                ++i, sub_ = false;
                continue;
            }
            const bool runSub = subscripted(s, i, sub_);
            int j = i + 1;
            while (j < s.size() && s[j] != '\t' && s[j] != ' ' && subscripted(s, j, runSub) == runSub) ++j;
            sub_ = runSub;
            const QString run = s.mid(i, j - i);
            path.addText(t.pos.x() + x, runSub ? y + fm.capHeight() * 0.35 : y, runSub ? sub : f, run);
            x += (runSub ? sm : fm).horizontalAdvance(run);
            i = j;
        }
    }
    return path;
}

void paintDocument(QPainter& p, const Document& doc, const RenderStyle& style) {
    const DrawingStyle& st = drawingStyle(doc.style);
    const double lineWidth = style.lineWidth > 0 ? style.lineWidth : st.lineWidth;
    p.save();
    p.setRenderHint(QPainter::Antialiasing);
    const BondsAt bondsAt = doc.bondsAt();
    std::vector<int> degree(doc.atoms.size(), 0);
    for (const auto& b : doc.bonds) ++degree[b.a], ++degree[b.b];
    std::vector<bool> labeled(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) labeled[i] = hasLabel(doc, int(i), degree);
    auto info = chem::atomInfo(doc);

    p.setPen(Qt::NoPen);
    for (const auto& f : doc.fills) {  // under everything else
        QPolygonF poly;
        for (int i : f.atoms) poly << doc.atoms[i].pos;
        p.setBrush(f.color);
        p.drawPolygon(poly);
    }
    p.setBrush(Qt::NoBrush);
    // A colour the user gave an item wins; otherwise the ink (the theme's on screen, black in exports).
    auto ink = [&](const QColor& own) { return own.isValid() ? own : style.ink; };
    // Aromatic circles: ring bonds drawn single, with a circle inside each ring.
    std::vector<std::vector<int>> circles;
    QSet<int> inCircle;  // bonds drawn single because of a circle
    if (doc.aromaticCircles || !doc.aromaticCircleOverrides.empty()) {
        for (const auto& ring : chem::aromaticRings(doc)) {
            auto ids = ring;
            std::sort(ids.begin(), ids.end());
            const bool override = std::find(doc.aromaticCircleOverrides.begin(),
                                            doc.aromaticCircleOverrides.end(), ids) != doc.aromaticCircleOverrides.end();
            if (doc.aromaticCircles != override) circles.push_back(ring);
        }
        for (const auto& r : circles)
            for (size_t k = 0; k < r.size(); ++k) inCircle.insert(doc.bondBetween(r[k], r[(k + 1) % r.size()]));
    }
    // Crossings: where a later bond (in front; Bring to Front moves one last) crosses an
    // earlier one that shares no atom with it, the earlier one gets a gap.
    // Sweep and prune: bonds sorted by left edge, so only pairs overlapping in x are tested.
    const int nb = int(doc.bonds.size());
    std::vector<QRectF> box(nb);
    for (int i = 0; i < nb; ++i)
        box[i] = QRectF(doc.atoms[doc.bonds[i].a].pos, doc.atoms[doc.bonds[i].b].pos).normalized();
    std::vector<int> byLeft(nb);
    std::iota(byLeft.begin(), byLeft.end(), 0);
    std::sort(byLeft.begin(), byLeft.end(), [&](int i, int j) { return box[i].left() < box[j].left(); });
    std::vector<std::vector<QPointF>> gaps(nb);
    for (int k = 0; k < nb; ++k)
        for (int m = k + 1; m < nb && box[byLeft[m]].left() <= box[byLeft[k]].right(); ++m) {
            const int i = std::min(byLeft[k], byLeft[m]), j = std::max(byLeft[k], byLeft[m]);
            if (box[i].top() > box[j].bottom() || box[j].top() > box[i].bottom()) continue;
            const Bond &x = doc.bonds[i], &y = doc.bonds[j];
            if (x.a == y.a || x.a == y.b || x.b == y.a || x.b == y.b) continue;
            QPointF at;
            if (QLineF(doc.atoms[x.a].pos, doc.atoms[x.b].pos).intersects(QLineF(doc.atoms[y.a].pos, doc.atoms[y.b].pos), &at) ==
                QLineF::BoundedIntersection)
                gaps[i].push_back(at);  // the earlier bond is behind
        }
    for (int bi = 0; bi < int(doc.bonds.size()); ++bi) {
        const Bond& b = doc.bonds[bi];
        p.setPen(QPen(ink(b.color), lineWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        if (inCircle.contains(bi) && b.order == 2 && b.stereo == BondStereo::None) {
            Bond single = b;
            single.order = 1;
            drawBond(p, doc, single, st, degree, labeled, bondsAt, gaps[bi]);
        } else {
            drawBond(p, doc, b, st, degree, labeled, bondsAt, gaps[bi]);
        }
    }
    p.setPen(QPen(style.ink, lineWidth));
    for (const auto& r : circles) {
        QPointF c;
        for (int i : r) c += doc.atoms[i].pos / double(r.size());
        const double apothem = kBondLength / (2 * std::tan(std::numbers::pi / r.size()));
        p.drawEllipse(c, 0.62 * apothem, 0.62 * apothem);
    }

    // H goes on the side away from the bonds; if bonds leave from both
    // left and right, above or below, whichever is clear.
    auto hSide = [&](int i) {
        double dx = 0;
        bool left = false, right = false, up = false;
        for (int nb : neighbors(doc, bondsAt, i)) {
            const QPointF d = doc.atoms[nb].pos - doc.atoms[i].pos;
            dx += d.x();
            left |= d.x() < -0.3 * kBondLength, right |= d.x() > 0.3 * kBondLength, up |= d.y() < -0.3 * kBondLength;
        }
        return left && right ? (up ? HSide::Below : HSide::Above) : dx > 0.1 ? HSide::Left : HSide::Right;
    };
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        p.setPen(QPen(info[i].valenceError ? style.error : ink(a.color), lineWidth));
        if (labeled[i]) {
            drawLabel(p, doc, int(i), doc.hideImplicitH ? 0 : info[i].hydrogens, hSide(int(i)), st);
        } else if (a.charge) {
            QFont sub = labelFont(st, 0.7);
            QString c = QString(a.charge > 0 ? "+" : "−");
            if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
            drawText(p, c, a.pos + QPointF(2, -3), sub);
        }
        if (info[i].valenceError && !labeled[i]) p.drawEllipse(a.pos, 3, 3);
    }
    for (const auto& a : doc.arrows) {
        p.setPen(QPen(ink(a.color), lineWidth, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
        drawArrow(p, a);
    }
    for (const auto& t : doc.texts) p.fillPath(textPath(t, st), ink(t.color));
    // Atom-map numbers (:n) and, when shown, indices from 1: small, in the atom's widest gap,
    // ties (a ring atom's equal gaps) going to the side facing out of the drawing.
    QPointF middle;
    for (const auto& a : doc.atoms) middle += a.pos / double(doc.atoms.size());
    auto numberDirection = [&](int i) {
        const QPointF p = doc.atoms[i].pos;
        std::vector<double> ang;
        for (int nb : neighbors(doc, bondsAt, i)) ang.push_back(std::atan2(doc.atoms[nb].pos.y() - p.y(), doc.atoms[nb].pos.x() - p.x()));
        if (ang.size() < 2) return doc.awayDirection(i);
        std::sort(ang.begin(), ang.end());
        QPointF best;
        double bestGap = -1, bestOut = 0;
        for (size_t k = 0; k < ang.size(); ++k) {
            const double next = k + 1 < ang.size() ? ang[k + 1] : ang[0] + 2 * std::numbers::pi, mid = (ang[k] + next) / 2;
            const QPointF d(std::cos(mid), std::sin(mid));
            const double out = QPointF::dotProduct(d, p - middle);
            if (next - ang[k] > bestGap + 0.1 || (next - ang[k] > bestGap - 0.1 && out > bestOut))
                best = d, bestGap = std::max(bestGap, next - ang[k]), bestOut = out;
        }
        return best;
    };
    // Electron marks and δ±, each in the widest gap left by the bonds, the H and the charge.
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const Atom& a = doc.atoms[i];
        const int marks = a.lonePairs + a.radicals + (a.partial ? 1 : 0);
        if (!marks) continue;
        std::vector<double> taken;
        for (int nb : neighbors(doc, bondsAt, int(i)))
            taken.push_back(std::atan2(doc.atoms[nb].pos.y() - a.pos.y(), doc.atoms[nb].pos.x() - a.pos.x()));
        if (labeled[i] && info[i].hydrogens && !doc.hideImplicitH) {
            const HSide s = hSide(int(i));
            taken.push_back(s == HSide::Right ? 0 : s == HSide::Left ? std::numbers::pi : s == HSide::Below ? std::numbers::pi / 2 : -std::numbers::pi / 2);
        }
        if (a.charge) taken.push_back(-std::numbers::pi / 4);  // up and to the right
        std::vector<QPointF> dirs;
        if (taken.empty()) dirs.push_back({0, -1}), taken.push_back(-std::numbers::pi / 2);  // a lone atom: on top first
        std::vector<std::pair<double, double>> gaps;  // (start, width)
        std::sort(taken.begin(), taken.end());
        for (size_t k = 0; k < taken.size(); ++k) {
            const double next = k + 1 < taken.size() ? taken[k + 1] : taken[0] + 2 * std::numbers::pi;
            gaps.push_back({taken[k], next - taken[k]});
        }
        while (int(dirs.size()) < marks) {
            auto widest = std::max_element(gaps.begin(), gaps.end(), [](auto x, auto y) { return x.second < y.second; });
            const auto [start, width] = *widest;
            const double mid = start + width / 2;
            dirs.push_back({std::cos(mid), std::sin(mid)});
            *widest = {start, mid - start};
            gaps.push_back({mid, start + width - mid});
        }
        const double r = (labeled[i] ? 0.66 : 0.3) * kBondLength, dot = 0.9 * lineWidth + 0.4;
        p.setPen(Qt::NoPen);
        p.setBrush(ink(a.color));
        size_t k = 0;
        for (int n = 0; n < a.radicals; ++n, ++k) p.drawEllipse(a.pos + dirs[k] * r, dot, dot);
        for (int n = 0; n < a.lonePairs; ++n, ++k)
            for (double side : {-1.0, 1.0}) p.drawEllipse(a.pos + dirs[k] * r + perp(dirs[k]) * side * 1.6, dot, dot);
        p.setBrush(Qt::NoBrush);
        if (a.partial) {
            const QFont f = labelFont(st, 0.7);
            const QString s = a.partial > 0 ? "δ+" : "δ−";
            QFontMetricsF fm(f);
            const QPointF at = a.pos + dirs[k] * (r + 2.5);
            p.setPen(QPen(ink(a.color), lineWidth));
            drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
        }
    }
    // Brackets: a pair of [ ] or ( ) just outside their atoms, the label at the bottom right.
    for (const Bracket& b : doc.brackets) {
        QPolygonF pts;
        for (int i : b.atoms) pts << doc.atoms[i].pos;
        const QRectF box = pts.boundingRect().adjusted(-0.55 * kBondLength, -0.6 * kBondLength, 0.55 * kBondLength, 0.6 * kBondLength);
        const double lip = 0.2 * kBondLength;
        p.setPen(QPen(style.ink, lineWidth));
        p.setBrush(Qt::NoBrush);
        for (double side : {-1.0, 1.0}) {
            const double x = side < 0 ? box.left() : box.right();
            QPainterPath path;
            if (b.square) {
                path.moveTo(x - side * lip, box.top());
                path.lineTo(x, box.top());
                path.lineTo(x, box.bottom());
                path.lineTo(x - side * lip, box.bottom());
            } else {
                path.moveTo(x - side * lip, box.top());
                path.quadTo(QPointF(x + side * lip * 0.6, box.center().y()), QPointF(x - side * lip, box.bottom()));
            }
            p.drawPath(path);
        }
        if (!b.label.isEmpty()) drawText(p, b.label, QPointF(box.right() + 1.5, box.bottom() + 1), labelFont(st, 0.7));
    }
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const Atom& a = doc.atoms[i];
        QStringList parts;
        if (doc.showAtomNumbers) parts << QString::number(i + 1);
        if (a.map) parts << ":" + QString::number(a.map);
        if (parts.isEmpty()) continue;
        const QFont f = labelFont(st, 0.6);
        const QString s = parts.join(' ');
        QFontMetricsF fm(f);
        const QPointF at = a.pos + numberDirection(int(i)) * ((labeled[i] ? 0.8 : 0.4) * kBondLength);
        p.setPen(QPen(ink(a.color), lineWidth));
        drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
    }
    if (doc.showStereo) {  // small italic (R)/(E), clear of the atom's bonds or the double bond's second line
        QFont f = labelFont(st, 0.7);
        f.setItalic(true);
        const double off = 0.55 * kBondLength;
        for (const auto& l : chem::stereoLabels(doc)) {
            QPointF at;
            if (l.atom >= 0) {
                at = doc.atoms[l.atom].pos + doc.awayDirection(l.atom) * off;
            } else {
                const Bond& b = doc.bonds[l.bond];
                const QPointF a = doc.atoms[b.a].pos, e = doc.atoms[b.b].pos, n = perp(unit(e - a));
                at = (a + e) / 2 - n * (doubleBondSide(doc, b, bondsAt) > 0 ? 1 : -1) * (0.45 * kBondLength);
            }
            const QString s = "(" + l.text + ")";
            QFontMetricsF fm(f);
            p.setPen(QPen(style.ink, lineWidth));
            drawText(p, s, at - QPointF(fm.horizontalAdvance(s) / 2, -fm.capHeight() / 2), f);
        }
    }
    p.restore();
}

QRectF documentBounds(const Document& doc) {
    const DrawingStyle& st = drawingStyle(doc.style);
    const double fs = st.fontSize;
    if (doc.empty()) return {};
    // Not QRectF::united: it ignores zero-size rects.
    double inf = std::numeric_limits<double>::infinity();
    QPointF lo(inf, inf), hi(-inf, -inf);
    auto grow = [&](QRectF r) {
        lo = {std::min(lo.x(), r.left()), std::min(lo.y(), r.top())};
        hi = {std::max(hi.x(), r.right()), std::max(hi.y(), r.bottom())};
    };
    std::vector<int> degree(doc.atoms.size(), 0);
    for (const auto& b : doc.bonds) ++degree[b.a], ++degree[b.b];
    for (size_t i = 0; i < doc.atoms.size(); ++i) {  // room for labels, which can run either way
        const Atom& a = doc.atoms[i];
        double w = fs * std::max(1.5, 0.7 * a.label.size());
        // A labelled atom with bonds on both sides may stack its H above or below.
        double h = fs * (degree[i] >= 2 && hasLabel(doc, int(i), degree) ? 2.2 : 1);
        grow(QRectF(a.pos, a.pos).adjusted(-w, -h, w, h));
    }
    for (const auto& a : doc.arrows) grow(arrowPath(a).boundingRect().adjusted(-4, -4, 4, 4));
    for (const auto& b : doc.brackets) {  // the brackets and their label
        QPolygonF pts;
        for (int i : b.atoms) pts << doc.atoms[i].pos;
        grow(pts.boundingRect().adjusted(-0.7 * kBondLength, -0.7 * kBondLength, 1.2 * kBondLength, 0.9 * kBondLength));
    }
    for (const auto& t : doc.texts) grow(textPath(t, st).boundingRect().adjusted(-2, -2, 2, 2));
    return QRectF(lo, hi);
}

const std::vector<PageSize>& pageSizes() {
    static const std::vector<PageSize> p{
        {"A4", {595.3, 841.9}, 72},
        {"US Letter", {612, 792}, 72},
        {"ACS single column", {240, 684}},  // 3.33 in × 9.5 in
        {"ACS double column", {504, 684}},  // 7 in
        {"RSC single column", {235.3, 660.5}},  // 8.3 cm × 23.3 cm
        {"RSC double column", {484.7, 660.5}},  // 17.1 cm
    };
    return p;
}

QRectF pageRect(const Document& doc) {
    for (const auto& p : pageSizes())
        if (p.name == doc.page) return QRectF(doc.pageOrigin, p.size / exportScale(doc));
    return {};
}

double exportScale(const Document& doc) { return drawingStyle(doc.style).bondLength / kBondLength; }

// The exported area (model units) and the points per model unit.
static std::pair<QRectF, double> exportFrame(const Document& doc, const ExportOptions& o) {
    const double s = exportScale(doc) * o.scale, m = o.margin / s;
    if (const QRectF page = pageRect(doc); !page.isEmpty()) return {page, s};  // the whole page, as laid out
    return {documentBounds(doc).adjusted(-m, -m, m, m), s};
}

// Paints the frame at `s` output units per point, background first.
static void paintFrame(QPainter& p, const Document& doc, const ExportOptions& o, double perPoint) {
    const auto [r, s] = exportFrame(doc, o);
    p.scale(s * perPoint, s * perPoint);
    p.translate(-r.topLeft());
    if (o.background.alpha()) p.fillRect(r, o.background);
    paintDocument(p, doc);
}

QImage renderImage(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QImage img((r.size() * s * o.dpi / 72.0).toSize().expandedTo({1, 1}), QImage::Format_ARGB32_Premultiplied);
    img.setDotsPerMeterX(int(o.dpi / 0.0254));
    img.setDotsPerMeterY(int(o.dpi / 0.0254));
    img.fill(Qt::transparent);
    QPainter p(&img);
    paintFrame(p, doc, o, o.dpi / 72.0);
    p.end();
    img.setText("penzene", QString::fromUtf8(doc.toJson()));  // reopens as an editable drawing
    return img;
}

QByteArray renderSvg(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QBuffer buf;
    QSvgGenerator gen;
    gen.setOutputDevice(&buf);
    gen.setSize((r.size() * s).toSize());
    gen.setViewBox(QRectF(QPointF(), r.size() * s));
    gen.setResolution(72);  // 1 unit == 1 pt
    gen.setTitle("Penzene");
    QPainter p(&gen);
    paintFrame(p, doc, o, 1);
    p.end();
    // The editable drawing rides along (base64: nothing in it can break the XML).
    QByteArray svg = buf.data();
    const qsizetype open = svg.indexOf('>', svg.indexOf("<svg"));
    if (open > 0) svg.insert(open + 1, "\n<metadata id=\"penzene\">" + doc.toJson().toBase64() + "</metadata>");
    return svg;
}

QByteArray renderPdf(const Document& doc, const ExportOptions& o) {
    const auto [r, s] = exportFrame(doc, o);
    QBuffer buf;
    buf.open(QIODevice::WriteOnly);
    {
        QPdfWriter pdf(&buf);
        pdf.setResolution(72);
        pdf.setPageSize(QPageSize(r.size() * s, QPageSize::Point));
        pdf.setPageMargins({});
        pdf.setCreator("Penzene");
        // The editable drawing rides along as attachments: ours in full, and a MOL for anything else.
        pdf.addFileAttachment("drawing.penz", doc.toJson(), "application/x-penzene");
        if (!doc.atoms.empty())
            if (const std::string mol = chem::toMolBlock(doc); !mol.empty())
                pdf.addFileAttachment("structure.mol", QByteArray::fromStdString(mol), "chemical/x-mdl-molfile");
        QPainter p(&pdf);
        paintFrame(p, doc, o, 1);
    }
    return buf.data();
}

bool exportDocument(const Document& doc, const QString& path, const ExportOptions& o) {
    if (doc.empty()) return false;
    const QString ext = QFileInfo(path).suffix().toLower();
    if (ext == "png") return renderImage(doc, o).save(path);
    if (ext == "svg") {
        QFile f(path);
        return f.open(QIODevice::WriteOnly) && f.write(renderSvg(doc, o)) > 0;
    }
    if (ext == "pdf") {
        QFile f(path);
        return f.open(QIODevice::WriteOnly) && f.write(renderPdf(doc, o)) > 0;
    }
    return false;
}

const std::vector<Theme>& themes() {
    auto cat = [](const char* name, bool dark, const char* base, const char* mantle, const char* surface,
                  const char* text, const char* red, const char* blue, const char* green) {
        return Theme{name, dark, QColor(base), QColor(text), QColor(red), QColor(blue), QColor(green),
                     QColor(mantle), QColor(surface), QColor(text)};
    };
    static const std::vector<Theme> t{
        {"System"},
        {"Light", false, QColor("#FBF8F1"), QColor("#2C2C2A"), QColor(220, 40, 40),
         QColor("#0F6E56"), QColor("#0F6E56"), QColor("#FBF8F1"), QColor("#FFFFFF"), QColor("#2C2C2A")},
        {"Dark", true, QColor("#22211F"), QColor("#F1EFE8"), QColor(255, 105, 97),
         QColor("#5DCAA5"), QColor("#5DCAA5"), QColor("#22211F"), QColor("#2C2C2A"), QColor("#F1EFE8")},
        // https://catppuccin.com/palette: base, mantle, surface0, text, red, blue, green (hotspot)
        // Latte's teal, not its green, for the hotspot: the green is under 3:1 on its base.
        cat("Catppuccin Latte", false, "#eff1f5", "#e6e9ef", "#ccd0da", "#4c4f69", "#d20f39", "#1e66f5", "#179299"),
        cat("Catppuccin Frappé", true, "#303446", "#292c3c", "#414559", "#c6d0f5", "#e78284", "#8caaee", "#a6d189"),
        cat("Catppuccin Macchiato", true, "#24273a", "#1e2030", "#363a4f", "#cad3f5", "#ed8796", "#8aadf4", "#a6da95"),
        cat("Catppuccin Mocha", true, "#1e1e2e", "#181825", "#313244", "#cdd6f4", "#f38ba8", "#89b4fa", "#a6e3a1"),
    };
    return t;
}

const Theme& theme(const QString& name) {
    for (const auto& t : themes())
        if (t.name == name) return t;
    return themes()[0];
}

Chrome chrome(const Theme& t) {
    if (t.name == "Light") return {QColor("#E4E1D6"), QColor("#5F5E5A"), QColor("#E1F5EE")};
    if (t.name == "Dark") return {QColor("#444441"), QColor("#B4B2A9"), QColor("#0B3B30")};
    // The Catppuccin themes keep their own colours.
    return {t.surface.lighter(125), t.text, t.dark ? t.surface.lighter(145) : t.surface.darker(110)};
}
