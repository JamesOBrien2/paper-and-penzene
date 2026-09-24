#include "Render.h"
#include "Chem.h"
#include "Geometry.h"

#include <QBuffer>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QPdfWriter>
#include <QSvgGenerator>
#include <algorithm>
#include <limits>

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
static bool hasLabel(const Document& doc, int i, const std::vector<int>& degree) {
    return doc.atoms[i].z != 6 || degree[i] == 0 || !doc.atoms[i].label.isEmpty();
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

static void drawLabel(QPainter& p, const Document& doc, int i, int hydrogens, bool hLeft, const DrawingStyle& st) {
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
        double hx = hLeft ? x - hw - nw : right;
        drawText(p, "H", {hx, base}, f);
        if (!n.isEmpty()) drawText(p, n, {hx + hw, base + fm.capHeight() * 0.35}, sub);
        if (!hLeft) right += hw + nw;
    }
    if (a.charge) {
        QString c = QString(a.charge > 0 ? "+" : "−");
        if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
        drawText(p, c, {right, base - fm.capHeight() * 0.7}, sub);
    }
}

static void drawBond(QPainter& p, const Document& doc, const Bond& b, const DrawingStyle& st, const std::vector<int>& degree,
                     const std::vector<bool>& labeled) {
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
    // Bold styles the main line, Dashed the other one (or the only one).
    const QPen pen = p.pen();
    auto line = [&](QPointF x, QPointF y, bool main) {
        QPen q = pen;
        if (b.stereo == BondStereo::Bold && main) q.setWidthF(st.boldWidth), q.setCapStyle(Qt::FlatCap);
        if (b.stereo == BondStereo::Dashed && (!main || b.order == 1)) q.setDashPattern({2.5, 2.5});
        p.setPen(q);
        p.drawLine(x, y);
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
        const int side = doubleBondSide(doc, b);
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

// Filled head at `tip` pointing along `dir`; `sides` +1/-1 for a half head.
static void drawHead(QPainter& p, QPointF tip, QPointF dir, int sides = 0) {
    QPointF d = unit(dir), n = perp(d), base = tip - d * kHeadLength;
    QPolygonF head{tip, base + n * (sides >= 0 ? kHeadWidth : 0), tip - d * (kHeadLength * 0.8),
                   base - n * (sides <= 0 ? kHeadWidth : 0)};
    p.setBrush(p.pen().color());
    p.drawPolygon(head);
    p.setBrush(Qt::NoBrush);
}

static void drawArrow(QPainter& p, const Arrow& a) {
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
    std::vector<QPointF> pts = arrowPoints(a);
    const QPointF endDir = pts.back() - pts[pts.size() - 2], startDir = pts.front() - pts[1];
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
    QPen pen(style.ink, lineWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    p.setPen(pen);
    for (const auto& b : doc.bonds) drawBond(p, doc, b, st, degree, labeled);

    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        p.setPen(QPen(info[i].valenceError ? style.error : style.ink, lineWidth));
        if (labeled[i]) {
            // H goes on the side away from the bonds.
            double dx = 0;
            for (int nb : doc.neighbors(int(i))) dx += doc.atoms[nb].pos.x() - a.pos.x();
            drawLabel(p, doc, int(i), info[i].hydrogens, dx > 0.1, st);
        } else if (a.charge) {
            QFont sub = labelFont(st, 0.7);
            QString c = QString(a.charge > 0 ? "+" : "−");
            if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
            drawText(p, c, a.pos + QPointF(2, -3), sub);
        }
        if (info[i].valenceError && !labeled[i]) p.drawEllipse(a.pos, 3, 3);
    }
    p.setPen(QPen(style.ink, lineWidth, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
    for (const auto& a : doc.arrows) drawArrow(p, a);
    for (const auto& t : doc.texts) p.fillPath(textPath(t, st), style.ink);
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
    for (const auto& a : doc.atoms) {  // room for labels, which can run either way
        double w = fs * std::max(1.5, 0.7 * a.label.size());
        grow(QRectF(a.pos, a.pos).adjusted(-w, -fs, w, fs));
    }
    for (const auto& a : doc.arrows) grow(arrowPath(a).boundingRect().adjusted(-4, -4, 4, 4));
    for (const auto& t : doc.texts) grow(textPath(t, st).boundingRect().adjusted(-2, -2, 2, 2));
    return QRectF(lo, hi);
}

double exportScale(const Document& doc) { return drawingStyle(doc.style).bondLength / kBondLength; }

QImage renderImage(const Document& doc, double dpi) {
    QRectF r = documentBounds(doc);
    double s = dpi / 72.0 * exportScale(doc);  // model units to pixels
    QImage img((r.size() * s).toSize().expandedTo({1, 1}), QImage::Format_ARGB32_Premultiplied);
    img.setDotsPerMeterX(int(dpi / 0.0254));
    img.setDotsPerMeterY(int(dpi / 0.0254));
    img.fill(Qt::transparent);
    QPainter p(&img);
    p.scale(s, s);
    p.translate(-r.topLeft());
    paintDocument(p, doc);
    return img;
}

QByteArray renderSvg(const Document& doc) {
    QRectF r = documentBounds(doc);
    const double s = exportScale(doc);
    QBuffer buf;
    QSvgGenerator gen;
    gen.setOutputDevice(&buf);
    gen.setSize((r.size() * s).toSize());
    gen.setViewBox(QRectF(QPointF(), r.size() * s));
    gen.setResolution(72);  // 1 unit == 1 pt
    gen.setTitle("Penzene");
    QPainter p(&gen);
    p.scale(s, s);
    p.translate(-r.topLeft());
    paintDocument(p, doc);
    p.end();
    return buf.data();
}

bool exportDocument(const Document& doc, const QString& path) {
    QRectF r = documentBounds(doc);
    if (doc.empty()) return false;
    const QString ext = QFileInfo(path).suffix().toLower();
    if (ext == "png") return renderImage(doc).save(path);
    if (ext == "svg") {
        QFile f(path);
        return f.open(QIODevice::WriteOnly) && f.write(renderSvg(doc)) > 0;
    }
    if (ext == "pdf") {
        QPdfWriter pdf(path);
        pdf.setResolution(72);
        const double s = exportScale(doc);
        pdf.setPageSize(QPageSize(r.size() * s, QPageSize::Point));
        pdf.setPageMargins({});
        pdf.setCreator("Penzene");
        QPainter p(&pdf);
        p.scale(s, s);
        p.translate(-r.topLeft());
        paintDocument(p, doc);
        return true;
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
        {"Light"},
        {"Dark", true, QColor(0x1e, 0x1e, 0x1e), QColor(0xe6, 0xe6, 0xe6), QColor(255, 105, 97),
         QColor(90, 160, 255), QColor(80, 200, 120)},
        // https://catppuccin.com/palette: base, mantle, surface0, text, red, blue, green
        cat("Catppuccin Latte", false, "#eff1f5", "#e6e9ef", "#ccd0da", "#4c4f69", "#d20f39", "#1e66f5", "#40a02b"),
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

