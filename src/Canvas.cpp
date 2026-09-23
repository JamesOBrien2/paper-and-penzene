#include "Canvas.h"
#include "Chem.h"

#include <QBuffer>
#include <QFile>
#include <QFileInfo>
#include <QImage>
#include <QInputDialog>
#include <QKeyEvent>
#include <QPdfWriter>
#include <QSvgGenerator>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPicture>
#include <QPlainTextEdit>
#include <QRegularExpression>
#include <QScrollBar>
#include <QUndoStack>
#include <QtMath>
#include <algorithm>
#include <functional>
#include <limits>

// Presets. Values come from the ChemDraw stationery (.cds) of the same name;
// wedge width, hash spacing and label gap keep ACS's proportions to ours.
// ponytail: bond length stays 14.4 pt for every preset (JDP's is 14.17, 1.6% off);
// a per-style bond length needs the drawing tools to read it too.
constexpr double kRsc = kBondLength / 12.2;

const std::vector<DrawingStyle>& drawingStyles() {
    static const std::vector<DrawingStyle> styles{
        {"ACS 1996", 0.6, 2.0, 4.5, 2.2, 0.18, 5.5, "Arial", QFont::Normal, 10},
        // JDPReport.cds: line 0.879, bold 1.814, hash 1.814, margin 1.162, IBM Plex Sans Light 10 pt.
        {"JDP", 0.879, 1.814, 4.5 * 1.814 / 2.0, 2.2 * 1.814 / 2.5, 0.18, 5.5 - 1.6 + 1.162, "IBM Plex Sans",
         QFont::Light, 10},
        // RSC (1 Column).cds (2 Column only differs in page size): bond 12.2, line 0.449, bold 1.602,
        // hash 1.75, margin 1.25, spacing 20%, Helvetica 7 pt; scaled by 14.4/12.2 to our bond length.
        {"RSC", 0.449 * kRsc, 1.602 * kRsc, 4.5 * 1.602 / 2.0, 2.2 * 1.75 / 2.5, 0.20, 3.9 * 0.826 + 1.25 * kRsc,
         "Helvetica", QFont::Normal, 7 * kRsc},
    };
    return styles;
}

const DrawingStyle& drawingStyle(const QString& name) {
    for (const auto& s : drawingStyles())
        if (s.name == name) return s;
    return drawingStyles()[0];
}
constexpr double kMergeRadius = 0.3 * kBondLength;

static double len(QPointF v) { return std::hypot(v.x(), v.y()); }
static QPointF unit(QPointF v) { double l = len(v); return l > 1e-9 ? v / l : QPointF(1, 0); }
static QPointF perp(QPointF v) { return {-v.y(), v.x()}; }
static double cross(QPointF a, QPointF b) { return a.x() * b.y() - a.y() * b.x(); }
static QPointF dirAt(double deg) { return {std::cos(qDegreesToRadians(deg)), std::sin(qDegreesToRadians(deg))}; }

static bool hasLabel(const Document& doc, int i, const std::vector<int>& degree) {
    return doc.atoms[i].z != 6 || degree[i] == 0 || !doc.atoms[i].label.isEmpty();
}

// sp centre: a triple bond, or two double bonds (allene). Its bonds are collinear.
static bool isSp(const Document& doc, int atom) {
    int doubles = 0;
    for (const auto& b : doc.bonds)
        if (b.a == atom || b.b == atom) {
            if (b.order == 3) return true;
            doubles += b.order == 2;
        }
    return doubles >= 2;
}

// ---------------------------------------------------------------- rendering

static QFont labelFont(const DrawingStyle& s, double scale = 1) {
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
        double side = 0;
        for (int end : {b.a, b.b})
            for (int nb : doc.neighbors(end))
                if (nb != b.a && nb != b.b) side += cross(d, doc.atoms[nb].pos - pa) > 0 ? 1 : -1;
        // Neighbours on opposite sides (trans chain) tie at 0: still offset, or
        // both lines would cross into the adjoining single bonds.
        // Also centred at an sp centre, so cumulated C=C=C lines meet.
        bool centred = degree[b.a] == 1 || degree[b.b] == 1 || isSp(doc, b.a) || isSp(doc, b.b);
        if (b.position == BondPosition::Centre) centred = true;
        else if (b.position != BondPosition::Auto) centred = false, side = b.position == BondPosition::Right ? 1 : -1;
        if (centred) {
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
static QPointF control(const Arrow& a) { return (a.from + a.to) / 2 - perp(unit(a.to - a.from)) * (2 * a.bend); }

QPainterPath arrowPath(const Arrow& a) {
    QPainterPath path(a.from);
    if (a.bend) path.quadTo(control(a), a.to);
    else path.lineTo(a.to);
    return path;
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
    // Stop the shaft inside the head so it doesn't poke through the tip.
    QPointF endDir = a.bend ? a.to - control(a) : a.to - a.from;
    QPointF startDir = a.bend ? a.from - control(a) : a.from - a.to;
    Arrow shaft = a;
    shaft.to -= unit(endDir) * kHeadLength * 0.7;
    if (a.kind == ArrowKind::Resonance) shaft.from -= unit(startDir) * kHeadLength * 0.7;
    p.drawPath(arrowPath(shaft));
    // Fishhook: the barb sits on the outside of the curve.
    drawHead(p, a.to, endDir, a.kind == ArrowKind::Fishhook ? (a.bend >= 0 ? -1 : 1) : 0);
    if (a.kind == ArrowKind::Resonance) drawHead(p, a.from, startDir);
}

static bool subscripted(const QString& s, int i, bool prevSub) {
    if (!s[i].isDigit() || i == 0) return false;
    QChar c = s[i - 1];
    return prevSub || c.isLetter() || c == ')' || c == ']';
}

// Tab stops every kTabSpaces spaces, on the canvas and in the text dialog alike.
constexpr int kTabSpaces = 8;

// Text as outlines, formula-style subscripts, one line per '\n'. Laid out in
// runs (not per letter) so kerning and spaces match ordinary text.
QPainterPath textPath(const Text& t, const DrawingStyle& st) {
    QFont f = labelFont(st), sub = labelFont(st, 0.7);
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

QImage renderImage(const Document& doc, double dpi) {
    QRectF r = documentBounds(doc);
    double s = dpi / 72.0;  // scene units are points
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
    QBuffer buf;
    QSvgGenerator gen;
    gen.setOutputDevice(&buf);
    gen.setSize(r.size().toSize());
    gen.setViewBox(QRectF(QPointF(), r.size()));
    gen.setResolution(72);  // 1 unit == 1 pt
    gen.setTitle("Penzene");
    QPainter p(&gen);
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
        pdf.setPageSize(QPageSize(r.size(), QPageSize::Point));
        pdf.setPageMargins({});
        pdf.setCreator("Penzene");
        QPainter p(&pdf);
        p.translate(-r.topLeft());
        paintDocument(p, doc);
        return true;
    }
    return false;
}

// ---------------------------------------------------------------- editing helpers

namespace {

// Direction for a new bond from `atom` that avoids existing bonds.
// `newOrder` is the order of the bond about to be added: it makes the atom sp
// (straight on) after a triple bond, or when it cumulates two double bonds.
QPointF freeDirection(const Document& doc, int atom, int newOrder = 1) {
    auto nbs = doc.neighbors(atom);
    QPointF p = doc.atoms[atom].pos;
    if (nbs.empty()) return dirAt(-30);
    if (nbs.size() == 1) {
        QPointF back = unit(doc.atoms[nbs[0]].pos - p);
        int have = doc.bonds[doc.bondBetween(atom, nbs[0])].order;
        if (have == 3 || newOrder == 3 || (have == 2 && newOrder == 2)) return -back;
        double base = qRadiansToDegrees(std::atan2(back.y(), back.x()));
        // Zig-zag: of the two 120° options, take the one farther from everything else.
        QPointF best;
        double bestScore = -1;
        for (double off : {120.0, -120.0}) {
            QPointF cand = p + dirAt(base + off) * kBondLength;
            double score = 1e9;
            for (size_t j = 0; j < doc.atoms.size(); ++j)
                if (int(j) != atom) score = std::min(score, len(doc.atoms[j].pos - cand));
            if (score > bestScore) bestScore = score, best = dirAt(base + off);
        }
        return best;
    }
    return doc.awayDirection(atom);
}

QPointF snapped(QPointF from, QPointF to) {
    QPointF v = to - from;
    double deg = qRadiansToDegrees(std::atan2(v.y(), v.x()));
    return dirAt(std::round(deg / 30) * 30);
}

int atomNear(const Document& doc, QPointF p, double r, int skip = -1) {
    int best = -1;
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        double d = len(doc.atoms[i].pos - p);
        if (int(i) != skip && d < r) r = d, best = int(i);
    }
    return best;
}

// Returns the atom at `p`, creating one if nothing is close enough.
int atomAtOrNew(Document& doc, QPointF p, int z = 6) {
    int i = atomNear(doc, p, kMergeRadius);
    return i >= 0 ? i : doc.addAtom(p, z);
}

void link(Document& doc, int a, int b, int order = 1, BondStereo stereo = BondStereo::None) {
    if (a == b || doc.bondBetween(a, b) >= 0) return;
    doc.bonds.push_back({a, b, order, stereo});
}

bool hasDouble(const Document& doc, int atom) {
    return std::any_of(doc.bonds.begin(), doc.bonds.end(),
                       [&](const Bond& b) { return b.order > 1 && (b.a == atom || b.b == atom); });
}

// Adds a ring through `verts` (merging with existing atoms) and optionally
// alternates double bonds where the atoms are still free for one.
std::vector<int> addRing(Document& doc, const std::vector<QPointF>& verts, bool aromatic) {
    std::vector<int> ids;
    for (QPointF v : verts) ids.push_back(atomAtOrNew(doc, v));
    const int n = int(ids.size());
    for (int k = 0; k < n; ++k) link(doc, ids[k], ids[(k + 1) % n]);
    if (!aromatic) return ids;
    for (int k = 1; k <= n; ++k) {  // start after the (possibly shared) first edge
        int a = ids[k % n], b = ids[(k + 1) % n];
        int bi = doc.bondBetween(a, b);
        if (bi >= 0 && doc.bonds[bi].order == 1 && !hasDouble(doc, a) && !hasDouble(doc, b))
            doc.bonds[bi].order = 2;
    }
    return ids;
}

std::vector<QPointF> polygon(QPointF centre, QPointF firstVertex, int n) {
    std::vector<QPointF> out;
    QPointF r = firstVertex - centre;
    for (int k = 0; k < n; ++k) {
        double t = 2 * M_PI * k / n;
        out.push_back(centre + QPointF(r.x() * std::cos(t) - r.y() * std::sin(t),
                                       r.x() * std::sin(t) + r.y() * std::cos(t)));
    }
    return out;
}

double circumradius(int n) { return kBondLength / (2 * std::sin(M_PI / n)); }

void ringAt(Document& doc, QPointF centre, int n, bool aromatic) {
    addRing(doc, polygon(centre, centre + QPointF(0, -circumradius(n)), n), aromatic);
}

// Ring through the atom, pointing away from its bonds so they bisect the ring's outside angle.
std::vector<int> ringOnAtom(Document& doc, int atom, int n, bool aromatic) {
    QPointF p = doc.atoms[atom].pos;
    return addRing(doc, polygon(p + doc.awayDirection(atom) * circumradius(n), p, n), aromatic);
}

// Ring fused onto the bond, on the side away from the other neighbours.
void ringOnBond(Document& doc, int bond, int n, bool aromatic) {
    const Bond& b = doc.bonds[bond];
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos, d = unit(pb - pa);
    double side = 0;
    for (int end : {b.a, b.b})
        for (int nb : doc.neighbors(end))
            if (nb != b.a && nb != b.b) side += cross(d, doc.atoms[nb].pos - pa);
    double apothem = kBondLength / (2 * std::tan(M_PI / n));
    QPointF centre = (pa + pb) / 2 + perp(d) * (side > 0 ? -apothem : apothem);
    auto verts = polygon(centre, pa, n);
    if (len(verts[1] - pb) > 1) verts = polygon(centre, pb, n);  // wind the right way
    addRing(doc, verts, aromatic);
}

// Chair cyclohexane fused onto the bond, built on template edge `edge` (0 or 1,
// the ChemDraw 9 / 0 keys), mirrored to the side away from the other neighbours.
void chairOnBond(Document& doc, int bond, int edge) {
    // Opposite edges parallel; roughly unit bonds.
    static const QPointF chair[6] = {{0, 0}, {0.95, 0.35}, {1.95, 0.05}, {2.55, 0.75}, {1.6, 0.4}, {0.6, 0.7}};
    const Bond& b = doc.bonds[bond];
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos, d = pb - pa;
    double side = 0;
    for (int end : {b.a, b.b})
        for (int nb : doc.neighbors(end))
            if (nb != b.a && nb != b.b) side += cross(d, doc.atoms[nb].pos - pa);
    QPointF t0 = chair[edge], t1 = chair[edge + 1], td = t1 - t0;
    const double scale = len(d) / len(td);
    std::vector<QPointF> best;
    for (int mirror : {1, -1}) {
        std::vector<QPointF> verts;
        for (int k = 0; k < 6; ++k) {
            QPointF r = chair[(edge + k) % 6] - t0;
            r.setY(r.y() * mirror);
            QPointF td2(td.x(), td.y() * mirror);
            double rr = std::atan2(d.y(), d.x()) - std::atan2(td2.y(), td2.x());
            verts.push_back(pa + QPointF(r.x() * std::cos(rr) - r.y() * std::sin(rr),
                                         r.x() * std::sin(rr) + r.y() * std::cos(rr)) * scale);
        }
        QPointF c;
        for (QPointF v : verts) c += v / 6;
        if (best.empty() || (cross(d, c - pa) > 0) != (side > 0)) best = verts;
    }
    addRing(doc, best, false);
}

// After a bond order change: if an end became an sp centre with two neighbours,
// swing a terminal neighbour into line (Clean handles the general case).
void straightenSp(Document& doc, int bond) {
    for (int e : {doc.bonds[bond].a, doc.bonds[bond].b}) {
        auto nbs = doc.neighbors(e);
        if (nbs.size() != 2 || !isSp(doc, e)) continue;
        for (int k : {0, 1}) {
            int mover = nbs[k], anchor = nbs[1 - k];
            if (doc.neighbors(mover).size() != 1) continue;
            QPointF c = doc.atoms[e].pos;
            doc.atoms[mover].pos = c + unit(c - doc.atoms[anchor].pos) * len(doc.atoms[mover].pos - c);
            break;
        }
    }
}

class Snapshot : public QUndoCommand {
public:
    Snapshot(Canvas* c, Document before, Document after, const QString& text)
        : QUndoCommand(text), c_(c), before_(std::move(before)), after_(std::move(after)) {}
    void undo() override { c_->setDocumentSilently(before_); }
    void redo() override { c_->setDocumentSilently(after_); }

private:
    Canvas* c_;
    Document before_, after_;
};

}  // namespace

// ---------------------------------------------------------------- canvas

Canvas::Canvas(QUndoStack* undo, QWidget* parent) : QGraphicsView(parent), undo_(undo) {
    setScene(new QGraphicsScene(-5000, -5000, 10000, 10000, this));
    setMouseTracking(true);
    setRenderHint(QPainter::Antialiasing);
    setViewportUpdateMode(FullViewportUpdate);
    setTransformationAnchor(AnchorUnderMouse);
    scale(2.5, 2.5);
    centerOn(0, 0);
}

void Canvas::commit(const Document& next, const QString& text) {
    undo_->push(new Snapshot(this, doc_, next, text));
}

void Canvas::setDocumentSilently(const Document& doc) {
    doc_ = doc;
    auto clamp = [](QSet<int>& sel, size_t n) { sel.removeIf([n](int i) { return i >= int(n); }); };
    clamp(selectedAtoms_, doc_.atoms.size());
    clamp(selectedArrows_, doc_.arrows.size());
    clamp(selectedTexts_, doc_.texts.size());
    if (hoverAtom_ >= int(doc_.atoms.size())) hoverAtom_ = -1;
    if (hoverBond_ >= int(doc_.bonds.size())) hoverBond_ = -1;
    refresh();
    emit documentChanged();
}

void Canvas::setSelection(QSet<int> atoms, QSet<int> arrows, QSet<int> texts) {
    selectedAtoms_ = std::move(atoms), selectedArrows_ = std::move(arrows), selectedTexts_ = std::move(texts);
    viewport()->update();
    emit selectionChanged();
}

static QSet<int> range(int from, int to) {
    QSet<int> out;
    for (int i = from; i < to; ++i) out.insert(i);
    return out;
}

void Canvas::selectAll() {
    setSelection(range(0, int(doc_.atoms.size())), range(0, int(doc_.arrows.size())), range(0, int(doc_.texts.size())));
}

// Drops every item not in the given sets.
static Document keepOnly(const Document& doc, const QSet<int>& atoms, const QSet<int>& arrows, const QSet<int>& texts) {
    Document out = doc;
    std::vector<int> drop;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if (!atoms.contains(i)) drop.push_back(i);
    out.removeAtoms(drop);
    out.arrows.clear(), out.texts.clear();
    for (int i = 0; i < int(doc.arrows.size()); ++i)
        if (arrows.contains(i)) out.arrows.push_back(doc.arrows[i]);
    for (int i = 0; i < int(doc.texts.size()); ++i)
        if (texts.contains(i)) out.texts.push_back(doc.texts[i]);
    return out;
}

Document Canvas::selectedSubset() const {
    if (selectedAtoms_.isEmpty() && selectedArrows_.isEmpty() && selectedTexts_.isEmpty()) return doc_;
    return keepOnly(doc_, selectedAtoms_, selectedArrows_, selectedTexts_);
}

void Canvas::deleteSelection() {
    if (selectedAtoms_.isEmpty() && selectedArrows_.isEmpty() && selectedTexts_.isEmpty()) return;
    auto others = [](const QSet<int>& sel, size_t n) { return range(0, int(n)).subtract(sel); };
    Document next = keepOnly(doc_, others(selectedAtoms_, doc_.atoms.size()), others(selectedArrows_, doc_.arrows.size()),
                             others(selectedTexts_, doc_.texts.size()));
    hoverAtom_ = hoverBond_ = -1;
    setSelection({});
    commit(next, tr("Delete"));
}

void Canvas::insert(Document frag, const QString& text) {
    if (frag.empty()) return;
    Document next = doc_;
    next.append(frag, viewCenter() - documentBounds(frag).center());
    commit(next, text);
    setSelection(range(int(doc_.atoms.size() - frag.atoms.size()), int(doc_.atoms.size())),
                 range(int(doc_.arrows.size() - frag.arrows.size()), int(doc_.arrows.size())),
                 range(int(doc_.texts.size() - frag.texts.size()), int(doc_.texts.size())));
}

QPointF Canvas::viewCenter() const { return mapToScene(viewport()->rect().center()); }

void Canvas::zoomBy(double factor) {
    double s = transform().m11() * factor;
    if (s > 0.2 && s < 40) scale(factor, factor);
}

void Canvas::fitToDocument() {
    if (doc_.empty()) return;
    fitInView(documentBounds(doc_).adjusted(-20, -20, 20, 20), Qt::KeepAspectRatio);
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

// Cache the drawing as a QPicture; hover/selection repaints just replay it.
void Canvas::refresh() {
    picture_ = QPicture();
    QPainter p(&picture_);
    paintDocument(p, doc_, {theme_.ink, theme_.error});
    p.end();
    viewport()->update();
}

void Canvas::drawBackground(QPainter* p, const QRectF& rect) {
    p->fillRect(rect, theme_.paper);
    picture_.play(p);
}

void Canvas::drawForeground(QPainter* p, const QRectF&) {
    QColor sel = theme_.accent, hover = theme_.accent, line = theme_.accent;
    sel.setAlpha(90), hover.setAlpha(60);
    p->setRenderHint(QPainter::Antialiasing);
    p->setPen(Qt::NoPen);
    for (const auto& b : doc_.bonds)
        if (selectedAtoms_.contains(b.a) && selectedAtoms_.contains(b.b)) {
            p->setPen(QPen(sel, 3, Qt::SolidLine, Qt::RoundCap));
            p->drawLine(doc_.atoms[b.a].pos, doc_.atoms[b.b].pos);
        }
    p->setPen(Qt::NoPen);
    p->setBrush(sel);
    for (int i : selectedAtoms_) p->drawEllipse(doc_.atoms[i].pos, 4, 4);

    for (int i : selectedArrows_) p->strokePath(arrowPath(doc_.arrows[i]), QPen(sel, 4, Qt::SolidLine, Qt::RoundCap));
    for (int i : selectedTexts_) p->drawRect(textPath(doc_.texts[i], drawingStyle(doc_.style)).boundingRect().adjusted(-1.5, -1.5, 1.5, 1.5));

    p->setBrush(hover);
    if (hoverAtom_ >= 0) {
        p->drawEllipse(doc_.atoms[hoverAtom_].pos, 5, 5);
        p->setBrush(theme_.hotspot);
        p->drawEllipse(doc_.atoms[hoverAtom_].pos, 1.2, 1.2);
    } else if (hoverBond_ >= 0) {
        const auto& b = doc_.bonds[hoverBond_];
        p->setPen(QPen(hover, 5, Qt::SolidLine, Qt::RoundCap));
        p->drawLine(doc_.atoms[b.a].pos, doc_.atoms[b.b].pos);
    }

    p->setBrush(Qt::NoBrush);
    if (drag_ == Drag::Rubber) {
        p->setPen(QPen(line, 0, Qt::DashLine));
        p->drawRect(QRectF(pressPos_, curPos_).normalized());
    } else if (drag_ == Drag::Bond || drag_ == Drag::Chain) {
        p->setPen(QPen(line, 0.8));
        for (size_t k = 1; k < preview_.size(); ++k) p->drawLine(preview_[k - 1], preview_[k]);
    } else if (drag_ == Drag::Arrow) {
        Document preview;
        preview.arrows.push_back(draggedArrow());
        paintDocument(*p, preview, {line});
    }
}

int Canvas::arrowAt(QPointF p) const {
    double tol = 5 / transform().m11() + 1;
    QPainterPathStroker stroker;
    stroker.setWidth(2 * tol);
    for (int i = int(doc_.arrows.size()) - 1; i >= 0; --i)
        if (stroker.createStroke(arrowPath(doc_.arrows[i])).contains(p)) return i;
    return -1;
}

int Canvas::textAt(QPointF p) const {
    for (int i = int(doc_.texts.size()) - 1; i >= 0; --i)
        if (textPath(doc_.texts[i], drawingStyle(doc_.style)).boundingRect().adjusted(-2, -2, 2, 2).contains(p)) return i;
    return -1;
}

// The arrow being dragged out: straight ones snap to 15°, curved ones bow left.
Arrow Canvas::draggedArrow() const {
    Arrow a{pressPos_, curPos_, arrowKind_};
    if (arrowCurved_) {
        a.bend = 0.3 * len(curPos_ - pressPos_);
    } else {
        QPointF v = curPos_ - pressPos_;
        double deg = std::round(qRadiansToDegrees(std::atan2(v.y(), v.x())) / 15) * 15;
        a.to = pressPos_ + dirAt(deg) * len(v);
    }
    return a;
}

int Canvas::atomAt(QPointF p) const {
    return atomNear(doc_, p, 8 / transform().m11() + 2);
}

int Canvas::bondAt(QPointF p) const {
    double tol = 6 / transform().m11() + 1;
    for (size_t i = 0; i < doc_.bonds.size(); ++i) {
        QPointF a = doc_.atoms[doc_.bonds[i].a].pos, b = doc_.atoms[doc_.bonds[i].b].pos;
        QPointF ab = b - a;
        double t = std::clamp(QPointF::dotProduct(p - a, ab) / QPointF::dotProduct(ab, ab), 0.0, 1.0);
        if (len(a + ab * t - p) < tol) return int(i);
    }
    return -1;
}

// Points of the bond or chain the user is dragging out.
std::vector<QPointF> Canvas::dragPath() const {
    QPointF start = pressAtom_ >= 0 ? doc_.atoms[pressAtom_].pos : pressPos_;
    int target = atomNear(doc_, curPos_, kMergeRadius, pressAtom_);
    if (drag_ == Drag::Bond) {
        if (target >= 0) return {start, doc_.atoms[target].pos};
        return {start, start + snapped(start, curPos_) * kBondLength};
    }
    // Chain: zig-zag at ±30° around the drag direction.
    QPointF dir = snapped(start, curPos_), side = perp(dir);
    double step = kBondLength * std::cos(M_PI / 6), rise = kBondLength * std::sin(M_PI / 6);
    int n = std::max(1, int(std::round(QPointF::dotProduct(curPos_ - start, dir) / step)));
    double sign = cross(dir, curPos_ - start) >= 0 ? 1 : -1;
    std::vector<QPointF> pts{start};
    for (int k = 1; k <= n; ++k) pts.push_back(start + dir * (step * k) + side * (k % 2 ? rise * sign : 0));
    return pts;
}

void Canvas::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::MiddleButton) {
        drag_ = Drag::Pan;
        panLast_ = e->pos();
        return;
    }
    if (e->button() != Qt::LeftButton) return;
    pressPos_ = curPos_ = mapToScene(e->pos());
    pressAtom_ = atomAt(pressPos_);
    int bond = pressAtom_ < 0 ? bondAt(pressPos_) : -1;
    beforeDrag_ = doc_;

    switch (tool_) {
    case Tool::Select: {
        const int arrow = pressAtom_ < 0 && bond < 0 ? arrowAt(pressPos_) : -1;
        const int text = pressAtom_ < 0 && bond < 0 && arrow < 0 ? textAt(pressPos_) : -1;
        const bool shift = e->modifiers() & Qt::ShiftModifier;
        if (pressAtom_ >= 0 || bond >= 0 || arrow >= 0 || text >= 0) {
            QSet<int> atoms = pressAtom_ >= 0 ? QSet<int>{pressAtom_}
                              : bond >= 0     ? QSet<int>{doc_.bonds[bond].a, doc_.bonds[bond].b}
                                              : QSet<int>{};
            QSet<int> arrows = arrow >= 0 ? QSet<int>{arrow} : QSet<int>{};
            QSet<int> texts = text >= 0 ? QSet<int>{text} : QSet<int>{};
            bool already = selectedAtoms_.contains(atoms) && selectedArrows_.contains(arrows) &&
                           selectedTexts_.contains(texts);
            if (shift) selectedAtoms_ |= atoms, selectedArrows_ |= arrows, selectedTexts_ |= texts;
            else if (!already) selectedAtoms_ = atoms, selectedArrows_ = arrows, selectedTexts_ = texts;
            drag_ = (e->modifiers() & Qt::AltModifier) ? Drag::Rotate : Drag::Move;
        } else {
            if (!shift) selectedAtoms_.clear(), selectedArrows_.clear(), selectedTexts_.clear();
            drag_ = Drag::Rubber;
        }
        break;
    }
    case Tool::Bond: case Tool::Wedge: case Tool::Hash:
        drag_ = Drag::Bond;
        break;
    case Tool::Chain:
        drag_ = Drag::Chain;
        break;
    case Tool::Arrow:
        drag_ = Drag::Arrow;
        break;
    default:
        drag_ = Drag::None;  // click tools act on release
    }
    viewport()->update();
    emit selectionChanged();
}

void Canvas::mouseMoveEvent(QMouseEvent* e) {
    if (drag_ == Drag::Pan) {
        QPoint d = e->pos() - panLast_;
        panLast_ = e->pos();
        horizontalScrollBar()->setValue(horizontalScrollBar()->value() - d.x());
        verticalScrollBar()->setValue(verticalScrollBar()->value() - d.y());
        return;
    }
    curPos_ = mapToScene(e->pos());
    if (drag_ == Drag::Move || drag_ == Drag::Rotate) {
        Document next = beforeDrag_;
        std::vector<QPointF*> pts;
        for (int i : selectedAtoms_) pts.push_back(&next.atoms[i].pos);
        for (int i : selectedArrows_) pts.push_back(&next.arrows[i].from), pts.push_back(&next.arrows[i].to);
        for (int i : selectedTexts_) pts.push_back(&next.texts[i].pos);
        QPointF c;
        for (QPointF* p : pts) c += *p;
        c /= std::max<double>(1, pts.size());
        double ang = std::atan2(curPos_.y() - c.y(), curPos_.x() - c.x()) -
                     std::atan2(pressPos_.y() - c.y(), pressPos_.x() - c.x());
        for (QPointF* p : pts) {
            if (drag_ == Drag::Move) {
                *p += curPos_ - pressPos_;
            } else {
                QPointF r = *p - c;
                *p = c + QPointF(r.x() * std::cos(ang) - r.y() * std::sin(ang),
                                 r.x() * std::sin(ang) + r.y() * std::cos(ang));
            }
        }
        doc_ = next;
        refresh();
        return;
    }
    if (drag_ == Drag::Bond || drag_ == Drag::Chain) preview_ = dragPath();
    if (drag_ == Drag::None) {
        // The hotspot sticks until the cursor reaches another atom or bond, so
        // hotkeys and arrow keys keep working after the mouse drifts off.
        if (int a = atomAt(curPos_); a >= 0) hoverAtom_ = a, hoverBond_ = -1;
        else if (int b = bondAt(curPos_); b >= 0) hoverBond_ = b, hoverAtom_ = -1;
    }
    viewport()->update();
}

void Canvas::mouseReleaseEvent(QMouseEvent* e) {
    if (drag_ == Drag::Pan) {
        drag_ = Drag::None;
        return;
    }
    if (e->button() != Qt::LeftButton) return;
    curPos_ = mapToScene(e->pos());
    const bool click = len(curPos_ - pressPos_) < 3 / transform().m11();
    const int bond = pressAtom_ < 0 ? bondAt(pressPos_) : -1;
    const Drag drag = std::exchange(drag_, Drag::None);
    preview_.clear();
    Document next = beforeDrag_;
    QString what;

    const BondStereo stereo = tool_ == Tool::Wedge ? BondStereo::Wedge
                              : tool_ == Tool::Hash ? BondStereo::Hash
                                                    : BondStereo::None;
    const int order = stereo == BondStereo::None ? bondOrder_ : 1;

    if (drag == Drag::Move || drag == Drag::Rotate) {
        if (!click) {
            Document moved = doc_;
            doc_ = beforeDrag_;
            commit(moved, drag == Drag::Move ? tr("Move") : tr("Rotate"));
        }
        return;
    } else if (drag == Drag::Rubber) {
        QRectF r = QRectF(pressPos_, curPos_).normalized();
        for (int i = 0; i < int(doc_.atoms.size()); ++i)
            if (r.contains(doc_.atoms[i].pos)) selectedAtoms_.insert(i);
        for (int i = 0; i < int(doc_.arrows.size()); ++i)
            if (r.contains(doc_.arrows[i].from) && r.contains(doc_.arrows[i].to)) selectedArrows_.insert(i);
        for (int i = 0; i < int(doc_.texts.size()); ++i)
            if (r.intersects(textPath(doc_.texts[i], drawingStyle(doc_.style)).boundingRect())) selectedTexts_.insert(i);
    } else if (drag == Drag::Arrow) {
        if (int hit = arrowAt(pressPos_); click && hit >= 0) {  // click an arrow: restyle, or flip a curve
            Arrow& a = next.arrows[hit];
            if (arrowCurved_ && a.bend && a.kind == arrowKind_) a.bend = -a.bend;
            else a.kind = arrowKind_, a.bend = arrowCurved_ ? 0.3 * len(a.to - a.from) : 0;
        } else {
            if (click) curPos_ = pressPos_ + QPointF(3 * kBondLength, 0);  // default length
            next.arrows.push_back(draggedArrow());
        }
        what = tr("Arrow");
    } else if ((drag == Drag::Bond || drag == Drag::Chain) && click) {
        if (bond >= 0) {  // click on a bond: change it in place
            Bond& b = next.bonds[bond];
            if (stereo != BondStereo::None) {
                if (b.stereo == stereo) std::swap(b.a, b.b);  // flip direction
                b.stereo = stereo, b.order = 1;
            } else {
                b.stereo = BondStereo::None;
                b.order = (b.order != order && order > 1) ? order : b.order % 3 + 1;
                straightenSp(next, bond);
            }
            what = tr("Change bond");
        } else {
            int from = pressAtom_ >= 0 ? pressAtom_ : next.addAtom(pressPos_);
            QPointF to = next.atoms[from].pos + freeDirection(next, from, order) * kBondLength;
            link(next, from, atomAtOrNew(next, to), order, stereo);
            what = tr("Add bond");
        }
    } else if (drag == Drag::Bond || drag == Drag::Chain) {
        drag_ = drag;  // dragPath reads it
        auto pts = dragPath();
        drag_ = Drag::None;
        int prev = pressAtom_ >= 0 ? pressAtom_ : next.addAtom(pts[0]);
        for (size_t k = 1; k < pts.size(); ++k) {
            int cur = atomAtOrNew(next, pts[k]);
            link(next, prev, cur, drag == Drag::Bond ? order : 1, drag == Drag::Bond ? stereo : BondStereo::None);
            prev = cur;
        }
        what = drag == Drag::Bond ? tr("Add bond") : tr("Add chain");
    } else if (click) {
        switch (tool_) {
        case Tool::Atom:
            if (pressAtom_ >= 0) next.atoms[pressAtom_].z = element_;
            else if (bond < 0) next.addAtom(pressPos_, element_);
            what = tr("Set atom");
            break;
        case Tool::ChargePlus: case Tool::ChargeMinus:
            if (pressAtom_ >= 0) next.atoms[pressAtom_].charge += tool_ == Tool::ChargePlus ? 1 : -1;
            what = tr("Charge");
            break;
        case Tool::Erase:
            if (pressAtom_ >= 0) next.removeAtoms({pressAtom_});
            else if (bond >= 0) next.bonds.erase(next.bonds.begin() + bond);
            else if (int a = arrowAt(pressPos_); a >= 0) next.arrows.erase(next.arrows.begin() + a);
            else if (int t = textAt(pressPos_); t >= 0) next.texts.erase(next.texts.begin() + t);
            what = tr("Erase");
            break;
        case Tool::Text:
            return editText(textAt(pressPos_), pressPos_);
        case Tool::Fill: {
            // Smallest ring around the click; clicking a ring in the same colour clears it.
            std::vector<int> best;
            double bestArea = std::numeric_limits<double>::infinity();
            for (const auto& ring : chem::rings(next)) {
                QPolygonF poly;
                for (int i : ring) poly << next.atoms[i].pos;
                QRectF r = poly.boundingRect();
                if (poly.containsPoint(pressPos_, Qt::OddEvenFill) && r.width() * r.height() < bestArea)
                    bestArea = r.width() * r.height(), best = ring;
            }
            if (best.empty()) break;
            auto same = [&](const ::Fill& f) {
                return QSet<int>(f.atoms.begin(), f.atoms.end()) == QSet<int>(best.begin(), best.end());
            };
            auto it = std::find_if(next.fills.begin(), next.fills.end(), same);
            if (it == next.fills.end()) next.fills.push_back({best, fillColor_});
            else if (it->color == fillColor_) next.fills.erase(it);
            else it->color = fillColor_;
            what = tr("Ring fill");
            break;
        }
        case Tool::Ring:
            if (bond >= 0) ringOnBond(next, bond, ringSize_, ringAromatic_);
            else if (pressAtom_ >= 0) ringOnAtom(next, pressAtom_, ringSize_, ringAromatic_);
            else ringAt(next, pressPos_, ringSize_, ringAromatic_);
            what = tr("Add ring");
            break;
        default:
            break;
        }
    }
    if (!what.isEmpty() && !(next == beforeDrag_)) commit(next, what);
    viewport()->update();
    emit selectionChanged();
}

void Canvas::mouseDoubleClickEvent(QMouseEvent* e) {
    if (tool_ != Tool::Select) return QGraphicsView::mouseDoubleClickEvent(e);
    int start = atomAt(mapToScene(e->pos()));
    if (int t = textAt(mapToScene(e->pos())); start < 0 && t >= 0) return editText(t);
    if (start < 0) return;
    // Select the whole connected fragment.
    QSet<int> seen{start};
    std::vector<int> stack{start};
    while (!stack.empty()) {
        int i = stack.back();
        stack.pop_back();
        for (int nb : doc_.neighbors(i))
            if (!seen.contains(nb)) seen.insert(nb), stack.push_back(nb);
    }
    setSelection(seen);
}

// ---------------------------------------------------------------- hotkeys
// Follows ChemDraw 21's atom and bond hotkey tables (User Guide, ch. 5).

// Atom label hotkeys (the hotspot atom becomes this element or group).
static QString labelHotkey(const QString& key) {
    static const QHash<QString, QString> k{
        {"c", "C"},   {"n", "N"},    {"w", "N"},     {"o", "O"},    {"q", "O"},   {"s", "S"},
        {"p", "P"},   {"f", "F"},    {"l", "Cl"},    {"C", "Cl"},   {"b", "Br"},  {"i", "I"},
        {"h", "H"},   {"d", "H"},    {"B", "B"},     {"S", "Si"},   {"L", "Li"},  {"m", "Me"},
        {"e", "Et"},  {"A", "Ac"},   {"P", "Ph"},    {"F", "CF3"},  {"N", "NO2"}, {"O", "OMe"},
        {"E", "CO2Me"}, {"Z", "N3"}, {"M", "MgBr"},  {"Q", "Fmoc"}, {"H", "Cbz"}, {"Y", "Boc"},
    };
    return k.value(key);
}

// Element symbol, abbreviation (drawn as its label) or SMILES (drawn out).
bool Canvas::applyLabel(Document& doc, int at, const QString& label) {
    Atom& a = doc.atoms[at];
    // "OH", "NH2": the element; hydrogens are implicit.
    static const QRegularExpression hydride("^([A-Z][a-z]?)H\\d*$");
    QString element = hydride.match(label).hasMatch() ? hydride.match(label).captured(1) : label;
    if (int z = chem::atomicNumber(element.toStdString()); z > 0) {
        a.z = z, a.label.clear();
        return true;
    }
    if (auto head = chem::abbreviationHead(label)) {
        a.z = head->z, a.charge = head->charge, a.label = label;
        return true;
    }
    return chem::attach(doc, at, label.toStdString());
}

namespace {

enum class Site { Primary, Secondary, Tertiary, Aromatic };

bool inRing(const Document& doc, int at) {
    auto nbs = doc.neighbors(at);
    for (int start : nbs) {  // can another neighbour be reached without passing through `at`?
        std::vector<bool> seen(doc.atoms.size());
        seen[at] = seen[start] = true;
        std::vector<int> stack{start};
        while (!stack.empty()) {
            int i = stack.back();
            stack.pop_back();
            for (int nb : doc.neighbors(i)) {
                if (nb == at && i != start) return true;
                if (!seen[nb]) seen[nb] = true, stack.push_back(nb);
            }
        }
    }
    return false;
}

Site site(const Document& doc, int at) {
    size_t deg = doc.neighbors(at).size();
    if (deg <= 1) return Site::Primary;
    if (deg >= 3) return Site::Tertiary;
    return inRing(doc, at) && hasDouble(doc, at) ? Site::Aromatic : Site::Secondary;
}

// The two 120° directions open to an atom with at most one bond: the zig-zag
// ("linear mode") one and the other ("cyclic mode").
std::pair<QPointF, QPointF> openDirections(const Document& doc, int at) {
    QPointF lin = freeDirection(doc, at);
    auto nbs = doc.neighbors(at);
    if (nbs.empty()) return {lin, dirAt(90)};
    QPointF back = unit(doc.atoms[nbs[0]].pos - doc.atoms[at].pos);
    // Reflect `lin` across the bond axis.
    QPointF other = 2 * QPointF::dotProduct(lin, back) * back - lin;
    return {lin, other};
}

QPointF rotated(QPointF v, double deg) {
    double a = qDegreesToRadians(deg);
    return {v.x() * std::cos(a) - v.y() * std::sin(a), v.x() * std::sin(a) + v.y() * std::cos(a)};
}

int sprout(Document& doc, int from, QPointF dir, int order = 1, int z = 6,
           BondStereo stereo = BondStereo::None, double length = kBondLength) {
    int n = atomAtOrNew(doc, doc.atoms[from].pos + dir * length, z);
    link(doc, from, n, order, stereo);
    return n;
}

// "Adds a C-C bond and then ...": tertiary and aromatic sites grow from a new carbon.
int linker(Document& doc, int at) { return sprout(doc, at, doc.awayDirection(at)); }

int farthestFrom(const Document& doc, const std::vector<int>& ids, int from) {
    int best = from;
    double d = -1;
    for (int i : ids)
        if (double l = len(doc.atoms[i].pos - doc.atoms[from].pos); l > d) d = l, best = i;
    return best;
}

constexpr int kUnhandled = -2;

// Atom "sprout" hotkeys. Returns the new hotspot atom, or kUnhandled.
int sproutHotkey(Document& doc, int at, const QString& key) {
    const Site s = site(doc, at);
    const bool needsLinker = s == Site::Tertiary || s == Site::Aromatic;

    if (key == "1") return sprout(doc, at, freeDirection(doc, at));
    if (key == "0") {
        if (s == Site::Primary) return sprout(doc, at, openDirections(doc, at).second);
        return sprout(doc, at, doc.awayDirection(at), 1, 6, BondStereo::None, 1.5 * kBondLength);
    }
    if (key == "2") {  // carbonyl / acetyl
        if (s == Site::Secondary) {
            sprout(doc, at, doc.awayDirection(at), 2, 8);
            return at;
        }
        int c = needsLinker ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, other, 2, 8);
        return sprout(doc, c, lin);
    }
    if (key == "3" || key == "a") {  // phenyl
        int c = s == Site::Primary ? at : linker(doc, at);
        return farthestFrom(doc, ringOnAtom(doc, c, 6, true), c);
    }
    if (key == "4" || key == "5") {  // wedged / hashed methyl
        const BondStereo st = key == "4" ? BondStereo::Wedge : BondStereo::Hash;
        if (s == Site::Secondary || s == Site::Tertiary) {
            sprout(doc, at, doc.awayDirection(at), 1, 6, st);
            return at;
        }
        int c = s == Site::Aromatic ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, other, 1, 6, st);
        return sprout(doc, c, lin);
    }
    static const QHash<QString, int> rings{{"6", 6}, {"7", 5}, {"u", 4}, {"v", 3}};
    if (rings.contains(key)) {  // cycloalkyl; spiro on a secondary carbon
        int c = needsLinker ? linker(doc, at) : at;
        return farthestFrom(doc, ringOnAtom(doc, c, rings[key], false), c);
    }
    if (key == "8") {  // methylidene
        int c = needsLinker ? linker(doc, at) : at;
        QPointF dir = site(doc, c) == Site::Primary ? freeDirection(doc, c, 2) : doc.awayDirection(c);
        return sprout(doc, c, dir, 2);
    }
    if (key == "9") {  // dimethyl / gem-dimethyl / isopropyl
        if (s == Site::Secondary) {
            QPointF away = doc.awayDirection(at);
            sprout(doc, at, rotated(away, 60));
            sprout(doc, at, rotated(away, -60));
            return at;
        }
        int c = needsLinker ? linker(doc, at) : at;
        auto [lin, other] = openDirections(doc, c);
        sprout(doc, c, lin);
        sprout(doc, c, other);
        return c;
    }
    if (key == "z") {  // alkyne, linear
        QPointF dir = freeDirection(doc, at);
        int c1 = sprout(doc, at, dir);
        return sprout(doc, c1, dir, 3);
    }
    if (key == "k") {  // sulfonyl
        QPointF dir = freeDirection(doc, at);
        int sulfur = sprout(doc, at, dir, 1, 16);
        sprout(doc, sulfur, perp(dir), 2, 8);
        sprout(doc, sulfur, -perp(dir), 2, 8);
        return sulfur;
    }
    if (key == "K") {  // t-Bu at 90°
        QPointF dir = freeDirection(doc, at);
        int c = sprout(doc, at, dir);
        sprout(doc, c, dir);
        sprout(doc, c, perp(dir));
        sprout(doc, c, -perp(dir));
        return at;
    }
    return kUnhandled;
}

QPointF arrowDirection(int key) {
    switch (key) {
    case Qt::Key_Left: return {-1, 0};
    case Qt::Key_Right: return {1, 0};
    case Qt::Key_Up: return {0, -1};
    default: return {0, 1};
    }
}

// Candidate (score, index) with the best direction match above `minDot`.
template <class F>
int bestToward(int count, QPointF dir, double minDot, F vectorOf) {
    int best = -1;
    for (int i = 0; i < count; ++i) {
        QPointF v = vectorOf(i);
        if (len(v) < 1e-6) continue;
        if (double d = QPointF::dotProduct(unit(v), dir); d > minDot) minDot = d, best = i;
    }
    return best;
}

}  // namespace

void Canvas::rotateSelection(double degrees) {
    Document next = doc_;
    std::vector<QPointF*> pts;
    for (int i : selectedAtoms_) pts.push_back(&next.atoms[i].pos);
    for (int i : selectedArrows_) pts.push_back(&next.arrows[i].from), pts.push_back(&next.arrows[i].to);
    for (int i : selectedTexts_) pts.push_back(&next.texts[i].pos);
    if (pts.empty()) return;
    QPointF c;
    for (QPointF* p : pts) c += *p;
    c /= double(pts.size());
    for (QPointF* p : pts) *p = c + rotated(*p - c, degrees);
    commit(next, tr("Rotate"));
}

// Arrow keys walk atom -> bond -> atom; with Shift, atom -> atom or bond -> bond.
void Canvas::moveHotspot(QPointF dir, bool jump) {
    const auto& d = doc_;
    auto mid = [&](int b) { return (d.atoms[d.bonds[b].a].pos + d.atoms[d.bonds[b].b].pos) / 2; };
    auto bondsOf = [&](int atom) {
        std::vector<int> out;
        for (int i = 0; i < int(d.bonds.size()); ++i)
            if (d.bonds[i].a == atom || d.bonds[i].b == atom) out.push_back(i);
        return out;
    };
    if (hoverAtom_ >= 0) {
        QPointF p = d.atoms[hoverAtom_].pos;
        if (jump) {
            auto nbs = d.neighbors(hoverAtom_);
            int k = bestToward(int(nbs.size()), dir, 0.3, [&](int i) { return d.atoms[nbs[i]].pos - p; });
            if (k >= 0) hoverAtom_ = nbs[k];
        } else {
            auto bs = bondsOf(hoverAtom_);
            int k = bestToward(int(bs.size()), dir, 0.3, [&](int i) { return mid(bs[i]) - p; });
            if (k >= 0) hoverBond_ = bs[k], hoverAtom_ = -1;
        }
    } else if (hoverBond_ >= 0) {
        const Bond& b = d.bonds[hoverBond_];
        QPointF m = mid(hoverBond_);
        if (jump) {
            std::vector<int> adj;
            for (int end : {b.a, b.b})
                for (int bi : bondsOf(end))
                    if (bi != hoverBond_) adj.push_back(bi);
            int k = bestToward(int(adj.size()), dir, 0.3, [&](int i) { return mid(adj[i]) - m; });
            if (k >= 0) hoverBond_ = adj[k];
        } else {
            int ends[] = {b.a, b.b};
            int k = bestToward(2, dir, 0.1, [&](int i) { return d.atoms[ends[i]].pos - m; });
            if (k >= 0) hoverAtom_ = ends[k], hoverBond_ = -1;
        }
    }
    viewport()->update();
}

void Canvas::editLabel(int at) {
    bool ok = false;
    QString label = QInputDialog::getText(this, tr("Atom label"),
                                          tr("Element, group (OMe, CF3, Ph, Boc…) or SMILES:"), QLineEdit::Normal,
                                          doc_.atoms[at].label.isEmpty()
                                              ? QString::fromStdString(chem::symbol(doc_.atoms[at].z))
                                              : doc_.atoms[at].label,
                                          &ok)
                        .trimmed();
    Document next = doc_;
    if (ok && !label.isEmpty() && applyLabel(next, at, label)) commit(next, tr("Edit label"));
}

void Canvas::expandAbbreviations() {
    Document next = doc_;
    for (int i = 0; i < int(doc_.atoms.size()); ++i) {
        bool wanted = selectedAtoms_.isEmpty() ? (hoverAtom_ < 0 || hoverAtom_ == i) : selectedAtoms_.contains(i);
        if (wanted && !doc_.atoms[i].label.isEmpty()) chem::attach(next, i, doc_.atoms[i].label.toStdString());
    }
    if (!(next == doc_)) commit(next, tr("Expand"));
}

// Copies the selection to the far side of the next arrow in `dir` (ChemDraw's
// Ctrl+arrow), or just past the selection when there is no arrow that way.
void Canvas::duplicateSelection(QPointF dir) {
    Document copy = selectedSubset();
    if (selectedAtoms_.isEmpty() && selectedArrows_.isEmpty() && selectedTexts_.isEmpty()) return;
    const QRectF box = documentBounds(copy);
    const QPointF c = box.center();
    const double half = std::abs(QPointF::dotProduct(QPointF(box.width(), box.height()) / 2, dir));
    double shift = 2 * half + 2 * kBondLength;
    double nearest = std::numeric_limits<double>::infinity();
    for (const auto& a : doc_.arrows) {
        double from = QPointF::dotProduct(a.from - c, dir), to = QPointF::dotProduct(a.to - c, dir);
        double lo = std::min(from, to), hi = std::max(from, to);
        if (lo > half - 1 && lo < nearest) nearest = lo, shift = hi + (lo - half) + half;
    }
    Document next = doc_;
    next.append(copy, dir * shift);
    commit(next, tr("Duplicate"));
    setSelection(range(int(doc_.atoms.size() - copy.atoms.size()), int(doc_.atoms.size())),
                 range(int(doc_.arrows.size() - copy.arrows.size()), int(doc_.arrows.size())),
                 range(int(doc_.texts.size() - copy.texts.size()), int(doc_.texts.size())));
}

void Canvas::editText(int i, QPointF pos) {
    // The editor uses the canvas font and tab stops, so spacing looks the same on both.
    QInputDialog dialog(this);
    dialog.setWindowTitle(tr("Text"));
    dialog.setLabelText(tr("Text (digits after letters become subscripts):"));
    dialog.setOption(QInputDialog::UsePlainTextEditForTextInput);
    dialog.setTextValue(i >= 0 ? doc_.texts[i].text : QString());
    if (auto* edit = dialog.findChild<QPlainTextEdit*>()) {
        QFont f = labelFont(drawingStyle(doc_.style));
        f.setPixelSize(16);
        edit->setFont(f);
        edit->setTabStopDistance(kTabSpaces * QFontMetricsF(f).horizontalAdvance(' '));
    }
    if (dialog.exec() != QDialog::Accepted) return;
    // Keep leading spaces and tabs; they are deliberate indentation.
    QString s = dialog.textValue();
    s.remove(QRegularExpression("\\s+$"));
    Document next = doc_;
    if (i < 0 && !s.isEmpty()) next.texts.push_back({pos, s});
    else if (i >= 0 && s.isEmpty()) next.texts.erase(next.texts.begin() + i);
    else if (i >= 0) next.texts[i].text = s;
    if (!(next == doc_)) commit(next, tr("Text"));
}

void Canvas::keyPressEvent(QKeyEvent* e) {
    const int key = e->key();
    const bool arrow = key == Qt::Key_Left || key == Qt::Key_Right || key == Qt::Key_Up || key == Qt::Key_Down;
    if (arrow && (e->modifiers() & Qt::AltModifier)) {
        if (key == Qt::Key_Left || key == Qt::Key_Right) rotateSelection(key == Qt::Key_Left ? -15 : 15);
        return;
    }
    if (arrow && !(e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)))
        return moveHotspot(arrowDirection(key), e->modifiers() & Qt::ShiftModifier);
    if (arrow && (e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)))
        return duplicateSelection(arrowDirection(key));
    if (key == Qt::Key_Escape) {
        hoverAtom_ = hoverBond_ = -1;
        return setSelection({});
    }
    if (key == Qt::Key_Delete || key == Qt::Key_Backspace) {
        if (!selectedAtoms_.isEmpty() || !selectedArrows_.isEmpty() || !selectedTexts_.isEmpty())
            return deleteSelection();
        Document next = doc_;
        if (hoverAtom_ >= 0) {
            Atom& a = next.atoms[hoverAtom_];
            // ChemDraw: removes a label first; a plain carbon is deleted.
            if (a.z != 6 || a.charge || !a.label.isEmpty()) a.z = 6, a.charge = 0, a.label.clear();
            else next.removeAtoms({hoverAtom_}), hoverAtom_ = -1;
        } else if (hoverBond_ >= 0) {
            next.bonds.erase(next.bonds.begin() + hoverBond_);
            hoverBond_ = -1;
        } else {
            return;
        }
        return commit(next, tr("Delete"));
    }
    if (e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)) return QGraphicsView::keyPressEvent(e);

    const QString t = e->text();
    Document next = doc_;

    if (hoverAtom_ >= 0) {
        const int at = hoverAtom_;
        if (key == Qt::Key_Return || key == Qt::Key_Enter || t == "=") return editLabel(at);
        if (t == "+" || t == "-") {
            next.atoms[at].charge += t == "+" ? 1 : -1;
            return commit(next, tr("Charge"));
        }
        int hot = sproutHotkey(next, at, t);
        if (hot == kUnhandled) {
            QString label = labelHotkey(t);
            if (label.isEmpty() || !applyLabel(next, at, label)) return QGraphicsView::keyPressEvent(e);
            hot = at;
        }
        commit(next, tr("Hotkey %1").arg(t));
        hoverAtom_ = hot, hoverBond_ = -1;
        viewport()->update();
        return;
    }
    if (hoverAtom_ < 0 && hoverBond_ < 0) {  // no hotspot: tool keys
        static const QStringList tools{"x", "X", "j", "t", "e", " "};
        if (tools.contains(t)) return emit toolKey(t);
    }
    if (hoverBond_ >= 0) {
        Bond& b = next.bonds[hoverBond_];
        static const QHash<QString, std::pair<int, bool>> fuse{
            {"a", {6, true}}, {"z", {5, true}}, {"v", {3, false}}, {"4", {4, false}},
            {"5", {5, false}}, {"6", {6, false}}, {"7", {7, false}}, {"8", {8, false}}};
        if (t == "1" || t == "2" || t == "3") {
            b.order = t.toInt(), b.stereo = BondStereo::None;
            straightenSp(next, hoverBond_);
        } else if (t == "w" || t == "h" || t == "H") {
            BondStereo s = t == "w" ? BondStereo::Wedge : BondStereo::Hash;
            if (b.stereo == s) std::swap(b.a, b.b);  // again: flip which end is narrow
            b.stereo = s, b.order = 1;
        } else if (fuse.contains(t)) {
            ringOnBond(next, hoverBond_, fuse[t].first, fuse[t].second);
        } else if (t == "9" || t == "0") {
            chairOnBond(next, hoverBond_, t == "9" ? 0 : 1);
        } else if (t == "d" || t == "b" || t == "y" || t == "D" || t == "B") {
            static const QHash<QString, BondStereo> styles{{"d", BondStereo::Dashed}, {"b", BondStereo::Bold},
                                                           {"y", BondStereo::Wavy},   {"D", BondStereo::Dashed},
                                                           {"B", BondStereo::Bold}};
            b.stereo = styles[t];
            b.order = t == "D" || t == "B" ? 2 : 1;
        } else if (t == "l" || t == "c" || t == "r") {
            if (b.order != 2) b.order = 2, b.stereo = BondStereo::None;
            b.position = t == "l" ? BondPosition::Left : t == "c" ? BondPosition::Centre : BondPosition::Right;
        } else {
            return QGraphicsView::keyPressEvent(e);
        }
        const int keep = hoverBond_;
        commit(next, tr("Hotkey %1").arg(t));
        hoverBond_ = keep;
        viewport()->update();
        return;
    }
    QGraphicsView::keyPressEvent(e);
}

void Canvas::wheelEvent(QWheelEvent* e) {
    if (e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)) {
        zoomBy(std::pow(1.0015, e->angleDelta().y()));
        return;
    }
    QGraphicsView::wheelEvent(e);
}
