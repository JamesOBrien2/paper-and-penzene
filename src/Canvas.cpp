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
#include <QScrollBar>
#include <QUndoStack>
#include <QtMath>
#include <algorithm>
#include <functional>

// ACS 1996 document settings, in points.
constexpr double kLineWidth = 0.6;
constexpr double kBondSpacing = 0.18 * kBondLength;  // double-bond gap
constexpr double kWedgeWidth = 4.5;
constexpr double kHashSpacing = 2.2;
constexpr double kFontSize = 10;
constexpr double kLabelRadius = 5.5;  // bonds stop short of labels
constexpr double kMergeRadius = 0.3 * kBondLength;

static double len(QPointF v) { return std::hypot(v.x(), v.y()); }
static QPointF unit(QPointF v) { double l = len(v); return l > 1e-9 ? v / l : QPointF(1, 0); }
static QPointF perp(QPointF v) { return {-v.y(), v.x()}; }
static double cross(QPointF a, QPointF b) { return a.x() * b.y() - a.y() * b.x(); }
static QPointF dirAt(double deg) { return {std::cos(qDegreesToRadians(deg)), std::sin(qDegreesToRadians(deg))}; }

static bool hasLabel(const Document& doc, int i, const std::vector<int>& degree) {
    return doc.atoms[i].z != 6 || degree[i] == 0;
}

// ---------------------------------------------------------------- rendering

static QFont labelFont() {
    QFont f("Arial");
    f.setPixelSize(int(kFontSize));  // 1 px == 1 pt in our scene units
    return f;
}

// Text as outlines, so it scales identically on screen, SVG, PDF and PNG.
static void drawText(QPainter& p, const QString& s, QPointF baselineLeft, const QFont& f) {
    QPainterPath path;
    path.addText(baselineLeft, f, s);
    p.fillPath(path, p.pen().color());
}

static void drawLabel(QPainter& p, const Document& doc, int i, int hydrogens, bool hLeft) {
    const auto& a = doc.atoms[i];
    QFont f = labelFont(), sub = f;
    sub.setPixelSize(int(kFontSize * 0.7));
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

static void drawBond(QPainter& p, const Document& doc, const Bond& b, const std::vector<int>& degree,
                     const std::vector<bool>& labeled) {
    QPointF pa = doc.atoms[b.a].pos, pb = doc.atoms[b.b].pos;
    QPointF d = unit(pb - pa), n = perp(d);
    // Trim at labels.
    QPointF a = labeled[b.a] ? pa + d * kLabelRadius : pa;
    QPointF e = labeled[b.b] ? pb - d * kLabelRadius : pb;

    if (b.stereo == BondStereo::Wedge) {
        QPolygonF tri{a, e + n * kWedgeWidth / 2, e - n * kWedgeWidth / 2};
        p.setBrush(p.pen().color());
        p.drawPolygon(tri);
        p.setBrush(Qt::NoBrush);
        return;
    }
    if (b.stereo == BondStereo::Hash) {
        double L = len(e - a);
        int count = std::max(3, int(L / kHashSpacing));
        for (int k = 0; k <= count; ++k) {
            double t = double(k) / count;
            QPointF c = a + (e - a) * t;
            double w = kWedgeWidth / 2 * t;
            p.drawLine(c + n * w, c - n * w);
        }
        return;
    }

    if (b.order == 1) {
        p.drawLine(a, e);
    } else if (b.order == 3) {
        p.drawLine(a, e);
        p.drawLine(a + n * kBondSpacing, e + n * kBondSpacing);
        p.drawLine(a - n * kBondSpacing, e - n * kBondSpacing);
    } else {
        // Offset the second line toward the side where the neighbours are
        // (inside the ring); centre it for terminal bonds like C=O.
        double side = 0;
        for (int end : {b.a, b.b})
            for (int nb : doc.neighbors(end))
                if (nb != b.a && nb != b.b) side += cross(d, doc.atoms[nb].pos - pa) > 0 ? 1 : -1;
        // Neighbours on opposite sides (trans chain) tie at 0: still offset, or
        // both lines would cross into the adjoining single bonds.
        bool centred = degree[b.a] == 1 || degree[b.b] == 1;
        if (centred) {
            QPointF o = n * kBondSpacing / 2;
            p.drawLine(a + o, e + o);
            p.drawLine(a - o, e - o);
        } else {
            QPointF o = n * (side >= 0 ? kBondSpacing : -kBondSpacing);
            QPointF shrink = d * (0.15 * kBondLength);
            QPointF ia = labeled[b.a] ? a : a + shrink, ie = labeled[b.b] ? e : e - shrink;
            p.drawLine(a, e);
            p.drawLine(ia + o, ie + o);
        }
    }
}

void paintDocument(QPainter& p, const Document& doc, const RenderStyle& style) {
    p.save();
    p.setRenderHint(QPainter::Antialiasing);
    std::vector<int> degree(doc.atoms.size(), 0);
    for (const auto& b : doc.bonds) ++degree[b.a], ++degree[b.b];
    std::vector<bool> labeled(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) labeled[i] = hasLabel(doc, int(i), degree);
    auto info = chem::atomInfo(doc);

    QPen pen(style.ink, kLineWidth, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
    p.setPen(pen);
    for (const auto& b : doc.bonds) drawBond(p, doc, b, degree, labeled);

    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        p.setPen(QPen(info[i].valenceError ? style.error : style.ink, kLineWidth));
        if (labeled[i]) {
            // H goes on the side away from the bonds.
            double dx = 0;
            for (int nb : doc.neighbors(int(i))) dx += doc.atoms[nb].pos.x() - a.pos.x();
            drawLabel(p, doc, int(i), info[i].hydrogens, dx > 0.1);
        } else if (a.charge) {
            QFont sub = labelFont();
            sub.setPixelSize(int(kFontSize * 0.7));
            QString c = QString(a.charge > 0 ? "+" : "−");
            if (std::abs(a.charge) > 1) c.prepend(QString::number(std::abs(a.charge)));
            drawText(p, c, a.pos + QPointF(2, -3), sub);
        }
        if (info[i].valenceError && !labeled[i]) p.drawEllipse(a.pos, 3, 3);
    }
    p.restore();
}

QRectF documentBounds(const Document& doc) {
    if (doc.atoms.empty()) return {};
    // Not QRectF::united: it ignores zero-size rects.
    QPointF lo = doc.atoms[0].pos, hi = lo;
    for (const auto& a : doc.atoms) {
        lo = {std::min(lo.x(), a.pos.x()), std::min(lo.y(), a.pos.y())};
        hi = {std::max(hi.x(), a.pos.x()), std::max(hi.y(), a.pos.y())};
    }
    return QRectF(lo, hi).adjusted(-kFontSize * 1.5, -kFontSize, kFontSize * 1.5, kFontSize);  // room for labels
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
    gen.setTitle("Paper & Penzene");
    QPainter p(&gen);
    p.translate(-r.topLeft());
    paintDocument(p, doc);
    p.end();
    return buf.data();
}

bool exportDocument(const Document& doc, const QString& path) {
    QRectF r = documentBounds(doc);
    if (r.isEmpty()) return false;
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
        pdf.setCreator("Paper & Penzene");
        QPainter p(&pdf);
        p.translate(-r.topLeft());
        paintDocument(p, doc);
        return true;
    }
    return false;
}

// ---------------------------------------------------------------- editing helpers

namespace {

QPointF awayDirection(const Document& doc, int atom);

// Direction for a new bond from `atom` that avoids existing bonds.
QPointF freeDirection(const Document& doc, int atom) {
    auto nbs = doc.neighbors(atom);
    QPointF p = doc.atoms[atom].pos;
    if (nbs.empty()) return dirAt(-30);
    if (nbs.size() == 1) {
        QPointF back = unit(doc.atoms[nbs[0]].pos - p);
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
    return awayDirection(doc, atom);
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

// Direction pointing away from all of the atom's bonds: the bisector of the
// widest gap between them (straight on for a terminal atom).
QPointF awayDirection(const Document& doc, int atom) {
    auto nbs = doc.neighbors(atom);
    if (nbs.empty()) return {0, -1};
    QPointF p = doc.atoms[atom].pos;
    std::vector<double> ang;
    for (int nb : nbs) ang.push_back(std::atan2(doc.atoms[nb].pos.y() - p.y(), doc.atoms[nb].pos.x() - p.x()));
    std::sort(ang.begin(), ang.end());
    double bestGap = -1, bestMid = 0;
    for (size_t i = 0; i < ang.size(); ++i) {
        double next = i + 1 < ang.size() ? ang[i + 1] : ang[0] + 2 * M_PI;
        if (next - ang[i] > bestGap + 1e-6) bestGap = next - ang[i], bestMid = (ang[i] + next) / 2;
    }
    return {std::cos(bestMid), std::sin(bestMid)};
}

double circumradius(int n) { return kBondLength / (2 * std::sin(M_PI / n)); }

void ringAt(Document& doc, QPointF centre, int n, bool aromatic) {
    addRing(doc, polygon(centre, centre + QPointF(0, -circumradius(n)), n), aromatic);
}

// Ring through the atom, pointing away from its bonds so they bisect the ring's outside angle.
std::vector<int> ringOnAtom(Document& doc, int atom, int n, bool aromatic) {
    QPointF p = doc.atoms[atom].pos;
    return addRing(doc, polygon(p + awayDirection(doc, atom) * circumradius(n), p, n), aromatic);
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
    QSet<int> keep;
    for (int i : selectedAtoms_)
        if (i < int(doc_.atoms.size())) keep.insert(i);
    selectedAtoms_ = keep;
    if (hoverAtom_ >= int(doc_.atoms.size())) hoverAtom_ = -1;
    if (hoverBond_ >= int(doc_.bonds.size())) hoverBond_ = -1;
    refresh();
    emit documentChanged();
}

void Canvas::setSelection(QSet<int> atoms) {
    selectedAtoms_ = std::move(atoms);
    viewport()->update();
}

void Canvas::selectAll() {
    QSet<int> all;
    for (int i = 0; i < int(doc_.atoms.size()); ++i) all.insert(i);
    setSelection(all);
}

Document Canvas::selectedSubset() const {
    if (selectedAtoms_.isEmpty()) return doc_;
    std::vector<int> drop;
    for (int i = 0; i < int(doc_.atoms.size()); ++i)
        if (!selectedAtoms_.contains(i)) drop.push_back(i);
    Document out = doc_;
    out.removeAtoms(drop);
    return out;
}

void Canvas::deleteSelection() {
    if (selectedAtoms_.isEmpty()) return;
    Document next = doc_;
    next.removeAtoms({selectedAtoms_.begin(), selectedAtoms_.end()});
    selectedAtoms_.clear();
    commit(next, tr("Delete"));
}

void Canvas::insert(Document frag, const QString& text) {
    if (frag.atoms.empty()) return;
    QPointF c;
    for (const auto& a : frag.atoms) c += a.pos;
    QPointF shift = viewCenter() - c / double(frag.atoms.size());
    Document next = doc_;
    const int base = int(next.atoms.size());
    QSet<int> added;
    for (auto a : frag.atoms) {
        a.pos += shift;
        added.insert(int(next.atoms.size()));
        next.atoms.push_back(a);
    }
    for (auto b : frag.bonds) {
        b.a += base, b.b += base;
        next.bonds.push_back(b);
    }
    commit(next, text);
    setSelection(added);
}

QPointF Canvas::viewCenter() const { return mapToScene(viewport()->rect().center()); }

void Canvas::zoomBy(double factor) {
    double s = transform().m11() * factor;
    if (s > 0.2 && s < 40) scale(factor, factor);
}

void Canvas::fitToDocument() {
    if (doc_.atoms.empty()) return;
    fitInView(documentBounds(doc_).adjusted(-20, -20, 20, 20), Qt::KeepAspectRatio);
}

// Cache the drawing as a QPicture; hover/selection repaints just replay it.
void Canvas::refresh() {
    picture_ = QPicture();
    QPainter p(&picture_);
    paintDocument(p, doc_);
    p.end();
    viewport()->update();
}

void Canvas::drawBackground(QPainter* p, const QRectF& rect) {
    p->fillRect(rect, Qt::white);
    picture_.play(p);
}

void Canvas::drawForeground(QPainter* p, const QRectF&) {
    const QColor sel(40, 120, 255, 90), hover(40, 120, 255, 60);
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

    p->setBrush(hover);
    if (hoverAtom_ >= 0) {
        p->drawEllipse(doc_.atoms[hoverAtom_].pos, 5, 5);
        p->setBrush(QColor(40, 170, 60));
        p->drawEllipse(doc_.atoms[hoverAtom_].pos, 1.2, 1.2);
    } else if (hoverBond_ >= 0) {
        const auto& b = doc_.bonds[hoverBond_];
        p->setPen(QPen(hover, 5, Qt::SolidLine, Qt::RoundCap));
        p->drawLine(doc_.atoms[b.a].pos, doc_.atoms[b.b].pos);
    }

    p->setBrush(Qt::NoBrush);
    if (drag_ == Drag::Rubber) {
        p->setPen(QPen(QColor(40, 120, 255), 0, Qt::DashLine));
        p->drawRect(QRectF(pressPos_, curPos_).normalized());
    } else if (drag_ == Drag::Bond || drag_ == Drag::Chain) {
        p->setPen(QPen(QColor(40, 120, 255), 0.8));
        for (size_t k = 1; k < preview_.size(); ++k) p->drawLine(preview_[k - 1], preview_[k]);
    }
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
    case Tool::Select:
        if (pressAtom_ >= 0 || bond >= 0) {
            QSet<int> hit = pressAtom_ >= 0 ? QSet<int>{pressAtom_}
                                            : QSet<int>{doc_.bonds[bond].a, doc_.bonds[bond].b};
            bool already = std::all_of(hit.begin(), hit.end(), [&](int i) { return selectedAtoms_.contains(i); });
            if (e->modifiers() & Qt::ShiftModifier) selectedAtoms_ |= hit;
            else if (!already) selectedAtoms_ = hit;
            drag_ = (e->modifiers() & Qt::AltModifier) ? Drag::Rotate : Drag::Move;
        } else {
            if (!(e->modifiers() & Qt::ShiftModifier)) selectedAtoms_.clear();
            drag_ = Drag::Rubber;
        }
        break;
    case Tool::Bond: case Tool::Wedge: case Tool::Hash:
        drag_ = Drag::Bond;
        break;
    case Tool::Chain:
        drag_ = Drag::Chain;
        break;
    default:
        drag_ = Drag::None;  // click tools act on release
    }
    viewport()->update();
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
        QPointF c;
        for (int i : selectedAtoms_) c += beforeDrag_.atoms[i].pos;
        c /= std::max<qsizetype>(1, selectedAtoms_.size());
        double ang = std::atan2(curPos_.y() - c.y(), curPos_.x() - c.x()) -
                     std::atan2(pressPos_.y() - c.y(), pressPos_.x() - c.x());
        for (int i : selectedAtoms_) {
            QPointF& p = next.atoms[i].pos;
            if (drag_ == Drag::Move) {
                p += curPos_ - pressPos_;
            } else {
                QPointF r = p - c;
                p = c + QPointF(r.x() * std::cos(ang) - r.y() * std::sin(ang),
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
    } else if ((drag == Drag::Bond || drag == Drag::Chain) && click) {
        if (bond >= 0) {  // click on a bond: change it in place
            Bond& b = next.bonds[bond];
            if (stereo != BondStereo::None) {
                if (b.stereo == stereo) std::swap(b.a, b.b);  // flip direction
                b.stereo = stereo, b.order = 1;
            } else {
                b.stereo = BondStereo::None;
                b.order = (b.order != order && order > 1) ? order : b.order % 3 + 1;
            }
            what = tr("Change bond");
        } else {
            int from = pressAtom_ >= 0 ? pressAtom_ : next.addAtom(pressPos_);
            QPointF to = next.atoms[from].pos + freeDirection(next, from) * kBondLength;
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
            what = tr("Erase");
            break;
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
}

void Canvas::mouseDoubleClickEvent(QMouseEvent* e) {
    if (tool_ != Tool::Select) return QGraphicsView::mouseDoubleClickEvent(e);
    int start = atomAt(mapToScene(e->pos()));
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

// Groups for label hotkeys and the Enter dialog, as SMILES whose first atom
// replaces the hotspot atom. ponytail: drawn out in full until abbreviations (#29).
const QHash<QString, QString>& groups() {
    static const QHash<QString, QString> g{
        {"Me", "C"},        {"Et", "CC"},          {"iPr", "C(C)C"},   {"tBu", "C(C)(C)C"},
        {"Ph", "c1ccccc1"}, {"OH", "O"},           {"OMe", "OC"},      {"NH2", "N"},
        {"NO2", "[N+](=O)[O-]"}, {"SH", "S"},      {"CF3", "C(F)(F)F"}, {"CN", "C#N"},
        {"CO2Me", "C(=O)OC"}, {"CO2H", "C(=O)O"},  {"CHO", "C=O"},     {"Ac", "C(C)=O"},
        {"OAc", "OC(C)=O"}, {"N3", "N=[N+]=[N-]"}, {"Boc", "C(=O)OC(C)(C)C"},
        {"Cbz", "C(=O)OCc1ccccc1"}, {"Fmoc", "C(=O)OCC1c2ccccc2-c2ccccc21"},
        {"MgBr", "[Mg]Br"}, {"SO2Me", "S(=O)(=O)C"}, {"Ts", "S(=O)(=O)c1ccc(C)cc1"},
    };
    return g;
}

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

// Replaces atom `at` with the group's first atom and lays the rest out pointing
// away from the atom's existing bonds. Returns false if `label` is unknown.
bool Canvas::applyLabel(Document& doc, int at, const QString& label) {
    if (int z = chem::atomicNumber(label.toStdString()); z > 0) {
        doc.atoms[at].z = z;
        return true;
    }
    QString smiles = groups().value(label, label);  // anything else: try it as SMILES
    auto frag = chem::fromSmiles(smiles.toStdString());
    if (!frag || frag->atoms.empty()) return false;
    doc.atoms[at].z = frag->atoms[0].z;
    doc.atoms[at].charge = frag->atoms[0].charge;
    if (frag->atoms.size() == 1) return true;

    // Rotate so the fragment's bulk points away from `at`'s bonds.
    QPointF origin = frag->atoms[0].pos, bulk;
    for (size_t i = 1; i < frag->atoms.size(); ++i) bulk += frag->atoms[i].pos;
    bulk = bulk / double(frag->atoms.size() - 1) - origin;
    QPointF want = doc.neighbors(at).empty() ? QPointF(1, 0) : awayDirection(doc, at);
    double ang = std::atan2(want.y(), want.x()) - std::atan2(bulk.y(), bulk.x());
    const int base = int(doc.atoms.size()) - 1;  // frag atom i -> base + i
    for (size_t i = 1; i < frag->atoms.size(); ++i) {
        QPointF r = frag->atoms[i].pos - origin;
        Atom a = frag->atoms[i];
        a.pos = doc.atoms[at].pos + QPointF(r.x() * std::cos(ang) - r.y() * std::sin(ang),
                                            r.x() * std::sin(ang) + r.y() * std::cos(ang));
        doc.atoms.push_back(a);
    }
    for (auto b : frag->bonds) {
        b.a = b.a == 0 ? at : base + b.a;
        b.b = b.b == 0 ? at : base + b.b;
        doc.bonds.push_back(b);
    }
    return true;
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
int linker(Document& doc, int at) { return sprout(doc, at, awayDirection(doc, at)); }

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
        return sprout(doc, at, awayDirection(doc, at), 1, 6, BondStereo::None, 1.5 * kBondLength);
    }
    if (key == "2") {  // carbonyl / acetyl
        if (s == Site::Secondary) {
            sprout(doc, at, awayDirection(doc, at), 2, 8);
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
            sprout(doc, at, awayDirection(doc, at), 1, 6, st);
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
        QPointF dir = site(doc, c) == Site::Primary ? freeDirection(doc, c) : awayDirection(doc, c);
        return sprout(doc, c, dir, 2);
    }
    if (key == "9") {  // dimethyl / gem-dimethyl / isopropyl
        if (s == Site::Secondary) {
            QPointF away = awayDirection(doc, at);
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
    if (selectedAtoms_.isEmpty()) return;
    Document next = doc_;
    QPointF c;
    for (int i : selectedAtoms_) c += doc_.atoms[i].pos;
    c /= double(selectedAtoms_.size());
    for (int i : selectedAtoms_) next.atoms[i].pos = c + rotated(doc_.atoms[i].pos - c, degrees);
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
                                          QString::fromStdString(chem::symbol(doc_.atoms[at].z)), &ok)
                        .trimmed();
    Document next = doc_;
    if (ok && !label.isEmpty() && applyLabel(next, at, label)) commit(next, tr("Edit label"));
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
    if (key == Qt::Key_Escape) {
        hoverAtom_ = hoverBond_ = -1;
        return setSelection({});
    }
    if (key == Qt::Key_Delete || key == Qt::Key_Backspace) {
        if (!selectedAtoms_.isEmpty()) return deleteSelection();
        Document next = doc_;
        if (hoverAtom_ >= 0) {
            Atom& a = next.atoms[hoverAtom_];
            // ChemDraw: removes a label first; a plain carbon is deleted.
            if (a.z != 6 || a.charge) a.z = 6, a.charge = 0;
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
    if (hoverBond_ >= 0) {
        Bond& b = next.bonds[hoverBond_];
        static const QHash<QString, std::pair<int, bool>> fuse{
            {"a", {6, true}}, {"z", {5, true}}, {"v", {3, false}}, {"4", {4, false}},
            {"5", {5, false}}, {"6", {6, false}}, {"7", {7, false}}, {"8", {8, false}}};
        if (t == "1" || t == "2" || t == "3") {
            b.order = t.toInt(), b.stereo = BondStereo::None;
        } else if (t == "w" || t == "h" || t == "H") {
            BondStereo s = t == "w" ? BondStereo::Wedge : BondStereo::Hash;
            if (b.stereo == s) std::swap(b.a, b.b);  // again: flip which end is narrow
            b.stereo = s, b.order = 1;
        } else if (fuse.contains(t)) {
            ringOnBond(next, hoverBond_, fuse[t].first, fuse[t].second);
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
