#include "Canvas.h"
#include "Chem.h"
#include "Edit.h"
#include "Geometry.h"

#include <QImage>
#include <QInputDialog>
#include <QApplication>
#include <QClipboard>
#include <QContextMenuEvent>
#include <QKeyEvent>
#include <QMenu>
#include <QLineEdit>
#include <QComboBox>
#include <QSpinBox>
#include <QFormLayout>
#include <QDialogButtonBox>
#include <QDialog>
#include <array>
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

using namespace edit;

namespace {

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

// ---------------------------------------------------------------- canvas

Canvas::Canvas(QUndoStack* undo, QWidget* parent) : QGraphicsView(parent), undo_(undo) {
    setScene(new QGraphicsScene(-5000, -5000, 10000, 10000, this));
    setMouseTracking(true);
    setAccessibleName(tr("Drawing"));
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
    announceHotspot();
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

// Positions only: labels, text and line widths keep their size. Arrows keep
// their curve's shape (its bend scales with the drawing, and flips with a mirror).
static void applyTransform(Document& d, const QSet<int>& atoms, const QSet<int>& arrows, const QSet<int>& texts,
                           const QTransform& t) {
    for (int i : atoms) d.atoms[i].pos = t.map(d.atoms[i].pos);
    for (int i : arrows) {
        Arrow& a = d.arrows[i];
        const double oldLen = len(a.to - a.from);
        a.from = t.map(a.from), a.to = t.map(a.to);
        if (oldLen > 1e-9) a.bend *= len(a.to - a.from) / oldLen * (t.determinant() < 0 ? -1 : 1);
    }
    for (int i : texts) d.texts[i].pos = t.map(d.texts[i].pos);
}

// How far a structure is turned from the way it's drawn (bonds on the 30° grid), in degrees,
// -15 to 15: the circular mean of its bond (and arrow) angles, modulo 30°. 0 if it has none.
static double orientation(const Document& d, const QSet<int>& atoms, const QSet<int>& arrows) {
    double s = 0, c = 0;
    auto add = [&](QPointF v) {
        const double a = 12 * std::atan2(v.y(), v.x());  // 30° is a whole turn here
        s += std::sin(a), c += std::cos(a);
    };
    for (const Bond& b : d.bonds)
        if (atoms.contains(b.a) && atoms.contains(b.b)) add(d.atoms[b.b].pos - d.atoms[b.a].pos);
    for (int i : arrows) add(d.arrows[i].to - d.arrows[i].from);
    return std::hypot(s, c) < 1e-9 ? 0 : qRadiansToDegrees(std::atan2(s, c)) / 12;
}

static std::array<QPointF, 8> handlePoints(const QRectF& r) {
    const QPointF c = r.center();
    return {r.topLeft(), {c.x(), r.top()}, r.topRight(), {r.right(), c.y()},
            r.bottomRight(), {c.x(), r.bottom()}, r.bottomLeft(), {r.left(), c.y()}};
}


void Canvas::selectAll() {
    setSelection(range(0, int(doc_.atoms.size())), range(0, int(doc_.arrows.size())), range(0, int(doc_.texts.size())));
}

// Drops every item not in the given sets.
static Document keepOnly(const Document& doc, const QSet<int>& atoms, const QSet<int>& arrows, const QSet<int>& texts) {
    Document out = doc;
    out.page.clear();  // a selection exports as just the drawing
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
    if (const QRectF page = pageRect(doc_); !page.isEmpty()) {  // the page's edge, and its margins dashed
        QColor edge = theme_.ink;
        edge.setAlphaF(0.35);
        p->setPen(QPen(edge, 0));
        p->setBrush(Qt::NoBrush);
        p->drawRect(page);
        for (const auto& s : pageSizes())
            if (s.name == doc_.page && s.margin > 0) {
                const double m = s.margin / exportScale(doc_);
                p->setPen(QPen(edge, 0, Qt::DashLine));
                p->drawRect(page.adjusted(m, m, -m, -m));
            }
    }
    picture_.play(p);
}

// Screen readers hear the hotspot as the canvas's description, updated after every key and
// mouse event (the only things that move it) and document change.
bool Canvas::event(QEvent* e) {
    const bool done = QGraphicsView::event(e);
    announceHotspot();
    return done;
}

bool Canvas::viewportEvent(QEvent* e) {
    const bool done = QGraphicsView::viewportEvent(e);
    announceHotspot();
    return done;
}

void Canvas::announceHotspot() {
    auto name = [this](int i) {
        const Atom& a = doc_.atoms[i];
        return (a.label.isEmpty() ? QString::fromStdString(chem::symbol(a.z)) : a.label) + QString::number(i + 1);
    };
    QString text;
    if (hoverAtom_ >= 0 && hoverAtom_ < int(doc_.atoms.size()))
    {
        const int n = int(doc_.neighbors(hoverAtom_).size());
        text = n == 1 ? tr("Hotspot: atom %1, 1 bond").arg(name(hoverAtom_))
                      : tr("Hotspot: atom %1, %2 bonds").arg(name(hoverAtom_)).arg(n);
    }
    else if (hoverBond_ >= 0 && hoverBond_ < int(doc_.bonds.size())) {
        const Bond& b = doc_.bonds[hoverBond_];
        const QString order = b.order == 3 ? tr("triple") : b.order == 2 ? tr("double") : tr("single");
        text = tr("Hotspot: %1 bond, %2 to %3").arg(order, name(b.a), name(b.b));
    }
    if (text != accessibleDescription()) setAccessibleDescription(text);  // also after an edit in place
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
    for (int i : selectedTexts_) p->drawRect(textPath(doc_.texts[i], documentStyle(doc_)).boundingRect().adjusted(-1.5, -1.5, 1.5, 1.5));

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

    if (tool_ == Tool::Select && (drag_ == Drag::None || drag_ == Drag::Scale)) {
        if (const QRectF box = selectionBox(); !box.isNull()) {
            p->setBrush(Qt::NoBrush);
            p->setPen(QPen(line, 0, Qt::DotLine));
            p->drawRect(box);
            p->setPen(QPen(line, 0));
            p->setBrush(theme_.paper);
            const double s = 2.5 / transform().m11();
            for (QPointF h : handlePoints(box)) p->drawRect(QRectF(h - QPointF(s, s), h + QPointF(s, s)));
            const QPointF knob = *rotateHandle();
            p->drawLine(QPointF(box.center().x(), box.top()), knob + QPointF(0, 1.4 * s));
            p->drawEllipse(knob, 1.4 * s, 1.4 * s);
        }
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
    } else if (drag_ == Drag::Ring && len(curPos_ - pressPos_) >= 3 / transform().m11()) {
        Document ring = doc_;
        addDraggedRing(ring);
        p->setPen(QPen(line, 0.8));
        for (const auto& b : ring.bonds)
            if (doc_.bondBetween(b.a, b.b) < 0 || b.a >= int(doc_.atoms.size()) || b.b >= int(doc_.atoms.size()))
                p->drawLine(ring.atoms[b.a].pos, ring.atoms[b.b].pos);
        QFont f = font();
        f.setPixelSize(8);
        p->setFont(f);
        p->drawText(curPos_ + QPointF(6, -6), QString::number(draggedRingSize()));  // the size, by the cursor
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
        if (textPath(doc_.texts[i], documentStyle(doc_)).boundingRect().adjusted(-2, -2, 2, 2).contains(p)) return i;
    return -1;
}

// A dragged ring grows by one atom per half bond length of drag, from 3 to 12.
int Canvas::draggedRingSize() const {
    return std::clamp(3 + int(len(curPos_ - pressPos_) / (0.5 * kBondLength)), 3, 12);
}

// Fused onto the bond or attached to the atom where the drag began; on empty
// space the ring hangs off the press point in the drag direction. Shift: aromatic.
void Canvas::addDraggedRing(Document& doc) const {
    const int n = draggedRingSize();
    const bool aromatic = shift_;
    if (pressAtom_ >= 0) ringOnAtom(doc, pressAtom_, n, aromatic);
    else if (int bond = bondAt(pressPos_); bond >= 0) ringOnBond(doc, bond, n, aromatic);
    else addRing(doc, polygon(pressPos_ + unit(curPos_ - pressPos_) * circumradius(n), pressPos_, n), aromatic);
}

// The arrow being dragged out: straight ones snap to 15°, curved ones bow left.
Arrow Canvas::draggedArrow() const {
    Arrow a{pressPos_, curPos_, arrowKind_};
    a.dashed = arrowDashed_;
    if (isShape(arrowKind_) && arrowKind_ != ArrowKind::Line) {  // boxes and ellipses: any corner, no snapping
        if (shift_) {  // Shift: a square or circle
            const QPointF v = curPos_ - pressPos_;
            const double side = std::max(std::abs(v.x()), std::abs(v.y()));
            a.to = pressPos_ + QPointF(v.x() < 0 ? -side : side, v.y() < 0 ? -side : side);
        }
    } else if (arrowCurved_) {
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
        const QPointF dir = shift_ ? unit(curPos_ - start) : snapped(start, curPos_);  // Shift: any angle
        return {start, start + dir * kBondLength};
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
        if (const auto knob = rotateHandle(); knob && len(*knob - pressPos_) < 5 / transform().m11()) {
            drag_ = Drag::Rotate;
            break;
        }
        if (int h = handleAt(pressPos_); h >= 0) {  // a scale handle wins over what's under it
            drag_ = Drag::Scale, scaleHandle_ = h, scaleBox_ = selectionBox();
            break;
        }
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
            drag_ = !(e->modifiers() & Qt::AltModifier) ? Drag::Move : shift ? Drag::Rotate3D : Drag::Rotate;
            if (drag_ == Drag::Rotate3D && !(pose_ = chem::pose3D(doc_, moleculesOfSelection()))) drag_ = Drag::Rotate;
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
    case Tool::Ring:  // a click adds the chosen ring; a drag sizes one
        drag_ = Drag::Ring;
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
    shift_ = e->modifiers() & Qt::ShiftModifier;
    if (drag_ == Drag::Move && shift_) {  // Shift: along whichever axis the drag mostly follows
        QPointF d = curPos_ - pressPos_;
        curPos_ = pressPos_ + (std::abs(d.x()) >= std::abs(d.y()) ? QPointF(d.x(), 0) : QPointF(0, d.y()));
    }
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
        if (drag_ == Drag::Rotate && e->modifiers() & Qt::ControlModifier) {
            // Ctrl: the structure lands square to the page or at a multiple of 45° from it, whatever
            // angle it started at (one drawn 1° off turns to 45°, not 46°).
            const double from = orientation(beforeDrag_, selectedAtoms_, selectedArrows_);
            ang = qDegreesToRadians(45 * std::round((from + qRadiansToDegrees(ang)) / 45) - from);
        } else if (drag_ == Drag::Rotate && shift_) {  // Shift: 15° steps from where it started
            ang = qDegreesToRadians(15 * std::round(qRadiansToDegrees(ang) / 15));
        }
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
    if (drag_ == Drag::Rotate3D) {  // four bond lengths of drag turn it half over
        const double perUnit = 180.0 / (4 * kBondLength);
        doc_ = chem::project3D(beforeDrag_, *pose_, (curPos_.y() - pressPos_.y()) * perUnit, (curPos_.x() - pressPos_.x()) * perUnit);
        refresh();
        return;
    }
    if (drag_ == Drag::Scale) {
        const auto h = handlePoints(scaleBox_);
        const QPointF anchor = h[(scaleHandle_ + 4) % 8], from = h[scaleHandle_] - anchor, to = curPos_ - anchor;
        auto factor = [](double want, double had) { return std::abs(had) < 1e-9 ? 1.0 : std::max(0.05, want / had); };
        double sx = 1, sy = 1;
        if (scaleHandle_ % 2 == 0) {  // corner: uniform, along the diagonal
            sx = sy = std::max(0.05, QPointF::dotProduct(to, from) / QPointF::dotProduct(from, from));
        } else if (scaleHandle_ == 1 || scaleHandle_ == 5) {
            sy = factor(to.y(), from.y());
        } else {
            sx = factor(to.x(), from.x());
        }
        Document next = beforeDrag_;
        applyTransform(next, selectedAtoms_, selectedArrows_, selectedTexts_,
                       QTransform::fromTranslate(-anchor.x(), -anchor.y()) * QTransform::fromScale(sx, sy) *
                           QTransform::fromTranslate(anchor.x(), anchor.y()));
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
    shift_ = e->modifiers() & Qt::ShiftModifier;
    const bool click = len(curPos_ - pressPos_) < 3 / transform().m11();
    const int bond = pressAtom_ < 0 ? bondAt(pressPos_) : -1;
    const Drag drag = std::exchange(drag_, Drag::None);
    preview_.clear();
    Document next = beforeDrag_;
    QString what;

    const BondStereo stereo = tool_ == Tool::Wedge ? BondStereo::Wedge
                              : tool_ == Tool::Hash ? BondStereo::Hash
                                                    : bondStyle_;
    const bool wedge = stereo == BondStereo::Wedge || stereo == BondStereo::Hash;
    const int order = wedge ? 1 : bondOrder_;

    if (drag == Drag::Rotate3D) {
        pose_.reset();
        if (!click) {
            Document turned = doc_;
            doc_ = beforeDrag_;
            commit(turned, tr("Rotate in 3D"));
        }
        return;
    }
    if (drag == Drag::Scale) {
        if (!click) {
            Document scaled = doc_;
            doc_ = beforeDrag_;
            commit(scaled, scaleHandle_ % 2 == 0 ? tr("Scale") : tr("Stretch"));
        }
        return;
    }
    if (drag == Drag::Move || drag == Drag::Rotate) {
        if (!click) {
            Document moved = doc_;
            doc_ = beforeDrag_;
            // A moved atom dropped onto another one fuses with it, as in ChemDraw.
            std::vector<std::pair<int, int>> fuse;
            for (int i : selectedAtoms_)
                if (int j = atomNear(moved, moved.atoms[i].pos, kMergeRadius, i); j >= 0 && !selectedAtoms_.contains(j))
                    fuse.push_back({j, i});
            if (!fuse.empty()) selectedAtoms_.clear(), hoverAtom_ = hoverBond_ = -1;
            mergeAtoms(moved, fuse);
            commit(moved, !fuse.empty() ? tr("Merge") : drag == Drag::Move ? tr("Move") : tr("Rotate"));
        }
        return;
    } else if (drag == Drag::Rubber) {
        QRectF r = QRectF(pressPos_, curPos_).normalized();
        for (int i = 0; i < int(doc_.atoms.size()); ++i)
            if (r.contains(doc_.atoms[i].pos)) selectedAtoms_.insert(i);
        for (int i = 0; i < int(doc_.arrows.size()); ++i)
            if (r.contains(doc_.arrows[i].from) && r.contains(doc_.arrows[i].to)) selectedArrows_.insert(i);
        for (int i = 0; i < int(doc_.texts.size()); ++i)
            if (r.intersects(textPath(doc_.texts[i], documentStyle(doc_)).boundingRect())) selectedTexts_.insert(i);
    } else if (drag == Drag::Arrow) {
        if (int hit = arrowAt(pressPos_); click && hit >= 0) {  // click an arrow: restyle, or flip a curve
            Arrow& a = next.arrows[hit];
            if (arrowCurved_ && a.bend && a.kind == arrowKind_) a.bend = -a.bend;
            else a.kind = arrowKind_, a.bend = arrowCurved_ ? 0.3 * len(a.to - a.from) : 0, a.dashed = arrowDashed_;
        } else {
            if (click)  // default size: an arrow's length, or a box 3 × 2 bonds
                curPos_ = pressPos_ + QPointF(3 * kBondLength, isShape(arrowKind_) && arrowKind_ != ArrowKind::Line ? 2 * kBondLength : 0);
            next.arrows.push_back(draggedArrow());
        }
        what = tr("Arrow");
    } else if ((drag == Drag::Bond || drag == Drag::Chain) && click) {
        if (bond >= 0) {  // click on a bond: change it in place
            Bond& b = next.bonds[bond];
            if (wedge) {
                if (b.stereo == stereo) std::swap(b.a, b.b);  // flip direction
                b.stereo = stereo, b.order = 1;
            } else if (stereo != BondStereo::None) {  // interaction / partial: restyle
                b.stereo = stereo, b.order = order, b.position = BondPosition::Auto;
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
    } else if (drag == Drag::Ring && !click) {
        addDraggedRing(next);
        what = tr("Add ring");
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
            else if (bond >= 0) next.removeBond(bond);
            else if (int a = arrowAt(pressPos_); a >= 0) next.arrows.erase(next.arrows.begin() + a);
            else if (int t = textAt(pressPos_); t >= 0) next.texts.erase(next.texts.begin() + t);
            what = tr("Erase");
            break;
        case Tool::Text:
            if (pressAtom_ >= 0) return editLabel(pressAtom_);
            return editText(textAt(pressPos_), pressPos_);
        case Tool::Colour: {
            // Paint what was clicked; clicking it again in the same colour clears it.
            auto paint = [&](QColor& c) { c = c == colour_ ? QColor() : colour_; };
            if (pressAtom_ >= 0) paint(next.atoms[pressAtom_].color);
            else if (bond >= 0) paint(next.bonds[bond].color);
            else if (int a = arrowAt(pressPos_); a >= 0) paint(next.arrows[a].color);
            else if (int t = textAt(pressPos_); t >= 0) paint(next.texts[t].color);
            what = tr("Colour");
            break;
        }
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

namespace {

// The selection (everything if nothing is selected) as separate objects: each
// molecule, arrow and text is one, so they align and distribute as wholes.
struct Piece {
    QSet<int> atoms, arrows, texts;
};

std::vector<Piece> pieces(const Document& doc, QSet<int> atoms, QSet<int> arrows, QSet<int> texts) {
    if (atoms.isEmpty() && arrows.isEmpty() && texts.isEmpty())
        for (int i = 0; i < int(doc.atoms.size()); ++i) atoms.insert(i);
    if (atoms.size() == int(doc.atoms.size()) && arrows.isEmpty() && texts.isEmpty()) {
        for (int i = 0; i < int(doc.arrows.size()); ++i) arrows.insert(i);
        for (int i = 0; i < int(doc.texts.size()); ++i) texts.insert(i);
    }
    std::vector<Piece> out;
    QSet<int> seen;
    for (int s : atoms) {
        if (seen.contains(s)) continue;
        Piece p;
        std::vector<int> stack{s};
        seen.insert(s);
        while (!stack.empty()) {
            int i = stack.back();
            stack.pop_back();
            p.atoms.insert(i);
            for (int nb : doc.neighbors(i))
                if (atoms.contains(nb) && !seen.contains(nb)) seen.insert(nb), stack.push_back(nb);
        }
        out.push_back(p);
    }
    for (int a : arrows) out.push_back({{}, {a}, {}});
    for (int t : texts) out.push_back({{}, {}, {t}});
    return out;
}

QRectF bounds(const Document& doc, const Piece& p) {
    Document part;
    for (int i : p.atoms) part.atoms.push_back(doc.atoms[i]);
    for (int i : p.arrows) part.arrows.push_back(doc.arrows[i]);
    for (int i : p.texts) part.texts.push_back(doc.texts[i]);
    part.style = doc.style;
    return documentBounds(part);
}

void shiftPiece(Document& doc, const Piece& p, QPointF by) {
    for (int i : p.atoms) doc.atoms[i].pos += by;
    for (int i : p.arrows) doc.arrows[i].from += by, doc.arrows[i].to += by;
    for (int i : p.texts) doc.texts[i].pos += by;
}

}  // namespace

void Canvas::flipSelection(bool horizontal) {
    Piece all;
    for (const Piece& p : pieces(doc_, selectedAtoms_, selectedArrows_, selectedTexts_))
        all.atoms |= p.atoms, all.arrows |= p.arrows, all.texts |= p.texts;
    if (all.atoms.isEmpty() && all.arrows.isEmpty() && all.texts.isEmpty()) return;
    const QPointF c = bounds(doc_, all).center();
    auto mirror = [&](QPointF q) {
        return horizontal ? QPointF(2 * c.x() - q.x(), q.y()) : QPointF(q.x(), 2 * c.y() - q.y());
    };
    Document next = doc_;
    for (int i : all.atoms) next.atoms[i].pos = mirror(doc_.atoms[i].pos);
    for (auto& b : next.bonds)  // a mirror swaps which side a double bond's second line is on
        if (all.atoms.contains(b.a) && all.atoms.contains(b.b))
            b.position = b.position == BondPosition::Left    ? BondPosition::Right
                         : b.position == BondPosition::Right ? BondPosition::Left
                                                             : b.position;
    for (int i : all.arrows) {
        Arrow& a = next.arrows[i];
        a.from = mirror(a.from), a.to = mirror(a.to), a.bend = -a.bend;
    }
    for (int i : all.texts) {  // text moves but reads the right way round
        const QRectF box = textPath(doc_.texts[i], documentStyle(doc_)).boundingRect();
        next.texts[i].pos += mirror(box.center()) - box.center();
    }
    commit(next, horizontal ? tr("Flip Horizontal") : tr("Flip Vertical"));
}

void Canvas::alignSelection(Align edge) {
    const auto ps = pieces(doc_, selectedAtoms_, selectedArrows_, selectedTexts_);
    if (ps.size() < 2) return;
    QRectF all;
    for (const Piece& p : ps) all |= bounds(doc_, p);
    Document next = doc_;
    for (const Piece& p : ps) {
        const QRectF r = bounds(doc_, p);
        QPointF by;
        switch (edge) {
        case Align::Left: by.setX(all.left() - r.left()); break;
        case Align::HCentre: by.setX(all.center().x() - r.center().x()); break;
        case Align::Right: by.setX(all.right() - r.right()); break;
        case Align::Top: by.setY(all.top() - r.top()); break;
        case Align::VCentre: by.setY(all.center().y() - r.center().y()); break;
        case Align::Bottom: by.setY(all.bottom() - r.bottom()); break;
        }
        shiftPiece(next, p, by);
    }
    commit(next, tr("Align"));
}

// Equal gaps between neighbouring objects; the outermost two stay put.
void Canvas::distributeSelection(bool horizontal) {
    auto ps = pieces(doc_, selectedAtoms_, selectedArrows_, selectedTexts_);
    if (ps.size() < 3) return;
    auto lo = [&](const Piece& p) { QRectF r = bounds(doc_, p); return horizontal ? r.left() : r.top(); };
    auto size = [&](const Piece& p) { QRectF r = bounds(doc_, p); return horizontal ? r.width() : r.height(); };
    std::sort(ps.begin(), ps.end(), [&](const Piece& a, const Piece& b) { return lo(a) < lo(b); });
    double used = 0;
    for (const Piece& p : ps) used += size(p);
    const double span = lo(ps.back()) + size(ps.back()) - lo(ps.front());
    const double gap = (span - used) / double(ps.size() - 1);
    Document next = doc_;
    double at = lo(ps.front());
    for (const Piece& p : ps) {
        const double d = at - lo(p);
        shiftPiece(next, p, horizontal ? QPointF(d, 0) : QPointF(0, d));
        at += size(p) + gap;
    }
    commit(next, tr("Distribute"));
}

// A reaction scheme on one baseline: molecules, "+" and straight arrows left to
// right, evenly spaced; whatever sits over or under an arrow (agents, conditions)
// is centred on it, and the arrow grows to fit.
void Canvas::arrangeScheme() {
    const auto ps = pieces(doc_, selectedAtoms_, selectedArrows_, selectedTexts_);
    auto arrowOf = [&](const Piece& p) -> const Arrow* {
        if (!p.atoms.isEmpty() || !p.texts.isEmpty() || p.arrows.size() != 1) return nullptr;
        const Arrow& a = doc_.arrows[*p.arrows.begin()];
        return isShape(a.kind) || a.bend ? nullptr : &a;
    };
    std::vector<int> owner(ps.size(), -1);  // the arrow piece an agent belongs to
    for (size_t i = 0; i < ps.size(); ++i) {
        if (arrowOf(ps[i])) continue;
        const QPointF c = bounds(doc_, ps[i]).center();
        double best = 3 * kBondLength;
        for (size_t j = 0; j < ps.size(); ++j)
            if (const Arrow* a = arrowOf(ps[j])) {
                const double y = (a->from.y() + a->to.y()) / 2, off = std::abs(c.y() - y);
                const bool over = c.x() > std::min(a->from.x(), a->to.x()) && c.x() < std::max(a->from.x(), a->to.x());
                if (over && off > 0.2 * kBondLength && off < best) best = off, owner[i] = int(j);
            }
    }
    std::vector<int> row;
    for (size_t i = 0; i < ps.size(); ++i)
        if (owner[i] < 0) row.push_back(int(i));
    if (row.size() < 2) return;
    std::sort(row.begin(), row.end(), [&](int a, int b) { return bounds(doc_, ps[a]).center().x() < bounds(doc_, ps[b]).center().x(); });
    double baseline = 0;
    int molecules = 0;
    for (int i : row)
        if (!ps[i].atoms.isEmpty()) baseline += bounds(doc_, ps[i]).center().y(), ++molecules;
    baseline = molecules ? baseline / molecules : bounds(doc_, ps[row[0]]).center().y();
    const double gap = kBondLength;
    Document next = doc_;
    double x = bounds(doc_, ps[row[0]]).left();
    for (int i : row) {
        if (const Arrow* a = arrowOf(ps[i])) {
            std::vector<int> above, below;
            for (size_t j = 0; j < ps.size(); ++j)
                if (owner[j] == i) (bounds(doc_, ps[j]).center().y() < (a->from.y() + a->to.y()) / 2 ? above : below).push_back(int(j));
            auto width = [&](const std::vector<int>& items) {
                double w = 0;
                for (int j : items) w += bounds(doc_, ps[j]).width() + (w ? gap / 2 : 0);
                return w;
            };
            const double len = std::max({3 * kBondLength, width(above) + gap, width(below) + gap});
            Arrow& out = next.arrows[*ps[i].arrows.begin()];
            const bool leftward = a->to.x() < a->from.x();
            out.from = {leftward ? x + len : x, baseline}, out.to = {leftward ? x : x + len, baseline};
            for (const auto* items : {&above, &below}) {  // side by side, centred on the arrow
                double at = x + (len - width(*items)) / 2;
                for (int j : *items) {
                    const QRectF r = bounds(doc_, ps[j]);
                    const double y = items == &above ? baseline - 0.3 * kBondLength - r.height() / 2
                                                     : baseline + 0.3 * kBondLength + r.height() / 2;
                    shiftPiece(next, ps[j], QPointF(at - r.left(), y - r.center().y()));
                    at += r.width() + gap / 2;
                }
            }
            x += len + gap;
        } else {
            const QRectF r = bounds(doc_, ps[i]);
            shiftPiece(next, ps[i], QPointF(x - r.left(), baseline - r.center().y()));
            x += r.width() + gap;
        }
    }
    commit(next, tr("Arrange scheme"));
}

void Canvas::bracketSelection(bool square, const QString& label) {
    if (selectedAtoms_.isEmpty()) return;
    std::vector<int> atoms(selectedAtoms_.begin(), selectedAtoms_.end());
    std::sort(atoms.begin(), atoms.end());
    Document next = doc_;
    next.brackets.push_back({atoms, square, label});
    commit(next, tr("Brackets"));
}

void Canvas::removeBrackets() {
    Document next = doc_;
    std::erase_if(next.brackets, [&](const Bracket& b) {  // those around any selected atom (all, with none selected)
        return selectedAtoms_.isEmpty() ||
               std::any_of(b.atoms.begin(), b.atoms.end(), [&](int i) { return selectedAtoms_.contains(i); });
    });
    if (!(next == doc_)) commit(next, tr("Remove brackets"));
}

std::vector<int> Canvas::moleculesOfSelection() const {
    std::vector<int> out;
    std::vector<bool> seen(doc_.atoms.size());
    std::vector<int> stack;
    for (int i = 0; i < int(doc_.atoms.size()); ++i)
        if (selectedAtoms_.isEmpty() || selectedAtoms_.contains(i)) stack.push_back(i), seen[i] = true;
    while (!stack.empty()) {
        const int i = stack.back();
        stack.pop_back();
        out.push_back(i);
        for (int nb : doc_.neighbors(i))
            if (!seen[nb]) seen[nb] = true, stack.push_back(nb);
    }
    std::sort(out.begin(), out.end());
    return out;
}

void Canvas::rotate3D(double aboutX, double aboutY) {
    if (auto pose = chem::pose3D(doc_, moleculesOfSelection()))
        commit(chem::project3D(doc_, *pose, aboutX, aboutY), tr("Rotate in 3D"));
}

void Canvas::rotateSelection(double degrees) {
    if (selectedAtoms_.isEmpty() && selectedArrows_.isEmpty() && selectedTexts_.isEmpty()) return;
    transformSelection(QTransform().rotate(degrees), tr("Rotate"));
}

void Canvas::transformSelection(const QTransform& t, const QString& what) {
    QSet<int> atoms = selectedAtoms_, arrows = selectedArrows_, texts = selectedTexts_;
    if (atoms.isEmpty() && arrows.isEmpty() && texts.isEmpty()) {
        atoms = range(0, int(doc_.atoms.size())), arrows = range(0, int(doc_.arrows.size()));
        texts = range(0, int(doc_.texts.size()));
    }
    std::vector<QPointF> pts;
    for (int i : atoms) pts.push_back(doc_.atoms[i].pos);
    for (int i : arrows) pts.push_back(doc_.arrows[i].from), pts.push_back(doc_.arrows[i].to);
    for (int i : texts) pts.push_back(doc_.texts[i].pos);
    if (pts.empty()) return;
    QPointF c;
    for (QPointF p : pts) c += p / double(pts.size());
    Document next = doc_;
    applyTransform(next, atoms, arrows, texts, QTransform::fromTranslate(-c.x(), -c.y()) * t * QTransform::fromTranslate(c.x(), c.y()));
    commit(next, what);
}

QRectF Canvas::selectionBox() const {
    QPolygonF pts;
    for (int i : selectedAtoms_) pts << doc_.atoms[i].pos;
    for (int i : selectedArrows_) pts << doc_.arrows[i].from << doc_.arrows[i].to;
    for (int i : selectedTexts_) pts << doc_.texts[i].pos;
    const QRectF r = pts.boundingRect();
    if (pts.size() < 2 || (r.width() < 1e-6 && r.height() < 1e-6)) return {};
    return r.adjusted(-6, -6, 6, 6);  // clear of the atoms' labels' centres
}

std::optional<QPointF> Canvas::rotateHandle() const {
    const QRectF box = selectionBox();
    if (box.isNull()) return std::nullopt;
    return QPointF(box.center().x(), box.top() - 16 / transform().m11());
}

int Canvas::handleAt(QPointF p) const {
    const QRectF r = selectionBox();
    if (r.isNull()) return -1;
    const auto h = handlePoints(r);
    for (int k = 0; k < 8; ++k)
        if (len(h[k] - p) < 5 / transform().m11()) return k;
    return -1;
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
                                          tr("Element, group (OMe, CF3, Ph, Boc…), SMILES or any text (R, X, MgEt):"), QLineEdit::Normal,
                                          doc_.atoms[at].label.isEmpty()
                                              ? QString::fromStdString(chem::symbol(doc_.atoms[at].z))
                                              : doc_.atoms[at].label,
                                          &ok)
                        .trimmed();
    Document next = doc_;
    if (ok && !label.isEmpty() && applyLabel(next, at, label, true)) commit(next, tr("Edit label"));
}

void Canvas::editAtomProperties(int at) {
    const Atom& a = doc_.atoms[at];
    QDialog dialog(this);
    dialog.setWindowTitle(tr("Atom Properties"));
    auto* form = new QFormLayout(&dialog);
    auto* label = new QLineEdit(a.label.isEmpty() ? QString::fromStdString(chem::symbol(a.z)) : a.label);
    auto spin = [&](int lo, int hi, int value) {
        auto* s = new QSpinBox;
        s->setRange(lo, hi), s->setValue(value);
        return s;
    };
    auto* charge = spin(-8, 8, a.charge);
    auto* map = spin(0, 999, a.map);
    auto* pairs = spin(0, 3, a.lonePairs);
    auto* radicals = spin(0, 2, a.radicals);
    auto* partial = new QComboBox;
    partial->addItems({tr("none"), "δ+", "δ−"});
    partial->setCurrentIndex(a.partial > 0 ? 1 : a.partial < 0 ? 2 : 0);
    form->addRow(tr("Element or label:"), label);
    form->addRow(tr("Charge:"), charge);
    form->addRow(tr("Atom-map number (0 = none):"), map);
    form->addRow(tr("Lone pairs:"), pairs);
    form->addRow(tr("Radical electrons:"), radicals);
    form->addRow(tr("Partial charge:"), partial);
    auto* buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    form->addRow(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() != QDialog::Accepted) return;
    Document next = doc_;
    const QString was = a.label.isEmpty() ? QString::fromStdString(chem::symbol(a.z)) : a.label;
    if (label->text().trimmed() != was && !applyLabel(next, at, label->text().trimmed(), true)) return;
    Atom& out = next.atoms[at];
    out.charge = charge->value(), out.map = map->value(), out.lonePairs = pairs->value(), out.radicals = radicals->value();
    out.partial = partial->currentIndex() == 1 ? 1 : partial->currentIndex() == 2 ? -1 : 0;
    if (!(next == doc_)) commit(next, tr("Atom properties"));
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

void Canvas::colourSelection() {
    Document next = doc_;
    for (int i : selectedAtoms_) next.atoms[i].color = colour_;
    for (auto& b : next.bonds)
        if (selectedAtoms_.contains(b.a) && selectedAtoms_.contains(b.b)) b.color = colour_;
    for (int i : selectedArrows_) next.arrows[i].color = colour_;
    for (int i : selectedTexts_) next.texts[i].color = colour_;
    if (!(next == doc_)) commit(next, tr("Colour"));
}

void Canvas::editText(int i, QPointF pos) {
    // The editor uses the canvas font and tab stops, so spacing looks the same on both.
    QInputDialog dialog(this);
    dialog.setWindowTitle(tr("Text"));
    dialog.setLabelText(tr("Text (digits after letters become subscripts):"));
    dialog.setOption(QInputDialog::UsePlainTextEditForTextInput);
    dialog.setTextValue(i >= 0 ? doc_.texts[i].text : QString());
    if (auto* edit = dialog.findChild<QPlainTextEdit*>()) {
        QFont f = labelFont(documentStyle(doc_));
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
    if (arrow && (e->modifiers() & Qt::AltModifier) && (e->modifiers() & Qt::ShiftModifier)) {  // out of the page
        if (key == Qt::Key_Left || key == Qt::Key_Right) rotate3D(0, key == Qt::Key_Left ? -15 : 15);
        else rotate3D(key == Qt::Key_Up ? -15 : 15, 0);
        return;
    }
    if (arrow && (e->modifiers() & Qt::AltModifier)) {
        if (key == Qt::Key_Left || key == Qt::Key_Right) rotateSelection(key == Qt::Key_Left ? -15 : 15);
        return;
    }
    const bool selection = !selectedAtoms_.isEmpty() || !selectedArrows_.isEmpty() || !selectedTexts_.isEmpty();
    if (arrow && selection && !(e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier))) {  // nudge: 1 pt, Shift 10
        const QPointF by = arrowDirection(key) * ((e->modifiers() & Qt::ShiftModifier) ? 10 : 1);
        return transformSelection(QTransform::fromTranslate(by.x(), by.y()), tr("Nudge"));
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
            next.removeBond(hoverBond_);
            hoverBond_ = -1;
        } else {
            return;
        }
        return commit(next, tr("Delete"));
    }
    if (e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)) return QGraphicsView::keyPressEvent(e);

    const QString t = e->text();
    Document next = doc_;

    // ChemDraw: Enter takes a selected molecule to a hotspot, Space takes a hotspot to its molecule.
    if (selection && (key == Qt::Key_Return || key == Qt::Key_Enter) && !selectedAtoms_.isEmpty()) {
        hoverAtom_ = *std::min_element(selectedAtoms_.begin(), selectedAtoms_.end()), hoverBond_ = -1;
        return setSelection({});
    }
    if ((hoverAtom_ >= 0 || hoverBond_ >= 0) && t == " ") {
        const int from = hoverAtom_ >= 0 ? hoverAtom_ : doc_.bonds[hoverBond_].a;
        QSet<int> molecule{from};
        std::vector<int> stack{from};
        while (!stack.empty()) {
            const int i = stack.back();
            stack.pop_back();
            for (int nb : doc_.neighbors(i))
                if (!molecule.contains(nb)) molecule.insert(nb), stack.push_back(nb);
        }
        hoverAtom_ = hoverBond_ = -1;
        emit toolKey(" ");  // the Select tool, as ChemDraw's marquee
        return setSelection(molecule);
    }
    if (t == "g" && (hoverAtom_ >= 0 || hoverBond_ >= 0)) {  // grab: the hotspot's atom or bond, selected
        const QSet<int> atoms = hoverAtom_ >= 0 ? QSet<int>{hoverAtom_} : QSet<int>{doc_.bonds[hoverBond_].a, doc_.bonds[hoverBond_].b};
        hoverAtom_ = hoverBond_ = -1;
        emit toolKey(" ");
        return setSelection(atoms);
    }
    if (hoverAtom_ >= 0 && (t == "/" || t == "?")) return editAtomProperties(hoverAtom_);
    if (hoverAtom_ >= 0 && (key == Qt::Key_Return || key == Qt::Key_Enter || t == "=" || t == "t"))
        return editLabel(hoverAtom_);
    if (hoverAtom_ < 0 && hoverBond_ < 0) {  // no hotspot: tool keys
        static const QStringList tools{"x", "X", "j", "t", "e", " "};
        if (tools.contains(t)) return emit toolKey(t);
    }
    if (const Hotspot h = hotkey(next, {hoverAtom_, hoverBond_}, t); h.valid()) {
        commit(next, tr("Hotkey %1").arg(t));
        hoverAtom_ = h.atom, hoverBond_ = h.bond;
        viewport()->update();
        return;
    }
    QGraphicsView::keyPressEvent(e);
}

// Every entry reuses an existing edit: hotkeys for bonds and charges,
// applyLabel for atoms, the selection commands for selections.
QMenu* Canvas::contextMenuAt(QPointF at) {
    auto* menu = new QMenu(this);
    auto keyOn = [this](Hotspot h, const QString& key, const QString& what) {
        return [this, h, key, what] {
            Document next = doc_;
            if (hotkey(next, h, key).valid()) commit(next, what);
        };
    };
    auto label = [this](int atom, const QString& text) {
        return [this, atom, text] {
            Document next = doc_;
            if (applyLabel(next, atom, text)) commit(next, tr("Set atom"));
        };
    };
    const int atom = atomAt(at), bond = atom < 0 ? bondAt(at) : -1;
    const bool onSelection = (atom >= 0 && selectedAtoms_.contains(atom)) ||
                             (bond >= 0 && selectedAtoms_.contains(doc_.bonds[bond].a) &&
                              selectedAtoms_.contains(doc_.bonds[bond].b));
    if (onSelection || (atom < 0 && bond < 0 && !selectedAtoms_.isEmpty())) {
        menu->addAction(tr("Clean"), this, [this] {
            commit(chem::clean2D(doc_, {selectedAtoms_.begin(), selectedAtoms_.end()}), tr("Clean"));
        });
        menu->addAction(tr("Flip Horizontal"), this, [this] { flipSelection(true); });
        menu->addAction(tr("Flip Vertical"), this, [this] { flipSelection(false); });
        menu->addAction(tr("Rotate 90°"), this, [this] { rotateSelection(90); });
        menu->addSeparator();
        menu->addAction(tr("Copy as SMILES"), this, [this] {
            QApplication::clipboard()->setText(QString::fromStdString(chem::toSmiles(selectedSubset())));
        });
        menu->addAction(tr("Copy as InChI"), this, [this] {
            QApplication::clipboard()->setText(QString::fromStdString(chem::toInchi(selectedSubset())));
        });
        menu->addSeparator();
        menu->addAction(tr("Delete"), this, &Canvas::deleteSelection);
    } else if (atom >= 0) {
        auto* elements = menu->addMenu(tr("Element"));
        for (auto sym : {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H", "B", "Si"})
            elements->addAction(sym, this, label(atom, sym));
        auto* groups = menu->addMenu(tr("Abbreviation"));
        for (auto g : {"Me", "Et", "iPr", "tBu", "Ph", "Bn", "OMe", "OAc", "Ac", "CO2Me", "CF3", "NO2", "CN",
                       "Boc", "Cbz", "Fmoc", "Ts", "TBS", "Bpin"})
            groups->addAction(g, this, label(atom, g));
        menu->addAction(tr("Edit Label…"), this, [this, atom] { editLabel(atom); });
        menu->addAction(tr("Atom Properties… — /"), this, [this, atom] { editAtomProperties(atom); });
        if (!doc_.atoms[atom].label.isEmpty())
            menu->addAction(tr("Expand Abbreviation"), this, [this, atom] {
                Document next = doc_;
                chem::attach(next, atom, doc_.atoms[atom].label.toStdString());
                commit(next, tr("Expand"));
            });
        menu->addSeparator();
        menu->addAction(tr("Increase Charge"), this, keyOn({atom, -1}, "+", tr("Charge")));
        menu->addAction(tr("Decrease Charge"), this, keyOn({atom, -1}, "-", tr("Charge")));
        auto* marks = menu->addMenu(tr("Electrons and δ"));
        marks->addAction(tr("Add Lone Pair — :"), this, keyOn({atom, -1}, ":", tr("Lone pairs")));
        marks->addAction(tr("Radical — *"), this, keyOn({atom, -1}, "*", tr("Radical")));
        for (int d : {1, -1, 0})
            marks->addAction(d > 0 ? tr("δ+") : d < 0 ? tr("δ−") : tr("No δ"), this, [this, atom, d] {
                Document next = doc_;
                next.atoms[atom].partial = d;
                commit(next, tr("Partial charge"));
            });
        menu->addSeparator();
        menu->addAction(tr("Delete Atom"), this, [this, atom] {
            Document next = doc_;
            next.removeAtoms({atom});
            hoverAtom_ = hoverBond_ = -1;
            commit(next, tr("Delete"));
        });
    } else if (bond >= 0) {
        const Hotspot h{-1, bond};
        for (auto [text, key] : {std::pair{tr("Single"), "1"}, {tr("Double"), "2"}, {tr("Triple"), "3"}})
            menu->addAction(text, this, keyOn(h, key, tr("Change bond")));
        auto* style = menu->addMenu(tr("Style"));
        for (auto [text, key] : {std::pair{tr("Wedge"), "w"}, {tr("Hash"), "h"}, {tr("Wavy"), "y"},
                                 {tr("Bold"), "b"}, {tr("Dashed"), "d"}})
            style->addAction(text, this, keyOn(h, key, tr("Bond style")));
        if (doc_.bonds[bond].order == 2) {
            auto* side = menu->addMenu(tr("Double Bond Position"));
            for (auto [text, key] : {std::pair{tr("Left"), "l"}, {tr("Centre"), "c"}, {tr("Right"), "r"}})
                side->addAction(text, this, keyOn(h, key, tr("Bond position")));
        }
        auto* fuse = menu->addMenu(tr("Fuse Ring"));
        for (auto [text, key] : {std::pair{tr("Benzene"), "a"}, {tr("Cyclopropane"), "v"}, {tr("Cyclobutane"), "4"},
                                 {tr("Cyclopentane"), "5"}, {tr("Cyclohexane"), "6"}, {tr("Chair"), "9"}})
            fuse->addAction(text, this, keyOn(h, key, tr("Fuse ring")));
        menu->addSeparator();
        menu->addAction(tr("Delete Bond"), this, [this, bond] {
            Document next = doc_;
            next.removeBond(bond);
            hoverAtom_ = hoverBond_ = -1;
            commit(next, tr("Delete"));
        });
    } else {
        menu->addAction(tr("Select All"), this, &Canvas::selectAll);
        menu->addAction(tr("Fit to Window"), this, &Canvas::fitToDocument);
    }
    return menu;
}

void Canvas::contextMenuEvent(QContextMenuEvent* e) {
    QMenu* menu = contextMenuAt(mapToScene(e->pos()));
    menu->setAttribute(Qt::WA_DeleteOnClose);
    menu->popup(e->globalPos());
}

void Canvas::wheelEvent(QWheelEvent* e) {
    if (e->modifiers() & (Qt::ControlModifier | Qt::MetaModifier)) {
        zoomBy(std::pow(1.0015, e->angleDelta().y()));
        return;
    }
    QGraphicsView::wheelEvent(e);
}
