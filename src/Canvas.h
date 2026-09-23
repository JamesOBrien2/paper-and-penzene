#pragma once
#include "Document.h"

#include <QGraphicsView>
#include <QPainterPath>
#include <QPicture>
#include <QSet>
#include <vector>

class QPainter;
class QUndoStack;

struct RenderStyle {
    QColor ink = Qt::black;
    QColor error = QColor(220, 40, 40);
};

// Paints a document with ACS 1996 proportions. Shared by the canvas and export.
void paintDocument(QPainter& p, const Document& doc, const RenderStyle& style = {});
QRectF documentBounds(const Document& doc);
// Writes .svg, .png or .pdf (by extension), cropped to the drawing.
bool exportDocument(const Document& doc, const QString& path);
QImage renderImage(const Document& doc, double dpi = 300);
QByteArray renderSvg(const Document& doc);
QPainterPath arrowPath(const Arrow& a);
QPainterPath textPath(const Text& t);

class Canvas : public QGraphicsView {
    Q_OBJECT
public:
    enum class Tool { Select, Atom, Bond, Wedge, Hash, Chain, Ring, ChargePlus, ChargeMinus, Erase, Arrow, Text };

    explicit Canvas(QUndoStack* undo, QWidget* parent = nullptr);

    const Document& document() const { return doc_; }
    // Every edit goes through here so it can be undone.
    void commit(const Document& next, const QString& text);
    void setDocumentSilently(const Document& doc);  // e.g. undo/redo, file open

    const QSet<int>& selection() const { return selectedAtoms_; }
    const QSet<int>& selectedArrows() const { return selectedArrows_; }
    const QSet<int>& selectedTexts() const { return selectedTexts_; }
    void setSelection(QSet<int> atoms, QSet<int> arrows = {}, QSet<int> texts = {});
    Document selectedSubset() const;  // selection (or everything) as a standalone doc
    void deleteSelection();
    void insert(Document fragment, const QString& text);  // centred in view, selected
    void selectAll();
    void rotateSelection(double degrees);
    void moveHotspot(QPointF dir, bool jump);
    void editLabel(int atom);
    void expandAbbreviations();  // selection, else hotspot atom, else everything
    void editText(int text, QPointF pos = {});  // text < 0: new text at pos
    int hotspotAtom() const { return hoverAtom_; }
    int hotspotBond() const { return hoverBond_; }
    // Element symbol, group (OMe, CF3, Ph…) or SMILES; first atom replaces `atom`.
    static bool applyLabel(Document& doc, int atom, const QString& label);
    QPointF viewCenter() const;
    void zoomBy(double factor);
    void fitToDocument();

    void setTool(Tool t) { tool_ = t; }
    void setElement(int z) { element_ = z; }
    void setBondOrder(int order) { bondOrder_ = order; }
    void setRing(int size, bool aromatic) { ringSize_ = size, ringAromatic_ = aromatic; }
    void setArrow(ArrowKind kind, bool curved) { arrowKind_ = kind, arrowCurved_ = curved; }

signals:
    void documentChanged();

protected:
    void drawBackground(QPainter* p, const QRectF& rect) override;
    void drawForeground(QPainter* p, const QRectF& rect) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;

private:
    int atomAt(QPointF p) const;
    int bondAt(QPointF p) const;
    int arrowAt(QPointF p) const;
    int textAt(QPointF p) const;
    Arrow draggedArrow() const;
    void refresh();
    std::vector<QPointF> dragPath() const;

    Document doc_;
    QPicture picture_;
    std::vector<QPointF> preview_;
    QUndoStack* undo_;
    Tool tool_ = Tool::Bond;
    int element_ = 6, bondOrder_ = 1, ringSize_ = 6;
    bool ringAromatic_ = true;
    ArrowKind arrowKind_ = ArrowKind::Reaction;
    bool arrowCurved_ = false;

    QSet<int> selectedAtoms_, selectedArrows_, selectedTexts_;
    int hoverAtom_ = -1, hoverBond_ = -1;

    // Drag state
    enum class Drag { None, Bond, Chain, Arrow, Move, Rotate, Rubber, Pan } drag_ = Drag::None;
    QPointF pressPos_, curPos_;
    int pressAtom_ = -1;
    Document beforeDrag_;
    QPoint panLast_;
};
