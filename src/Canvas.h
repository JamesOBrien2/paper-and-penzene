#pragma once
#include "Render.h"

#include <QFont>
#include <QGraphicsView>
#include <QPainterPath>
#include <QPicture>
#include <QSet>
#include <vector>

class QPainter;
class QMenu;
class QUndoStack;

class Canvas : public QGraphicsView {
    Q_OBJECT
public:
    enum class Tool { Select, Atom, Bond, Wedge, Hash, Chain, Ring, ChargePlus, ChargeMinus, Erase, Arrow, Text, Fill, Colour };

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
    // Brackets around the selected atoms (square or round, with a subscript such as "n"); none removes theirs.
    void bracketSelection(bool square, const QString& label);
    void removeBrackets();
    // `t` (rotate, scale, stretch) about the selection's centre; everything if nothing is selected.
    void transformSelection(const QTransform& t, const QString& what);
    void duplicateSelection(QPointF dir);
    void flipSelection(bool horizontal);  // mirror image, as ChemDraw's flip
    enum class Align { Left, HCentre, Right, Top, VCentre, Bottom };
    void alignSelection(Align edge);
    void distributeSelection(bool horizontal);
    // The right-click menu for whatever is at `scenePos` (public so tests can inspect it).
    QMenu* contextMenuAt(QPointF scenePos);
    void moveHotspot(QPointF dir, bool jump);
    void editLabel(int atom);
    void expandAbbreviations();  // selection, else hotspot atom, else everything
    void editText(int text, QPointF pos = {});  // text < 0: new text at pos
    int hotspotAtom() const { return hoverAtom_; }
    int hotspotBond() const { return hoverBond_; }
    QPointF viewCenter() const;
    void zoomBy(double factor);
    void fitToDocument();

    void setTool(Tool t) { tool_ = t; }
    void setElement(int z) { element_ = z; }
    // The Bond tool's order and style (None, Interaction or Partial).
    void setBondOrder(int order, BondStereo style = BondStereo::None) { bondOrder_ = order, bondStyle_ = style; }
    void setRing(int size, bool aromatic) { ringSize_ = size, ringAromatic_ = aromatic; }
    void setArrow(ArrowKind kind, bool curved) { arrowKind_ = kind, arrowCurved_ = curved; }
    void setTheme(const Theme& t) { theme_ = t, refresh(); }
    void setFillColor(QColor c) { fillColor_ = c; }
    QColor fillColor() const { return fillColor_; }
    void setColour(QColor c) { colour_ = c; }
    QColor colour() const { return colour_; }
    void colourSelection();  // the current colour on the selected atoms, bonds, arrows and text

signals:
    void documentChanged();
    void selectionChanged();
    void toolKey(const QString& key);  // x bond, X chain, j benzene, t text, e arrow, space select

protected:
    void drawBackground(QPainter* p, const QRectF& rect) override;
    void drawForeground(QPainter* p, const QRectF& rect) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;

private:
    int atomAt(QPointF p) const;
    int bondAt(QPointF p) const;
    int arrowAt(QPointF p) const;
    int textAt(QPointF p) const;
    Arrow draggedArrow() const;
    int draggedRingSize() const;
    void addDraggedRing(Document& doc) const;
    void refresh();
    std::vector<QPointF> dragPath() const;

    Document doc_;
    QPicture picture_;
    Theme theme_;
    std::vector<QPointF> preview_;
    QUndoStack* undo_;
    Tool tool_ = Tool::Bond;
    int element_ = 6, bondOrder_ = 1, ringSize_ = 6;
    BondStereo bondStyle_ = BondStereo::None;
    bool ringAromatic_ = true;
    ArrowKind arrowKind_ = ArrowKind::Reaction;
    bool arrowCurved_ = false;
    QColor fillColor_ = QColor(207, 227, 255);
    QColor colour_ = QColor(214, 39, 40);

    QSet<int> selectedAtoms_, selectedArrows_, selectedTexts_;
    int hoverAtom_ = -1, hoverBond_ = -1;

    // Drag state
    enum class Drag { None, Bond, Chain, Arrow, Ring, Move, Rotate, Rubber, Pan, Scale } drag_ = Drag::None;
    // Scale handles around the selection: corners scale, edges stretch along one axis.
    QRectF selectionBox() const;  // empty unless something with extent is selected
    int handleAt(QPointF p) const;  // 0..7 clockwise from the top-left corner, or -1
    int scaleHandle_ = -1;
    QRectF scaleBox_;
    QPointF pressPos_, curPos_;
    bool shift_ = false;  // held during the drag: free bond angle, or move along one axis
    int pressAtom_ = -1;
    Document beforeDrag_;
    QPoint panLast_;
};
