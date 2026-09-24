#pragma once
// Drawing a Document: styles, themes, painting and export. No widgets, so the
// command line and (later) the Python API render exactly like the app.
#include "Document.h"

#include <QColor>
#include <QFont>
#include <QImage>
#include <QPainterPath>
#include <vector>

class QPainter;

struct RenderStyle {
    QColor ink = Qt::black;
    QColor error = QColor(220, 40, 40);
    double lineWidth = 0;  // > 0 overrides the drawing style's (toolbar icons)
};

// Screen colours. Exports and copies always use black ink on a clear
// background, so a dark theme never leaks into a paper.
struct Theme {
    QString name;
    bool dark = false;
    QColor paper = Qt::white, ink = Qt::black, error = QColor(220, 40, 40);
    QColor accent = QColor(40, 120, 255), hotspot = QColor(40, 170, 60);
    QColor window, surface, text;  // UI palette; invalid = leave Qt's own
};
const std::vector<Theme>& themes();  // "System" first; "System" follows the OS light/dark
const Theme& theme(const QString& name);

// A document style preset, like ChemDraw stationery. Lengths in points.
// Geometry is always drawn with 14.4 pt bonds, so every length here is in
// those model units; exports scale by bondLength / 14.4 to the style's size.
struct DrawingStyle {
    QString name;
    double bondLength;  // the style's own bond length, in points
    double lineWidth, boldWidth, wedgeWidth, hashSpacing;
    double bondSpacing;  // double-bond gap, as a fraction of the bond length
    double labelRadius;  // bonds stop this short of a label's centre
    QString font;
    QFont::Weight weight;
    double fontSize;
};
const std::vector<DrawingStyle>& drawingStyles();  // ACS 1996 first
const DrawingStyle& drawingStyle(const QString& name);  // unknown or empty: ACS 1996

// Paints a document in its drawing style. Shared by the canvas and export.
void paintDocument(QPainter& p, const Document& doc, const RenderStyle& style = {});
QRectF documentBounds(const Document& doc);
double exportScale(const Document& doc);  // points per model unit in exports
// Writes .svg, .png or .pdf (by extension), cropped to the drawing.
bool exportDocument(const Document& doc, const QString& path);
QImage renderImage(const Document& doc, double dpi = 300);
QByteArray renderSvg(const Document& doc);
QPainterPath arrowPath(const Arrow& a);
QFont labelFont(const DrawingStyle& s, double scale = 1);
constexpr int kTabSpaces = 8;  // text tab stops, in spaces: the canvas and the text dialog agree
QPainterPath textPath(const Text& t, const DrawingStyle& s = drawingStyles()[0]);
