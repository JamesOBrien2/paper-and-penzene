#pragma once
// Structure editing without a UI: the geometry behind the drawing tools and
// every atom/bond hotkey. The canvas and (later) the Python API both use it.
#include "Document.h"

#include <QString>
#include <vector>

namespace edit {
constexpr double kMergeRadius = 0.3 * kBondLength;  // closer than this is the same atom

// Direction for a new bond from `atom` that avoids existing bonds; `newOrder`
// makes it straight on at an sp centre.
QPointF freeDirection(const Document& doc, int atom, int newOrder = 1);
QPointF snapped(QPointF from, QPointF to);  // unit vector rounded to 30°
int atomNear(const Document& doc, QPointF p, double r, int skip = -1);
int atomAtOrNew(Document& doc, QPointF p, int z = 6);
void link(Document& doc, int a, int b, int order = 1, BondStereo stereo = BondStereo::None);
std::vector<int> addRing(Document& doc, const std::vector<QPointF>& verts, bool aromatic);
void ringAt(Document& doc, QPointF centre, int n, bool aromatic);
std::vector<int> ringOnAtom(Document& doc, int atom, int n, bool aromatic);
void ringOnBond(Document& doc, int bond, int n, bool aromatic);
void chairOnBond(Document& doc, int bond, int edge);
void straightenSp(Document& doc, int bond);
// Fuses each (keep, drop) pair: drop's bonds and ring fills move to keep, then
// drop is removed. Duplicate and self bonds are dropped (the higher order wins).
void mergeAtoms(Document& doc, const std::vector<std::pair<int, int>>& keepDrop);

// Element symbol, abbreviation (drawn as its label) or SMILES (drawn out).
bool applyLabel(Document& doc, int atom, const QString& label);

// ChemDraw's hotkeys, typed with `h` as the hotspot. Returns the new hotspot,
// or an empty one ({-1, -1}) if the key means nothing there.
struct Hotspot {
    int atom = -1, bond = -1;
    bool valid() const { return atom >= 0 || bond >= 0; }
};
Hotspot hotkey(Document& doc, Hotspot h, const QString& key);
}  // namespace edit
