#include "Edit.h"
#include "Chem.h"
#include "Geometry.h"

#include <QHash>
#include <QRegularExpression>
#include <algorithm>

namespace edit {


// Direction for a new bond from `atom` that avoids existing bonds.
// `newOrder` is the order of the bond about to be added: it makes the atom sp
// (straight on) after a triple bond, or when it cumulates two double bonds.
QPointF freeDirection(const Document& doc, int atom, int newOrder) {
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

int atomNear(const Document& doc, QPointF p, double r, int skip) {
    int best = -1;
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        double d = len(doc.atoms[i].pos - p);
        if (int(i) != skip && d < r) r = d, best = int(i);
    }
    return best;
}

// Returns the atom at `p`, creating one if nothing is close enough.
int atomAtOrNew(Document& doc, QPointF p, int z) {
    int i = atomNear(doc, p, kMergeRadius);
    return i >= 0 ? i : doc.addAtom(p, z);
}

void link(Document& doc, int a, int b, int order, BondStereo stereo) {
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

}  // namespace

// Atom label hotkeys (the hotspot atom becomes this element or group).
QString labelHotkey(const QString& key) {
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
bool applyLabel(Document& doc, int at, const QString& label) {
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


Hotspot hotkey(Document& doc, Hotspot h, const QString& t) {
    if (h.atom >= 0) {
        const int at = h.atom;
        if (t == "+" || t == "-") {
            doc.atoms[at].charge += t == "+" ? 1 : -1;
            return h;
        }
        int hot = sproutHotkey(doc, at, t);
        if (hot == kUnhandled) {
            QString label = labelHotkey(t);
            if (label.isEmpty() || !applyLabel(doc, at, label)) return {};
            hot = at;
        }
        return {hot, -1};
    }
    if (h.bond < 0) return {};
    Bond& b = doc.bonds[h.bond];
    static const QHash<QString, std::pair<int, bool>> fuse{
        {"a", {6, true}}, {"z", {5, true}}, {"v", {3, false}}, {"4", {4, false}},
        {"5", {5, false}}, {"6", {6, false}}, {"7", {7, false}}, {"8", {8, false}}};
    if (t == "2" && b.order == 2 && b.stereo == BondStereo::None) {
        // Already double: move the second line to the other side (centred goes to one side).
        b.position = doubleBondSide(doc, b) > 0 ? BondPosition::Left : BondPosition::Right;
    } else if (t == "1" || t == "2" || t == "3") {
        b.order = t.toInt(), b.stereo = BondStereo::None, b.position = BondPosition::Auto;
        straightenSp(doc, h.bond);
    } else if (t == "w" || t == "h" || t == "H") {
        BondStereo s = t == "w" ? BondStereo::Wedge : BondStereo::Hash;
        if (b.stereo == s) std::swap(b.a, b.b);  // again: flip which end is narrow
        b.stereo = s, b.order = 1;
    } else if (fuse.contains(t)) {
        ringOnBond(doc, h.bond, fuse[t].first, fuse[t].second);
    } else if (t == "9" || t == "0") {
        chairOnBond(doc, h.bond, t == "9" ? 0 : 1);
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
        return {};
    }
    return h;
}

}  // namespace edit
