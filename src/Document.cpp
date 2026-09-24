#include "Document.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <algorithm>
#include <cmath>
#include <numbers>

static const char* kStereo[] = {"none", "wedge", "hash", "bold", "dashed", "wavy"};
static const char* kPosition[] = {"auto", "left", "centre", "right"};
static const char* kArrow[] = {"reaction", "equilibrium", "resonance", "retro", "fishhook"};

QByteArray Document::toJson() const {
    QJsonArray as, bs;
    for (const auto& a : atoms) {
        QJsonObject o{{"x", a.pos.x()}, {"y", a.pos.y()}, {"z", a.z}};
        if (a.charge) o["charge"] = a.charge;
        if (!a.label.isEmpty()) o["label"] = a.label;
        if (a.color.isValid()) o["color"] = a.color.name();
        as.append(o);
    }
    for (const auto& b : bonds) {
        QJsonObject o{{"a", b.a}, {"b", b.b}, {"order", b.order}};
        if (b.stereo != BondStereo::None) o["stereo"] = kStereo[int(b.stereo)];
        if (b.position != BondPosition::Auto) o["position"] = kPosition[int(b.position)];
        if (b.color.isValid()) o["color"] = b.color.name();
        bs.append(o);
    }
    QJsonObject root{{"format", "penzene"}, {"version", 1}, {"atoms", as}, {"bonds", bs}};
    QJsonArray ar, ts;
    for (const auto& a : arrows) {
        QJsonObject o{{"x1", a.from.x()}, {"y1", a.from.y()}, {"x2", a.to.x()}, {"y2", a.to.y()},
                      {"kind", kArrow[int(a.kind)]}};
        if (a.bend) o["bend"] = a.bend;
        if (a.color.isValid()) o["color"] = a.color.name();
        ar.append(o);
    }
    for (const auto& t : texts) {
        QJsonObject o{{"x", t.pos.x()}, {"y", t.pos.y()}, {"text", t.text}};
        if (t.scale != 1) o["scale"] = t.scale;
        if (t.color.isValid()) o["color"] = t.color.name();
        ts.append(o);
    }
    if (!ar.isEmpty()) root["arrows"] = ar;
    if (!ts.isEmpty()) root["texts"] = ts;
    QJsonArray fs;
    for (const auto& f : fills) {
        QJsonArray ids;
        for (int i : f.atoms) ids.append(i);
        fs.append(QJsonObject{{"atoms", ids}, {"color", f.color.name()}});
    }
    if (!fs.isEmpty()) root["fills"] = fs;
    if (!style.isEmpty()) root["style"] = style;
    if (carbonLabels != CarbonLabels::None) root["carbonLabels"] = carbonLabels == CarbonLabels::All ? "all" : "terminal";
    if (hideImplicitH) root["hideImplicitH"] = true;
    if (showStereo) root["showStereo"] = true;
    if (aromaticCircles) root["aromaticCircles"] = true;
    QJsonArray circleOverrides;
    for (const auto& ring : aromaticCircleOverrides) {
        QJsonArray ids;
        for (int i : ring) ids.append(i);
        circleOverrides.append(ids);
    }
    if (!circleOverrides.isEmpty()) root["aromaticCircleOverrides"] = circleOverrides;
    return QJsonDocument(root).toJson(QJsonDocument::Indented);
}

std::optional<Document> Document::fromJson(const QByteArray& data) {
    auto root = QJsonDocument::fromJson(data).object();
    if (root["format"].toString() != "penzene" || root["version"].toInt() != 1)
        return std::nullopt;
    Document doc;
    doc.style = root["style"].toString();
    const QString carbons = root["carbonLabels"].toString();
    doc.carbonLabels = carbons == "all"        ? Document::CarbonLabels::All
                       : carbons == "terminal" ? Document::CarbonLabels::Terminal
                                               : Document::CarbonLabels::None;
    doc.hideImplicitH = root["hideImplicitH"].toBool();
    doc.showStereo = root["showStereo"].toBool();
    doc.aromaticCircles = root["aromaticCircles"].toBool();
    for (const auto& v : root["atoms"].toArray()) {
        auto o = v.toObject();
        doc.atoms.push_back({QPointF(o["x"].toDouble(), o["y"].toDouble()),
                             o["z"].toInt(6), o["charge"].toInt(), o["label"].toString(),
                             QColor(o["color"].toString())});
    }
    const int n = int(doc.atoms.size());
    for (const auto& v : root["bonds"].toArray()) {
        auto o = v.toObject();
        Bond b{o["a"].toInt(-1), o["b"].toInt(-1), o["order"].toInt(1)};
        // Untrusted file: reject dangling or self bonds rather than crash later.
        if (b.a < 0 || b.a >= n || b.b < 0 || b.b >= n || b.a == b.b) return std::nullopt;
        b.order = std::clamp(b.order, 1, 3);
        auto index = [](const auto& names, const QString& s) {
            auto it = std::find(std::begin(names), std::end(names), s);
            return it == std::end(names) ? 0 : int(it - std::begin(names));  // unknown: default
        };
        b.stereo = BondStereo(index(kStereo, o["stereo"].toString()));
        b.position = BondPosition(index(kPosition, o["position"].toString()));
        b.color = QColor(o["color"].toString());
        doc.bonds.push_back(b);
    }
    auto finite = [](std::initializer_list<double> v) {
        return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
    };
    for (const auto& v : root["arrows"].toArray()) {
        auto o = v.toObject();
        Arrow a{{o["x1"].toDouble(), o["y1"].toDouble()}, {o["x2"].toDouble(), o["y2"].toDouble()}};
        auto k = std::find(std::begin(kArrow), std::end(kArrow), o["kind"].toString("reaction"));
        a.bend = o["bend"].toDouble();
        if (k == std::end(kArrow) || !finite({a.from.x(), a.from.y(), a.to.x(), a.to.y(), a.bend}))
            return std::nullopt;
        a.kind = ArrowKind(k - std::begin(kArrow));
        a.color = QColor(o["color"].toString());
        doc.arrows.push_back(a);
    }
    for (const auto& v : root["texts"].toArray()) {
        auto o = v.toObject();
        Text t{{o["x"].toDouble(), o["y"].toDouble()}, o["text"].toString(), o["scale"].toDouble(1),
               QColor(o["color"].toString())};
        if (!finite({t.pos.x(), t.pos.y(), t.scale}) || t.scale <= 0) return std::nullopt;
        doc.texts.push_back(t);
    }
    for (const auto& v : root["fills"].toArray()) {
        auto o = v.toObject();
        Fill f{{}, QColor(o["color"].toString())};
        for (const auto& i : o["atoms"].toArray()) f.atoms.push_back(i.toInt(-1));
        if (!f.color.isValid() || f.atoms.size() < 3 ||
            std::any_of(f.atoms.begin(), f.atoms.end(), [n](int i) { return i < 0 || i >= n; }))
            return std::nullopt;
        doc.fills.push_back(f);
    }
    for (const auto& v : root["aromaticCircleOverrides"].toArray()) {
        std::vector<int> ring;
        for (const auto& i : v.toArray()) ring.push_back(i.toInt(-1));
        std::sort(ring.begin(), ring.end());
        if (ring.size() < 3 || ring.front() < 0 || ring.back() >= n ||
            std::adjacent_find(ring.begin(), ring.end()) != ring.end()) return std::nullopt;
        doc.aromaticCircleOverrides.push_back(std::move(ring));
    }
    return doc;
}

void Document::append(const Document& o, QPointF shift) {
    const int base = int(atoms.size());
    for (auto a : o.atoms) a.pos += shift, atoms.push_back(a);
    for (auto b : o.bonds) b.a += base, b.b += base, bonds.push_back(b);
    for (auto a : o.arrows) a.from += shift, a.to += shift, arrows.push_back(a);
    for (auto t : o.texts) t.pos += shift, texts.push_back(t);
    for (auto f : o.fills) {
        for (int& i : f.atoms) i += base;
        fills.push_back(f);
    }
    for (auto ring : o.aromaticCircleOverrides) {
        for (int& i : ring) i += base;
        aromaticCircleOverrides.push_back(std::move(ring));
    }
}

int Document::addAtom(QPointF pos, int z) {
    atoms.push_back({pos, z});
    return int(atoms.size()) - 1;
}

int Document::bondBetween(int a, int b) const {
    for (size_t i = 0; i < bonds.size(); ++i)
        if ((bonds[i].a == a && bonds[i].b == b) || (bonds[i].a == b && bonds[i].b == a))
            return int(i);
    return -1;
}

std::vector<int> Document::neighbors(int atom) const {
    std::vector<int> out;
    for (const auto& b : bonds)
        if (b.a == atom) out.push_back(b.b);
        else if (b.b == atom) out.push_back(b.a);
    return out;
}

void Document::removeBond(int bond) {
    const int a = bonds[bond].a, b = bonds[bond].b;
    const bool dropA = neighbors(a).size() == 1;
    const bool dropB = neighbors(b).size() == 1;
    bonds.erase(bonds.begin() + bond);
    std::vector<int> drop;
    if (dropA) drop.push_back(a);
    if (dropB) drop.push_back(b);
    if (!drop.empty()) removeAtoms(drop);
}

void Document::removeAtoms(const std::vector<int>& drop) {
    std::vector<int> remap(atoms.size(), 0);
    for (int i : drop) remap[i] = -1;
    std::vector<Atom> kept;
    for (size_t i = 0; i < atoms.size(); ++i)
        if (remap[i] != -1) remap[i] = int(kept.size()), kept.push_back(atoms[i]);
    atoms = std::move(kept);
    std::erase_if(bonds, [&](const Bond& b) { return remap[b.a] < 0 || remap[b.b] < 0; });
    for (auto& b : bonds) b.a = remap[b.a], b.b = remap[b.b];
    std::erase_if(fills, [&](const Fill& f) {
        return std::any_of(f.atoms.begin(), f.atoms.end(), [&](int i) { return remap[i] < 0; });
    });
    for (auto& f : fills)
        for (int& i : f.atoms) i = remap[i];
    std::erase_if(aromaticCircleOverrides, [&](const auto& ring) {
        return std::any_of(ring.begin(), ring.end(), [&](int i) { return remap[i] < 0; });
    });
    for (auto& ring : aromaticCircleOverrides)
        for (int& i : ring) i = remap[i];
}

// Direction pointing away from all of the atom's bonds: the bisector of the
// widest gap between them (straight on for a terminal atom).
QPointF Document::awayDirection(int atom) const {
    auto nbs = neighbors(atom);
    if (nbs.empty()) return {0, -1};
    QPointF p = atoms[atom].pos;
    std::vector<double> ang;
    for (int nb : nbs) ang.push_back(std::atan2(atoms[nb].pos.y() - p.y(), atoms[nb].pos.x() - p.x()));
    std::sort(ang.begin(), ang.end());
    double bestGap = -1, bestMid = 0;
    for (size_t i = 0; i < ang.size(); ++i) {
        double next = i + 1 < ang.size() ? ang[i + 1] : ang[0] + 2 * std::numbers::pi;
        if (next - ang[i] > bestGap + 1e-6) bestGap = next - ang[i], bestMid = (ang[i] + next) / 2;
    }
    return {std::cos(bestMid), std::sin(bestMid)};
}
