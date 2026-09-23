#include "Document.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <algorithm>
#include <cmath>

static const char* kStereo[] = {"none", "wedge", "hash"};
static const char* kArrow[] = {"reaction", "equilibrium", "resonance", "retro", "fishhook"};

QByteArray Document::toJson() const {
    QJsonArray as, bs;
    for (const auto& a : atoms) {
        QJsonObject o{{"x", a.pos.x()}, {"y", a.pos.y()}, {"z", a.z}};
        if (a.charge) o["charge"] = a.charge;
        as.append(o);
    }
    for (const auto& b : bonds) {
        QJsonObject o{{"a", b.a}, {"b", b.b}, {"order", b.order}};
        if (b.stereo != BondStereo::None) o["stereo"] = kStereo[int(b.stereo)];
        bs.append(o);
    }
    QJsonObject root{{"format", "penzene"}, {"version", 1}, {"atoms", as}, {"bonds", bs}};
    QJsonArray ar, ts;
    for (const auto& a : arrows) {
        QJsonObject o{{"x1", a.from.x()}, {"y1", a.from.y()}, {"x2", a.to.x()}, {"y2", a.to.y()},
                      {"kind", kArrow[int(a.kind)]}};
        if (a.bend) o["bend"] = a.bend;
        ar.append(o);
    }
    for (const auto& t : texts) ts.append(QJsonObject{{"x", t.pos.x()}, {"y", t.pos.y()}, {"text", t.text}});
    if (!ar.isEmpty()) root["arrows"] = ar;
    if (!ts.isEmpty()) root["texts"] = ts;
    return QJsonDocument(root).toJson(QJsonDocument::Indented);
}

std::optional<Document> Document::fromJson(const QByteArray& data) {
    auto root = QJsonDocument::fromJson(data).object();
    if (root["format"].toString() != "penzene" || root["version"].toInt() != 1)
        return std::nullopt;
    Document doc;
    for (const auto& v : root["atoms"].toArray()) {
        auto o = v.toObject();
        doc.atoms.push_back({QPointF(o["x"].toDouble(), o["y"].toDouble()),
                             o["z"].toInt(6), o["charge"].toInt()});
    }
    const int n = int(doc.atoms.size());
    for (const auto& v : root["bonds"].toArray()) {
        auto o = v.toObject();
        Bond b{o["a"].toInt(-1), o["b"].toInt(-1), o["order"].toInt(1)};
        // Untrusted file: reject dangling or self bonds rather than crash later.
        if (b.a < 0 || b.a >= n || b.b < 0 || b.b >= n || b.a == b.b) return std::nullopt;
        b.order = std::clamp(b.order, 1, 3);
        auto s = o["stereo"].toString();
        b.stereo = s == "wedge" ? BondStereo::Wedge : s == "hash" ? BondStereo::Hash : BondStereo::None;
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
        doc.arrows.push_back(a);
    }
    for (const auto& v : root["texts"].toArray()) {
        auto o = v.toObject();
        Text t{{o["x"].toDouble(), o["y"].toDouble()}, o["text"].toString()};
        if (!finite({t.pos.x(), t.pos.y()})) return std::nullopt;
        doc.texts.push_back(t);
    }
    return doc;
}

void Document::append(const Document& o, QPointF shift) {
    const int base = int(atoms.size());
    for (auto a : o.atoms) a.pos += shift, atoms.push_back(a);
    for (auto b : o.bonds) b.a += base, b.b += base, bonds.push_back(b);
    for (auto a : o.arrows) a.from += shift, a.to += shift, arrows.push_back(a);
    for (auto t : o.texts) t.pos += shift, texts.push_back(t);
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

void Document::removeAtoms(const std::vector<int>& drop) {
    std::vector<int> remap(atoms.size(), 0);
    for (int i : drop) remap[i] = -1;
    std::vector<Atom> kept;
    for (size_t i = 0; i < atoms.size(); ++i)
        if (remap[i] != -1) remap[i] = int(kept.size()), kept.push_back(atoms[i]);
    atoms = std::move(kept);
    std::erase_if(bonds, [&](const Bond& b) { return remap[b.a] < 0 || remap[b.b] < 0; });
    for (auto& b : bonds) b.a = remap[b.a], b.b = remap[b.b];
}
