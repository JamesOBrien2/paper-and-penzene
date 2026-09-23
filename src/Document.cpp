#include "Document.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <algorithm>

static const char* kStereo[] = {"none", "wedge", "hash"};

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
    return doc;
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
