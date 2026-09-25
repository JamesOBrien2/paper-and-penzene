#include "Templates.h"
#include "Chem.h"
#include "Edit.h"

#include <QtMath>

#include <QDir>
#include <QFile>
#include <QRegularExpression>
#include <QStandardPaths>

namespace {

// A drawing in bond lengths: x right, y down.
struct Sketch {
    Document d;
    int atom(double x, double y, const QString& label = "C") {
        const int i = d.addAtom(QPointF(x, y) * kBondLength);
        if (label != "C") edit::applyLabel(d, i, label);
        return i;
    }
    void bond(int a, int b, BondStereo s = BondStereo::None, int order = 1) { d.bonds.push_back({a, b, order, s}); }
};

// Haworth: the flattened ring, front edge bold; `up[k]` is ring atom k's substituent
// straight up (true) or down, "" for none. Ring atom 0 is the ring O at the back.
Document haworth(const std::vector<QPointF>& ring, const std::vector<std::pair<QString, bool>>& subs) {
    Sketch s;
    std::vector<int> ids;
    for (size_t k = 0; k < ring.size(); ++k) ids.push_back(s.atom(ring[k].x(), ring[k].y(), k == 0 ? "O" : "C"));
    for (size_t k = 0; k < ids.size(); ++k) {
        const QPointF a = ring[k], b = ring[(k + 1) % ring.size()];
        const bool front = a.y() + b.y() > 0;  // the half nearer the viewer
        s.bond(ids[k], ids[(k + 1) % ids.size()], front ? BondStereo::Bold : BondStereo::None);
    }
    for (size_t k = 0; k < subs.size(); ++k) {
        const auto& [label, up] = subs[k];
        if (label.isEmpty()) continue;
        // Up from the front edge stays short, so its label sits inside the ring.
        const double y = ring[k].y() + (up ? (ring[k].y() > 0.1 ? -0.5 : -0.75) : 0.75);
        if (label == "CH2OH") {  // drawn out, so its C sits on the vertical
            const int c = s.atom(ring[k].x(), y);
            s.bond(ids[k], c);
            s.bond(c, s.atom(ring[k].x() - 0.65, y - 0.4, "OH"));
        } else {
            s.bond(ids[k], s.atom(ring[k].x(), y, label));
        }
    }
    return s.d;
}

const std::vector<QPointF> kPyranose{{0.65, -0.45}, {1.3, 0}, {0.65, 0.45}, {-0.65, 0.45}, {-1.3, 0}, {-0.65, -0.45}};
const std::vector<QPointF> kFuranose{{0, -0.55}, {1.0, -0.1}, {0.65, 0.45}, {-0.65, 0.45}, {-1.0, -0.1}};

// Fischer: the carbon chain straight down; `right[k]` is centre k's OH on the right (H left) or left.
Document fischer(const QString& top, const std::vector<bool>& right, const QString& bottom) {
    Sketch s;
    int prev = s.atom(0, -1.0, top);
    double y = 0;
    for (bool r : right) {
        const int c = s.atom(0, y);
        s.bond(prev, c);
        s.bond(c, s.atom(r ? 1 : -1, y, "OH"));
        s.bond(c, s.atom(r ? -1 : 1, y, "H"));
        prev = c, y += 1;
    }
    s.bond(prev, s.atom(0, y, bottom));
    return s.d;
}

// Newman projections are drawings (no chemistry: the two carbons coincide).
Document newman(double backTurn, const QStringList& front, const QStringList& back) {
    Document d;
    const double r = 0.5 * kBondLength, L = 1.15 * kBondLength;
    d.arrows.push_back({{-r, -r}, {r, r}, ArrowKind::Ellipse});
    auto spoke = [&](double deg, double from, const QString& label) {
        const QPointF u(std::cos(qDegreesToRadians(deg)), -std::sin(qDegreesToRadians(deg)));
        d.arrows.push_back({u * from, u * L, ArrowKind::Line});
        const QPointF at = u * (L + 5.5);
        d.texts.push_back({at + QPointF(-2.9 * label.size(), 3.5), label});
    };
    for (int k = 0; k < 3; ++k) spoke(90 + 120 * k, 0, front[k]);
    for (int k = 0; k < 3; ++k) spoke(90 + backTurn + 120 * k, r, back[k]);
    return d;
}

}  // namespace

Document templateDocument(const Template& t) {
    if (t.draw) return t.draw();
    auto doc = chem::fromSmiles(t.smiles.toStdString());
    return doc ? *doc : Document{};
}

// Structures from PubChem (isomeric SMILES), laid out by RDKit when inserted;
// projections drawn out, with their stereo checked against PubChem's SMILES.
const std::vector<Template>& builtinTemplates() {
    static const std::vector<Template> t{
        {"Amino acids", "Glycine", "C(C(=O)O)N"},
        {"Amino acids", "L-Alanine", "C[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Valine", "CC(C)[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Leucine", "CC(C)C[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Isoleucine", "CC[C@H](C)[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Proline", "C1C[C@H](NC1)C(=O)O"},
        {"Amino acids", "L-Phenylalanine", "C1=CC=C(C=C1)C[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Tryptophan", "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Methionine", "CSCC[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Serine", "C([C@@H](C(=O)O)N)O"},
        {"Amino acids", "L-Threonine", "C[C@H]([C@@H](C(=O)O)N)O"},
        {"Amino acids", "L-Cysteine", "C([C@@H](C(=O)O)N)S"},
        {"Amino acids", "L-Tyrosine", "C1=CC(=CC=C1C[C@@H](C(=O)O)N)O"},
        {"Amino acids", "L-Asparagine", "C([C@@H](C(=O)O)N)C(=O)N"},
        {"Amino acids", "L-Glutamine", "C(CC(=O)N)[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Aspartic acid", "C([C@@H](C(=O)O)N)C(=O)O"},
        {"Amino acids", "L-Glutamic acid", "C(CC(=O)O)[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Lysine", "C(CCN)C[C@@H](C(=O)O)N"},
        {"Amino acids", "L-Arginine", "C(C[C@@H](C(=O)O)N)CN=C(N)N"},
        {"Amino acids", "L-Histidine", "C1=C(NC=N1)C[C@@H](C(=O)O)N"},
        {"Sugars", "β-D-Glucopyranose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"},
        {"Sugars", "α-D-Glucopyranose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"},
        {"Sugars", "β-D-Galactopyranose", "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O"},
        {"Sugars", "α-D-Mannopyranose", "C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O"},
        {"Sugars", "β-D-Fructofuranose", "C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O"},
        {"Sugars", "β-D-Ribofuranose", "C([C@@H]1[C@H]([C@H]([C@@H](O1)O)O)O)O"},
        {"Sugars", "2-Deoxy-β-D-ribofuranose", "C1[C@@H]([C@H](O[C@H]1O)CO)O"},
        {"Sugars", "Sucrose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O"},
        {"Nucleobases", "Adenine", "C1=NC2=NC=NC(=C2N1)N"},
        {"Nucleobases", "Guanine", "C1=NC2=C(N1)C(=O)NC(=N2)N"},
        {"Nucleobases", "Cytosine", "C1=C(NC(=O)N=C1)N"},
        {"Nucleobases", "Thymine", "CC1=CNC(=O)NC1=O"},
        {"Nucleobases", "Uracil", "C1=CNC(=O)NC1=O"},
        {"Rings and scaffolds", "Naphthalene", "C1=CC=C2C=CC=CC2=C1"},
        {"Rings and scaffolds", "Anthracene", "C1=CC=C2C=C3C=CC=CC3=CC2=C1"},
        {"Rings and scaffolds", "Phenanthrene", "C1=CC=C2C(=C1)C=CC3=CC=CC=C32"},
        {"Rings and scaffolds", "Indole", "C1=CC=C2C(=C1)C=CN2"},
        {"Rings and scaffolds", "Quinoline", "C1=CC=C2C(=C1)C=CC=N2"},
        {"Rings and scaffolds", "Isoquinoline", "C1=CC=C2C=NC=CC2=C1"},
        {"Rings and scaffolds", "Purine", "C1=C2C(=NC=N1)N=CN2"},
        {"Rings and scaffolds", "Pyrimidine", "C1=CN=CN=C1"},
        {"Rings and scaffolds", "Pyridine", "C1=CC=NC=C1"},
        {"Rings and scaffolds", "Furan", "C1=COC=C1"},
        {"Rings and scaffolds", "Thiophene", "C1=CSC=C1"},
        {"Rings and scaffolds", "Pyrrole", "C1=CNC=C1"},
        {"Rings and scaffolds", "Imidazole", "C1=CN=CN1"},
        {"Rings and scaffolds", "Morpholine", "C1COCCN1"},
        {"Rings and scaffolds", "Piperidine", "C1CCNCC1"},
        {"Rings and scaffolds", "Piperazine", "C1CNCCN1"},
        {"Rings and scaffolds", "Adamantane", "C1C2CC3CC1CC(C2)C3"},
        {"Rings and scaffolds", "Norbornane", "C1CC2CCC1C2"},
        {"Rings and scaffolds", "Bicyclo[2.2.2]octane", "C1CC2CCC1CC2"},
        {"Rings and scaffolds", "Cubane", "C12C3C4C1C5C2C3C45"},
        {"Rings and scaffolds", "Steroid core (gonane)", "C1CCC2CC[C@H]3[C@@H]4CCC[C@H]4CC[C@@H]3[C@H]2C1"},
        {"Rings and scaffolds", "Porphine", "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2"},
        {"Rings and scaffolds", "18-Crown-6", "C1COCCOCCOCCOCCOCCO1"},
        {"Rings and scaffolds", "β-Lactam (azetidin-2-one)", "C1CNC1=O"},
        {"Projections", "Haworth: β-D-glucopyranose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
         [] { return haworth(kPyranose, {{"", false}, {"OH", true}, {"OH", false}, {"OH", true}, {"OH", false}, {"CH2OH", true}}); }},
        {"Projections", "Haworth: α-D-glucopyranose", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
         [] { return haworth(kPyranose, {{"", false}, {"OH", false}, {"OH", false}, {"OH", true}, {"OH", false}, {"CH2OH", true}}); }},
        {"Projections", "Haworth: β-D-ribofuranose", "C([C@@H]1[C@H]([C@H]([C@@H](O1)O)O)O)O",
         [] { return haworth(kFuranose, {{"", false}, {"OH", true}, {"OH", false}, {"OH", false}, {"CH2OH", true}}); }},
        {"Projections", "Fischer: D-glucose", "C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O",
         [] { return fischer("CHO", {true, false, true, true}, "CH2OH"); }},
        {"Projections", "Fischer: D-glyceraldehyde", "C([C@H](C=O)O)O", [] { return fischer("CHO", {true}, "CH2OH"); }},
        {"Projections", "Newman: butane, anti (staggered)", "",
         [] { return newman(60, {"CH3", "H", "H"}, {"H", "H", "CH3"}); }},
        {"Projections", "Newman: ethane, eclipsed", "", [] { return newman(24, {"H", "H", "H"}, {"H", "H", "H"}); }},
    };
    return t;
}

static QString userDir() { return QStandardPaths::writableLocation(QStandardPaths::AppDataLocation) + "/templates"; }

std::vector<std::pair<QString, Document>> userTemplates() {
    std::vector<std::pair<QString, Document>> out;
    for (const QFileInfo& f : QDir(userDir()).entryInfoList({"*.penz"}, QDir::Files, QDir::Name)) {
        QFile file(f.filePath());
        if (!file.open(QIODevice::ReadOnly)) continue;
        if (auto doc = Document::fromJson(file.readAll()); doc && !doc->empty()) out.push_back({f.completeBaseName(), *doc});
    }
    return out;
}

bool saveUserTemplate(const QString& name, const Document& doc) {
    QString safe = name.trimmed();
    safe.replace(QRegularExpression(R"([/\\:*?"<>|])"), "-");
    if (safe.isEmpty() || doc.empty() || !QDir().mkpath(userDir())) return false;
    QFile f(userDir() + "/" + safe + ".penz");
    return f.open(QIODevice::WriteOnly) && f.write(doc.toJson()) > 0;
}

bool removeUserTemplate(const QString& name) { return QFile::remove(userDir() + "/" + name + ".penz"); }
