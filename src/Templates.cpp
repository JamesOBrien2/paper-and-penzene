#include "Templates.h"

#include <QDir>
#include <QFile>
#include <QRegularExpression>
#include <QStandardPaths>

// Structures from PubChem (isomeric SMILES), laid out by RDKit when inserted.
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
