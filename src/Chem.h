#pragma once
// The only file that talks to RDKit. Everything else sees a Document.
#include "Document.h"
#include <optional>
#include <string>
#include <vector>

namespace chem {
std::optional<Document> fromSmiles(const std::string& smiles);
std::optional<Document> fromMolBlock(const std::string& block);
// ChemDraw .cdxml (molecules, arrows, text) or binary .cdx (molecules only,
// where RDKit was built with ChemDraw support).
std::optional<Document> fromChemDraw(const QByteArray& data);
// .penz, .cdxml/.cdx, or MOL/SDF, by extension.
std::optional<Document> readFile(const QString& path);
std::string toMolBlock(const Document& doc);
std::string toSmiles(const Document& doc);  // "" if the structure isn't valid
// New layout, same atom order; each molecule keeps its centroid. With `only`,
// just the molecules containing those atoms are touched.
Document clean2D(const Document& doc, const std::vector<int>& only = {});

struct Properties {
    std::string formula;  // Hill order, all fragments together
    double mw = 0, exactMass = 0;
};
std::optional<Properties> properties(const Document& doc);  // nullopt if empty or invalid
std::string toInchi(const Document& doc);                   // "" if invalid
std::string toInchiKey(const Document& doc);

// Smallest set of smallest rings, each in ring order (abbreviations excluded).
std::vector<std::vector<int>> rings(const Document& doc);

// Explicit hydrogens: add them where atoms have implicit ones (placed by RDKit),
// or remove plain terminal H atoms again (wedged/hashed ones carry stereo, so stay).
Document addHydrogens(const Document& doc);
Document removeHydrogens(const Document& doc);
// CIP descriptors: (R)/(S) (or r/s) on atoms, (E)/(Z) on double bonds.
struct StereoLabel {
    int atom = -1, bond = -1;
    QString text;
};
std::vector<StereoLabel> stereoLabels(const Document& doc);

struct AtomInfo {
    int hydrogens = 0;
    bool valenceError = false;
};
std::vector<AtomInfo> atomInfo(const Document& doc);

// Abbreviations (Me, OMe, Boc…): drawn as a label, expanded for chemistry.
std::optional<Atom> abbreviationHead(const QString& label);  // attaching atom; nullopt if unknown
QStringList abbreviations();
// Replaces `atom` with the first atom of `smiles` (or an abbreviation) and lays
// the rest out away from its bonds. New atoms are appended, so indices stay valid.
bool attach(Document& doc, int atom, const std::string& smilesOrAbbreviation);
Document expanded(const Document& doc);  // abbreviations drawn out in full

std::string symbol(int z);
int atomicNumber(const std::string& symbol);  // 0 if unknown
}  // namespace chem
