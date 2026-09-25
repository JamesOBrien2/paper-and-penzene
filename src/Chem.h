#pragma once
// The only file that talks to RDKit. Everything else sees a Document.
#include "Document.h"
#include <optional>
#include <string>
#include <vector>

namespace chem {
std::optional<Document> fromSmiles(const std::string& smiles);
std::optional<Document> fromMolBlock(const std::string& block);
std::optional<Document> fromInchi(const std::string& inchi);
// ChemDraw .cdxml (molecules, arrows, text) or binary .cdx (molecules only,
// where RDKit was built with ChemDraw support).
std::optional<Document> fromChemDraw(const QByteArray& data);
// .penz, .cdxml/.cdx, a Penzene SVG/PNG, MOL, or every record of an SDF,
// .smi or .inchi file laid out as a grid, by extension.
std::optional<Document> readFile(const QString& path);
// Each record of a multi-record file (SDF, .smi "SMILES name" lines, .inchi
// lines), named by the file (its title or name column) or as base-N.
struct Record {
    QString name;
    std::optional<Document> doc;
};
std::vector<Record> readRecords(const QString& path);
std::string toMolBlock(const Document& doc, bool v3000 = false);
QByteArray toCdxml(const Document& doc);  // molecules, arrows and text
QByteArray toCdx(const Document& doc);    // binary CDX: molecules only
std::string toSmiles(const Document& doc);  // "" if the structure isn't valid

// A drawn reaction: the molecules before, alongside and after its arrow.
struct Reaction {
    std::vector<Document> reactants, agents, products;
};
std::optional<Reaction> reactionOf(const Document& doc);  // nullopt without a reaction arrow
std::string toReactionSmiles(const Reaction& r);           // reactants>agents>products
std::string toRxn(const Reaction& r);                      // MDL Rxnfile (V2000)
Document layoutReaction(const Reaction& r);
std::optional<Document> fromReactionSmiles(const std::string& smiles);
std::optional<Document> fromRxn(const std::string& text);
// New layout, same atom order; each molecule keeps its centroid. With `only`,
// just the molecules containing those atoms are touched.
Document clean2D(const Document& doc, const std::vector<int>& only = {});

struct Properties {
    std::string formula;  // Hill order, all fragments together
    double mw = 0, exactMass = 0;
};
std::optional<Properties> properties(const Document& doc);  // nullopt if empty or invalid
// The properties panel: everything an SI table or a med-chem glance needs.
struct Profile {
    Properties basic;
    double logP = 0, tpsa = 0;
    int hbd = 0, hba = 0, rotatable = 0, heavyAtoms = 0;
    std::vector<std::pair<std::string, double>> elemental;  // symbol, mass %, in Hill order
    int lipinskiViolations = 0;  // Ro5: MW > 500, logP > 5, HBD > 5, HBA > 10
    bool veber = true;           // rotatable bonds <= 10 and TPSA <= 140
};
std::optional<Profile> profile(const Document& doc);  // nullopt if empty or invalid
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
// Things a chemist would want fixed before publishing: valence errors,
// stereocentres without a wedge, wedges on non-stereocentres, unknown labels,
// overlapping atoms. Each names the atoms involved, to select them.
struct Problem {
    QString message;
    std::vector<int> atoms;
};
std::vector<Problem> checkStructure(const Document& doc);
// Rings RDKit perceives as aromatic (every atom aromatic), in ring order.
std::vector<std::vector<int>> aromaticRings(const Document& doc);

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
// Fischer crossings and Haworth rings redrawn with the wedges they mean (chemistry uses this).
Document projectionsAsWedges(const Document& doc);

std::string symbol(int z);
int atomicNumber(const std::string& symbol);  // 0 if unknown
}  // namespace chem
