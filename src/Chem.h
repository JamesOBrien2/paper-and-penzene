#pragma once
// The only file that talks to RDKit. Everything else sees a Document.
#include "Document.h"
#include <optional>
#include <string>
#include <vector>

namespace chem {
std::optional<Document> fromSmiles(const std::string& smiles);
std::optional<Document> fromMolBlock(const std::string& block);
std::string toMolBlock(const Document& doc);
std::string toSmiles(const Document& doc);  // "" if the structure isn't valid
Document clean2D(const Document& doc);      // new layout, same atom order and centroid

struct Properties {
    std::string formula;  // Hill order, all fragments together
    double mw = 0, exactMass = 0;
};
std::optional<Properties> properties(const Document& doc);  // nullopt if empty or invalid
std::string toInchi(const Document& doc);                   // "" if invalid
std::string toInchiKey(const Document& doc);

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
