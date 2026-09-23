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

struct AtomInfo {
    int hydrogens = 0;
    bool valenceError = false;
};
std::vector<AtomInfo> atomInfo(const Document& doc);
}  // namespace chem
