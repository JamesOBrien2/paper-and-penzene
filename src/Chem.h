#pragma once
// The only file that talks to RDKit. Everything else sees a Document.
#include "Document.h"
#include <optional>
#include <string>

namespace chem {
std::optional<Document> fromSmiles(const std::string& smiles);
}
