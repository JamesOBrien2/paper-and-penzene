#include "Chem.h"

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <memory>

namespace chem {

// RDKit depicts with 1.5 Å bonds and y pointing up; scenes use points, y down.
constexpr double kScale = kBondLength / 1.5;

static Document fromRDKit(RDKit::RWMol& mol) {
    RDKit::MolOps::Kekulize(mol, true);  // draw explicit double bonds
    const auto& conf = mol.getConformer();
    Document doc;
    for (const auto* a : mol.atoms()) {
        const auto& p = conf.getAtomPos(a->getIdx());
        doc.atoms.push_back({QPointF(p.x * kScale, -p.y * kScale),
                             int(a->getAtomicNum()), a->getFormalCharge()});
    }
    for (const auto* b : mol.bonds()) {
        int order = b->getBondType() == RDKit::Bond::DOUBLE   ? 2
                    : b->getBondType() == RDKit::Bond::TRIPLE ? 3
                                                              : 1;
        doc.bonds.push_back({int(b->getBeginAtomIdx()), int(b->getEndAtomIdx()), order});
    }
    return doc;
}

std::optional<Document> fromSmiles(const std::string& smiles) {
    std::unique_ptr<RDKit::RWMol> mol;
    try {
        mol.reset(RDKit::SmilesToMol(smiles));
    } catch (...) {
        return std::nullopt;
    }
    if (!mol) return std::nullopt;
    RDDepict::compute2DCoords(*mol);
    return fromRDKit(*mol);
}

}  // namespace chem
