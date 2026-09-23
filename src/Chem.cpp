#include "Chem.h"

#include <GraphMol/Chirality.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <cmath>
#include <memory>

namespace chem {
using RDKit::RWMol;

// MOL files and RDKit use 1.5 Å bonds with y up; scenes use points with y down.
constexpr double kScale = kBondLength / 1.5;

static std::unique_ptr<RWMol> toRDKit(const Document& doc) {
    auto mol = std::make_unique<RWMol>();
    auto* conf = new RDKit::Conformer(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        auto* atom = new RDKit::Atom(a.z);
        atom->setFormalCharge(a.charge);
        mol->addAtom(atom, true, true);
        conf->setAtomPos(i, {a.pos.x() / kScale, -a.pos.y() / kScale, 0});
    }
    for (const auto& b : doc.bonds) {
        auto type = b.order == 2 ? RDKit::Bond::DOUBLE
                    : b.order == 3 ? RDKit::Bond::TRIPLE
                                   : RDKit::Bond::SINGLE;
        mol->addBond(b.a, b.b, type);
        if (b.stereo != BondStereo::None) {
            auto* bond = mol->getBondBetweenAtoms(b.a, b.b);
            bond->setBondDir(b.stereo == BondStereo::Wedge ? RDKit::Bond::BEGINWEDGE
                                                           : RDKit::Bond::BEGINDASH);
            bond->setProp(RDKit::common_properties::_MolFileBondStereo,
                          b.stereo == BondStereo::Wedge ? 1u : 6u);
        }
    }
    conf->set3D(false);
    mol->addConformer(conf, true);
    return mol;
}

// Best-effort sanitize + stereo from wedges. Drawings in progress are often
// invalid (pentavalent carbon), so failure just leaves a partly perceived mol.
static bool perceive(RWMol& mol) {
    try {
        RDKit::MolOps::sanitizeMol(mol);
    } catch (...) {
        mol.updatePropertyCache(false);
        RDKit::MolOps::fastFindRings(mol);
        return false;
    }
    RDKit::MolOps::assignChiralTypesFromBondDirs(mol);
    RDKit::MolOps::assignStereochemistry(mol, true, true);
    return true;
}

static Document fromRDKit(RWMol& mol) {
    try {
        RDKit::MolOps::Kekulize(mol, true);  // draw explicit double bonds
    } catch (...) {
    }
    const auto& conf = mol.getConformer();
    // Normalise whatever bond length the source used to ours.
    double sum = 0;
    for (const auto* b : mol.bonds())
        sum += (conf.getAtomPos(b->getBeginAtomIdx()) - conf.getAtomPos(b->getEndAtomIdx())).length();
    double scale = mol.getNumBonds() && sum > 1e-6 ? kBondLength * mol.getNumBonds() / sum : kScale;

    Document doc;
    for (const auto* a : mol.atoms()) {
        const auto& p = conf.getAtomPos(a->getIdx());
        doc.atoms.push_back({QPointF(p.x * scale, -p.y * scale), int(a->getAtomicNum()),
                             a->getFormalCharge()});
    }
    for (const auto* b : mol.bonds()) {
        Bond out{int(b->getBeginAtomIdx()), int(b->getEndAtomIdx())};
        out.order = b->getBondType() == RDKit::Bond::DOUBLE   ? 2
                    : b->getBondType() == RDKit::Bond::TRIPLE ? 3
                                                              : 1;
        out.stereo = b->getBondDir() == RDKit::Bond::BEGINWEDGE  ? BondStereo::Wedge
                     : b->getBondDir() == RDKit::Bond::BEGINDASH ? BondStereo::Hash
                                                                 : BondStereo::None;
        doc.bonds.push_back(out);
    }
    return doc;
}

std::optional<Document> fromSmiles(const std::string& smiles) {
    std::unique_ptr<RWMol> mol;
    try {
        mol.reset(RDKit::SmilesToMol(smiles));
    } catch (...) {
        return std::nullopt;
    }
    if (!mol) return std::nullopt;
    RDDepict::preferCoordGen = true;
    RDDepict::compute2DCoords(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    return fromRDKit(*mol);
}

std::optional<Document> fromMolBlock(const std::string& block) {
    std::unique_ptr<RWMol> mol;
    for (bool sanitize : {true, false}) {
        try {
            mol.reset(RDKit::MolBlockToMol(block, sanitize, false));
        } catch (...) {
        }
        if (mol) break;
    }
    if (!mol || !mol->getNumAtoms()) return std::nullopt;
    RDKit::Chirality::reapplyMolBlockWedging(*mol);
    return fromRDKit(*mol);
}

std::string toMolBlock(const Document& doc) {
    auto mol = toRDKit(doc);
    perceive(*mol);
    RDKit::Chirality::reapplyMolBlockWedging(*mol);  // keep the user's wedges
    return RDKit::MolToMolBlock(*mol, true, -1, false);
}

std::string toSmiles(const Document& doc) {
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return "";
    return RDKit::MolToSmiles(*mol);
}

Document clean2D(const Document& doc) {
    if (doc.atoms.empty()) return doc;
    auto mol = toRDKit(doc);
    perceive(*mol);
    RDDepict::preferCoordGen = true;
    RDDepict::compute2DCoords(*mol);
    RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
    Document out = fromRDKit(*mol);

    auto centroid = [](const Document& d) {
        QPointF c;
        for (const auto& a : d.atoms) c += a.pos;
        return c / double(d.atoms.size());
    };
    QPointF shift = centroid(doc) - centroid(out);
    for (auto& a : out.atoms) a.pos += shift;
    return out;
}

std::vector<AtomInfo> atomInfo(const Document& doc) {
    auto mol = toRDKit(doc);
    std::vector<AtomInfo> info(doc.atoms.size());
    for (auto* a : mol->atoms()) {
        auto& i = info[a->getIdx()];
        try {
            a->updatePropertyCache(true);
        } catch (...) {
            i.valenceError = true;
            continue;
        }
        i.hydrogens = int(a->getNumImplicitHs());
    }
    return info;
}

}  // namespace chem
