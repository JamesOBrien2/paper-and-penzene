#include "Chem.h"

#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/CDXMLParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/inchi.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QXmlStreamReader>
#include <QStringList>
#include <cmath>
#include <memory>

namespace chem {
using RDKit::RWMol;

// MOL files and RDKit use 1.5 Å bonds with y up; scenes use points with y down.
constexpr double kScale = kBondLength / 1.5;

// SMILES whose first atom is the attachment point.
static const QHash<QString, QString>& groups() {
    static const QHash<QString, QString> g{
        {"Me", "C"},        {"Et", "CC"},          {"iPr", "C(C)C"},   {"tBu", "C(C)(C)C"},
        {"Ph", "c1ccccc1"}, {"OMe", "OC"},         {"NO2", "[N+](=O)[O-]"}, {"CF3", "C(F)(F)F"},
        {"CN", "C#N"},      {"CO2Me", "C(=O)OC"},  {"CO2Et", "C(=O)OCC"}, {"CO2H", "C(=O)O"},
        {"CHO", "C=O"},     {"Ac", "C(C)=O"},      {"OAc", "OC(C)=O"},  {"N3", "N=[N+]=[N-]"},
        {"Boc", "C(=O)OC(C)(C)C"}, {"Cbz", "C(=O)OCc1ccccc1"}, {"Fmoc", "C(=O)OCC1c2ccccc2-c2ccccc21"},
        {"Bn", "Cc1ccccc1"}, {"Bz", "C(=O)c1ccccc1"}, {"MgBr", "[Mg]Br"}, {"SO2Me", "S(=O)(=O)C"},
        {"Ts", "S(=O)(=O)c1ccc(C)cc1"}, {"Ms", "S(=O)(=O)C"}, {"Tf", "S(=O)(=O)C(F)(F)F"},
        {"OTf", "OS(=O)(=O)C(F)(F)F"}, {"OTs", "OS(=O)(=O)c1ccc(C)cc1"}, {"TMS", "[Si](C)(C)C"},
        {"TBS", "[Si](C)(C)C(C)(C)C"}, {"OTBS", "O[Si](C)(C)C(C)(C)C"}, {"PMB", "Cc1ccc(OC)cc1"},
        {"Bpin", "B1OC(C)(C)C(C)(C)O1"}, {"nBu", "CCCC"}, {"Pr", "CCC"}, {"Cy", "C1CCCCC1"},
    };
    return g;
}

QStringList abbreviations() { return groups().keys(); }

std::optional<Atom> abbreviationHead(const QString& label) {
    auto it = groups().find(label);
    if (it == groups().end()) return std::nullopt;
    auto frag = fromSmiles(it->toStdString());
    return frag ? std::optional(frag->atoms[0]) : std::nullopt;
}

bool attach(Document& doc, int at, const std::string& what) {
    auto frag = fromSmiles(groups().value(QString::fromStdString(what), QString::fromStdString(what)).toStdString());
    if (!frag || frag->atoms.empty()) return false;
    doc.atoms[at].z = frag->atoms[0].z;
    doc.atoms[at].charge = frag->atoms[0].charge;
    doc.atoms[at].label.clear();
    if (frag->atoms.size() == 1) return true;

    // Rotate so the fragment's bulk points away from `at`'s bonds.
    QPointF origin = frag->atoms[0].pos, bulk;
    for (size_t i = 1; i < frag->atoms.size(); ++i) bulk += frag->atoms[i].pos;
    bulk = bulk / double(frag->atoms.size() - 1) - origin;
    QPointF want = doc.neighbors(at).empty() ? QPointF(1, 0) : doc.awayDirection(at);
    double ang = std::atan2(want.y(), want.x()) - std::atan2(bulk.y(), bulk.x());
    const int base = int(doc.atoms.size()) - 1;  // frag atom i -> base + i
    for (size_t i = 1; i < frag->atoms.size(); ++i) {
        QPointF r = frag->atoms[i].pos - origin;
        Atom a = frag->atoms[i];
        a.pos = doc.atoms[at].pos + QPointF(r.x() * std::cos(ang) - r.y() * std::sin(ang),
                                            r.x() * std::sin(ang) + r.y() * std::cos(ang));
        doc.atoms.push_back(a);
    }
    for (auto b : frag->bonds) {
        b.a = b.a == 0 ? at : base + b.a;
        b.b = b.b == 0 ? at : base + b.b;
        doc.bonds.push_back(b);
    }
    return true;
}

Document expanded(const Document& doc) {
    Document out = doc;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if (!doc.atoms[i].label.isEmpty() && !attach(out, i, doc.atoms[i].label.toStdString()))
            out.atoms[i].label.clear();  // unknown label (e.g. from a newer file): keep the atom as is
    return out;
}

// Abbreviations are expanded first; atom i of `doc` is atom i of the mol.
static std::unique_ptr<RWMol> toRDKit(const Document& in) {
    const Document doc = expanded(in);
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
        if (b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash || b.stereo == BondStereo::Wavy) {
            auto* bond = mol->getBondBetweenAtoms(b.a, b.b);
            bond->setBondDir(b.stereo == BondStereo::Wedge  ? RDKit::Bond::BEGINWEDGE
                             : b.stereo == BondStereo::Hash ? RDKit::Bond::BEGINDASH
                                                            : RDKit::Bond::UNKNOWN);
            bond->setProp(RDKit::common_properties::_MolFileBondStereo,
                          b.stereo == BondStereo::Wedge ? 1u : b.stereo == BondStereo::Hash ? 6u : 4u);
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

// ponytail: RDKit depictor + ring templates. CoordGen looks nicer for macrocycles,
// but its preferCoordGen switch is a global the Windows DLL doesn't export.
static void layout(RWMol& mol) {
    RDDepict::Compute2DCoordParameters params;
    params.canonOrient = true;
    params.useRingTemplates = true;
    RDDepict::compute2DCoords(mol, params);
}

// scale 0: normalise whatever bond length the source used to ours.
static Document fromRDKit(RWMol& mol, double scale = 0) {
    try {
        RDKit::MolOps::Kekulize(mol, true);  // draw explicit double bonds
    } catch (...) {
    }
    const auto& conf = mol.getConformer();
    if (!scale) {
        double sum = 0;
        for (const auto* b : mol.bonds())
            sum += (conf.getAtomPos(b->getBeginAtomIdx()) - conf.getAtomPos(b->getEndAtomIdx())).length();
        scale = mol.getNumBonds() && sum > 1e-6 ? kBondLength * mol.getNumBonds() / sum : kScale;
    }

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
                     : b->getBondDir() == RDKit::Bond::UNKNOWN   ? BondStereo::Wavy
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
    layout(*mol);
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

// Arrows and text from the CDXML itself; RDKit only reads the molecules.
// ponytail: plain lines, brackets, shapes and binary .cdx graphics are skipped.
static void chemDrawGraphics(const QByteArray& xml, Document& doc) {
    QXmlStreamReader r(xml);
    double scale = kBondLength / 30;  // CDXML's default BondLength
    QStringList stack;
    auto point = [&](QStringView s) {
        auto v = s.split(' ');
        return v.size() >= 2 ? QPointF(v[0].toDouble(), v[1].toDouble()) * scale : QPointF();
    };
    Text* text = nullptr;
    while (!r.atEnd()) {
        auto tok = r.readNext();
        if (tok == QXmlStreamReader::EndElement) {
            if (r.name() == u"t") text = nullptr;
            stack.removeLast();
            continue;
        }
        if (tok == QXmlStreamReader::Characters && text && stack.last() == "s") {
            text->text += r.text();
            continue;
        }
        if (tok != QXmlStreamReader::StartElement) continue;
        stack.append(r.name().toString());
        const auto at = r.attributes();
        if (r.name() == u"CDXML" && at.hasAttribute("BondLength")) {
            scale = kBondLength / std::max(1.0, at.value("BondLength").toDouble());
        } else if (r.name() == u"t" && !stack.contains("n") && !stack.contains("fragment")) {
            // p is the first baseline; the bounding box gives the left edge whatever the justification.
            QPointF p = point(at.value("p"));
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() == 4) p.setX(std::min(box[0].toDouble(), box[2].toDouble()) * scale);
            doc.texts.push_back({p, {}});
            text = &doc.texts.back();
        } else if (r.name() == u"arrow") {
            const auto head = at.value("ArrowheadHead"), tail = at.value("ArrowheadTail");
            if (head.isEmpty() && tail.isEmpty()) continue;  // a plain line
            Arrow a{point(at.value("Tail3D")), point(at.value("Head3D"))};
            if (!head.isEmpty() && !tail.isEmpty())
                a.kind = at.hasAttribute("ArrowShaftSpacing") ? ArrowKind::Equilibrium : ArrowKind::Resonance;
            else if (at.value("ArrowheadType") == u"Hollow")
                a.kind = ArrowKind::Retro;
            else if (head.startsWith(u"Half") || tail.startsWith(u"Half"))
                a.kind = ArrowKind::Fishhook;
            if (head.isEmpty()) std::swap(a.from, a.to);
            if (double deg = std::abs(at.value("AngularSize").toDouble()); deg > 1 && a.kind != ArrowKind::Equilibrium) {
                // Circular arc about Center3D: bend is how far its midpoint sits off the chord.
                QPointF c = point(at.value("Center3D")), mid = (a.from + a.to) / 2, d = a.to - a.from;
                double radius = std::hypot(point(at.value("MajorAxisEnd3D")).x() - c.x(),
                                           point(at.value("MajorAxisEnd3D")).y() - c.y());
                QPointF out = mid - c;
                double l = std::hypot(out.x(), out.y()), dl = std::hypot(d.x(), d.y());
                if (l > 1e-6 && dl > 1e-6) {
                    QPointF arcMid = c + out / l * radius * (deg <= 180 ? 1 : -1);
                    QPointF n(-d.y() / dl, d.x() / dl);
                    a.bend = -QPointF::dotProduct(arcMid - mid, n);
                }
            }
            doc.arrows.push_back(a);
        }
    }
    std::erase_if(doc.texts, [](const Text& t) { return t.text.trimmed().isEmpty(); });
    for (auto& t : doc.texts) t.text = t.text.trimmed().replace('\r', '\n');
}

std::optional<Document> fromChemDraw(const QByteArray& data) {
    std::vector<std::unique_ptr<RWMol>> mols;
    try {
        // Binary .cdx needs RDKit's ChemDraw library, which not every build has.
        if (!data.trimmed().startsWith('<') && !RDKit::v2::CDXMLParser::hasChemDrawCDXSupport()) return std::nullopt;
        mols = RDKit::v2::CDXMLParser::MolsFromCDXML(data.toStdString());
    } catch (...) {
        return std::nullopt;
    }
    Document doc;
    for (auto& mol : mols) {
        if (!mol->getNumConformers()) continue;
        if (!mol->getNumAtoms()) continue;
        RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
        doc.append(fromRDKit(*mol, kScale));  // RDKit scales CDXML to 1.5 Å bonds, like MOL
    }
    // Nicknames come back expanded plus stray unbonded atoms (a dummy, a lone
    // carbon). ponytail: this also drops a deliberately drawn methane.
    std::vector<int> strays;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if ((doc.atoms[i].z == 0 || (doc.atoms[i].z == 6 && !doc.atoms[i].charge)) && doc.neighbors(i).empty())
            strays.push_back(i);
    doc.removeAtoms(strays);
    if (data.trimmed().startsWith('<')) chemDrawGraphics(data, doc);
    if (doc.empty()) return std::nullopt;
    return doc;
}

std::optional<Document> readFile(const QString& path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return std::nullopt;
    const QByteArray data = f.readAll();
    const QString ext = QFileInfo(path).suffix().toLower();
    if (ext == "penz") return Document::fromJson(data);
    if (ext == "cdxml" || ext == "cdx") return fromChemDraw(data);
    return fromMolBlock(data.toStdString());
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

std::optional<Properties> properties(const Document& doc) {
    if (doc.atoms.empty()) return std::nullopt;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return std::nullopt;
    return Properties{RDKit::Descriptors::calcMolFormula(*mol), RDKit::Descriptors::calcAMW(*mol),
                      RDKit::Descriptors::calcExactMW(*mol)};
}

std::string toInchi(const Document& doc) {
    if (doc.atoms.empty()) return "";
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return "";
    RDKit::ExtraInchiReturnValues rv;
    return RDKit::MolToInchi(*mol, rv);
}

std::string toInchiKey(const Document& doc) {
    std::string inchi = toInchi(doc);
    return inchi.empty() ? "" : RDKit::InchiToInchiKey(inchi);
}

// One connected fragment, laid out around its old centroid.
static Document cleanFragment(const Document& doc) {
    auto mol = toRDKit(doc);
    perceive(*mol);
    layout(*mol);
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

// Each fragment is cleaned in place, so a reaction scheme keeps its layout;
// arrows and text pass through untouched.
Document clean2D(const Document& doc, const std::vector<int>& only) {
    const int n = int(doc.atoms.size());
    std::vector<int> comp(n, -1);
    int count = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        std::vector<int> stack{s};
        comp[s] = count;
        while (!stack.empty()) {
            int i = stack.back();
            stack.pop_back();
            for (int nb : doc.neighbors(i))
                if (comp[nb] < 0) comp[nb] = count, stack.push_back(nb);
        }
        ++count;
    }
    std::vector<bool> wanted(count, only.empty());
    for (int i : only)
        if (i >= 0 && i < n) wanted[comp[i]] = true;
    Document out = doc;
    out.bonds.clear();
    for (const auto& b : doc.bonds)  // molecules left alone keep their bonds as drawn
        if (!wanted[comp[b.a]]) out.bonds.push_back(b);
    for (int c = 0; c < count; ++c) {
        if (!wanted[c]) continue;
        std::vector<int> ids, drop;
        for (int i = 0; i < n; ++i) (comp[i] == c ? ids : drop).push_back(i);
        Document frag = doc;
        frag.arrows.clear(), frag.texts.clear();
        frag.removeAtoms(drop);  // keeps order: frag atom k is doc atom ids[k]
        Document clean = cleanFragment(frag);
        for (size_t k = 0; k < ids.size(); ++k) out.atoms[ids[k]].pos = clean.atoms[k].pos;
        const int m = int(ids.size());  // atoms past m are expanded abbreviations: dropped again
        for (auto b : clean.bonds)
            if (b.a < m && b.b < m) b.a = ids[b.a], b.b = ids[b.b], out.bonds.push_back(b);
    }
    return out;
}

std::vector<std::vector<int>> rings(const Document& doc) {
    auto mol = toRDKit(doc);
    std::vector<std::vector<int>> out;
    try {
        RDKit::MolOps::findSSSR(*mol);
    } catch (...) {
        return out;
    }
    const int n = int(doc.atoms.size());
    for (const auto& r : mol->getRingInfo()->atomRings())
        if (std::all_of(r.begin(), r.end(), [n](int i) { return i < n; })) out.push_back(r);
    return out;
}

std::vector<AtomInfo> atomInfo(const Document& doc) {
    auto mol = toRDKit(doc);
    std::vector<AtomInfo> info(doc.atoms.size());
    for (auto* a : mol->atoms()) {
        if (a->getIdx() >= info.size()) break;  // expanded abbreviation atoms
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

// Via Atom, not PeriodicTable: its inline methods reference a logger global
// that the Windows RDKit DLL doesn't export.
std::string symbol(int z) {
    return RDKit::Atom(z).getSymbol();
}

int atomicNumber(const std::string& sym) {
    try {
        return int(RDKit::Atom(sym).getAtomicNum());
    } catch (...) {
        return 0;
    }
}

}  // namespace chem
