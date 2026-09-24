#include "Chem.h"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/FileParsers/CDXMLParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <map>
#include <GraphMol/inchi.h>
#if __has_include(<GraphMol/chemdraw.h>)
#include <GraphMol/chemdraw.h>
#define PENZENE_CDX_WRITER 1
#endif
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QObject>
#include <QPolygonF>
#include <QRegularExpression>
#include <QLineF>
#include <QSet>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QStringList>
#include <cmath>
#include <numbers>
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

// Abbreviations are expanded first (unless `expand` is false, for layout);
// atom i of `doc` is atom i of the mol.
static std::unique_ptr<RWMol> toRDKit(const Document& in, bool expand = true) {
    const Document doc = expand ? expanded(in) : in;
    auto mol = std::make_unique<RWMol>();
    auto* conf = new RDKit::Conformer(doc.atoms.size());
    for (size_t i = 0; i < doc.atoms.size(); ++i) {
        const auto& a = doc.atoms[i];
        auto* atom = new RDKit::Atom(a.z);
        atom->setFormalCharge(a.charge);
        atom->setNumRadicalElectrons(a.radicals);
        // Not setAtomMapNum: it logs through rdErrorLog, which the Windows DLL doesn't export.
        if (a.map > 0) atom->setProp(RDKit::common_properties::molAtomMapNumber, a.map);
        // Generic atoms (expansion drops their labels, so they come from `in`): R1…Rn as
        // MDL R-groups (R# with RGP, [n*] in SMILES), any other text (X, Ar) as an MDL atom alias.
        // ponytail: V3000 has no alias block, so there X/Ar become plain * atoms.
        if (a.z == 0 && i < in.atoms.size() && !in.atoms[i].label.isEmpty()) {
            static const QRegularExpression rgroup("^R(\\d+)$");
            const QString label = in.atoms[i].label;
            if (const auto m = rgroup.match(label); m.hasMatch()) {
                const unsigned n = m.captured(1).toUInt();
                atom->setProp(RDKit::common_properties::_MolFileRLabel, n);
                atom->setIsotope(n);
            } else {
                atom->setProp(RDKit::common_properties::molFileAlias, label.toStdString());
            }
        }
        mol->addAtom(atom, true, true);
        conf->setAtomPos(i, {a.pos.x() / kScale, -a.pos.y() / kScale, 0});
    }
    for (const auto& b : doc.bonds) {
        const int order = chemicalOrder(b);
        if (order < 1) continue;  // drawn only
        auto type = order == 2 ? RDKit::Bond::DOUBLE
                    : order == 3 ? RDKit::Bond::TRIPLE
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
static bool perceive(RWMol& mol, bool aromatic = true) {
    try {
        unsigned failed = 0;
        RDKit::MolOps::sanitizeMol(mol, failed,
                                   aromatic ? RDKit::MolOps::SANITIZE_ALL
                                            : RDKit::MolOps::SANITIZE_ALL ^ RDKit::MolOps::SANITIZE_SETAROMATICITY);
    } catch (...) {
        mol.updatePropertyCache(false);
        RDKit::MolOps::fastFindRings(mol);
        return false;
    }
    RDKit::MolOps::assignChiralTypesFromBondDirs(mol);  // stereocentres from wedges
    if (mol.getNumConformers()) RDKit::MolOps::detectBondStereochemistry(mol);  // E/Z from the drawing
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
        QString label;  // a generic atom's: its alias (X, Ar), or R-group number
        unsigned r = 0;
        std::string alias;
        if (a->getAtomicNum() == 0 && a->getPropIfPresent(RDKit::common_properties::molFileAlias, alias) && !alias.empty())
            label = QString::fromStdString(alias);
        else if (a->getAtomicNum() == 0 && a->getPropIfPresent(RDKit::common_properties::_MolFileRLabel, r) && r)
            label = QString("R%1").arg(r);
        doc.atoms.push_back({QPointF(p.x * scale, -p.y * scale), int(a->getAtomicNum()),
                             a->getFormalCharge(), label, {}, int(a->getAtomMapNum()), 0,
                             int(std::min(2u, a->getNumRadicalElectrons()))});
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

std::optional<Document> fromInchi(const std::string& inchi) {
    RDKit::ExtraInchiReturnValues rv;
    std::unique_ptr<RWMol> mol;
    try {
        mol.reset(RDKit::InchiToMol(inchi, rv));
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

// A ChemDraw label node (nickname, generic R, or free text drawn as an atom).
struct LabelNode {
    int id = 0;
    QPointF pos, textPos;
    QString text;
    double textScale = 1;
    std::vector<QPointF> inner;  // positions of the atoms inside its fragment
};

// Arrows, free text and label nodes from the CDXML itself; RDKit only reads the
// molecules. Returns the label nodes and how many bonds each node id has.
// ponytail: plain lines, brackets, shapes and binary .cdx graphics are skipped.
static std::vector<LabelNode> chemDrawGraphics(const QByteArray& xml, Document& doc, QHash<int, int>& bondCount) {
    QXmlStreamReader r(xml);
    double scale = kBondLength / 30;  // CDXML's default BondLength
    struct Open { QString tag; int label = -1; };  // label: index into `labels` for label <n>s
    std::vector<Open> stack;
    std::vector<LabelNode> labels;
    auto point = [&](QStringView s) {
        auto v = s.split(' ');
        return v.size() >= 2 ? QPointF(v[0].toDouble(), v[1].toDouble()) * scale : QPointF();
    };
    // Where characters go: a free Text, or a label node's own text (not its inner fragment's).
    Text* text = nullptr;
    int labelText = -1;
    auto inside = [&](const char* tag) {
        return std::any_of(stack.begin(), stack.end(), [&](const Open& o) { return o.tag == tag; });
    };
    static const QStringList labelTypes{"Fragment",    "Nickname",  "GenericNickname", "Unspecified", "Anonymous",
                                        "AnonymousAlternativeGroup", "NamedAlternativeGroup", "Variable"};
    while (!r.atEnd()) {
        auto tok = r.readNext();
        if (tok == QXmlStreamReader::EndElement) {
            if (r.name() == u"t") text = nullptr, labelText = -1;
            stack.pop_back();
            continue;
        }
        if (tok == QXmlStreamReader::Characters && !stack.empty() && stack.back().tag == "s") {
            if (text) text->text += r.text();
            if (labelText >= 0) labels[labelText].text += r.text();
            continue;
        }
        if (tok != QXmlStreamReader::StartElement) continue;
        const QString tag = r.name().toString();
        const auto at = r.attributes();
        const int parentLabel = stack.empty() ? -1 : stack.back().label;
        stack.push_back({tag});
        if (tag == "CDXML" && at.hasAttribute("BondLength")) {
            scale = kBondLength / std::max(1.0, at.value("BondLength").toDouble());
        } else if (tag == "b") {
            ++bondCount[at.value("B").toInt()], ++bondCount[at.value("E").toInt()];
        } else if (tag == "n") {
            const int id = at.value("id").toInt();
            for (const Open& o : stack)  // an atom inside a label's fragment
                if (o.label >= 0) labels[o.label].inner.push_back(point(at.value("p")));
            if (labelTypes.contains(at.value("NodeType").toString())) {
                labels.push_back({id, point(at.value("p"))});
                stack.back().label = int(labels.size()) - 1;
            }
        } else if (tag == "t" && parentLabel >= 0) {
            labels[parentLabel].textPos = point(at.value("p"));
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() == 4) labels[parentLabel].textPos.setX(std::min(box[0].toDouble(), box[2].toDouble()) * scale);
            labelText = parentLabel;
        } else if (tag == "s" && (text || labelText >= 0) && at.hasAttribute("size")) {
            const double rel = at.value("size").toDouble() * scale / 10;  // 10 pt: the default (ACS) label size
            (text ? text->scale : labels[labelText].textScale) = rel;
        } else if (tag == "t" && !inside("n") && !inside("fragment")) {
            // p is the first baseline; the bounding box gives the left edge whatever the justification.
            QPointF p = point(at.value("p"));
            auto box = at.value("BoundingBox").split(' ');
            if (box.size() == 4) p.setX(std::min(box[0].toDouble(), box[2].toDouble()) * scale);
            doc.texts.push_back({p, {}});
            text = &doc.texts.back();
        } else if (tag == "arrow") {
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
    for (auto& l : labels) l.text = l.text.trimmed().replace('\r', '\n');
    return labels;
}

// RDKit expands label nodes itself, at their fragment's own coordinates (often
// far off the page), and leaves unbonded ones as stray atoms. Put them back the
// way ChemDraw shows them (#86):
//  - an unbonded label is text (reagents such as "LiBr, acetone");
//  - a bonded label is one labelled atom, as ChemDraw draws it: a known
//    abbreviation keeps its chemistry, an unknown one (SCoA) becomes an
//    unknown group (*) rather than pretending to be its attachment atom;
//  - a label RDKit kept as one atom (R, X) keeps its text as the atom's label.
// RDKit tags only some atoms with their node id, but keeps every atom at its
// CDXML position (even inside a nickname's own fragment), so match by position.
static void placeLabels(Document& doc, const std::vector<int>& nodeOf, const std::vector<LabelNode>& labels,
                        const QHash<int, int>& bondCount) {
    std::vector<int> drop;
    auto at = [](QPointF a, QPointF b) { return std::abs(a.x() - b.x()) < 0.6 && std::abs(a.y() - b.y()) < 0.6; };
    for (const LabelNode& l : labels) {
        std::vector<int> atoms;
        for (int i = 0; i < int(doc.atoms.size()); ++i) {
            const QPointF p = doc.atoms[i].pos;
            if (nodeOf[i] == l.id || at(p, l.pos) ||
                std::any_of(l.inner.begin(), l.inner.end(), [&](QPointF q) { return at(p, q); }))
                atoms.push_back(i);
        }
        if (bondCount.value(l.id) == 0) {
            if (!l.text.isEmpty()) doc.texts.push_back({l.textPos, l.text, l.textScale});
            drop.insert(drop.end(), atoms.begin(), atoms.end());
            continue;
        }
        if (atoms.empty()) continue;
        auto self = std::find_if(atoms.begin(), atoms.end(), [&](int i) {
            return nodeOf[i] == l.id || (at(doc.atoms[i].pos, l.pos) && l.inner.empty());
        });
        if (self != atoms.end()) {  // kept as a single atom
            doc.atoms[*self].pos = l.pos;
            if (!l.text.isEmpty()) doc.atoms[*self].label = l.text;
            continue;
        }
        const QSet<int> group(atoms.begin(), atoms.end());
        auto attach = std::find_if(atoms.begin(), atoms.end(), [&](int i) {
            auto n = doc.neighbors(i);
            return std::any_of(n.begin(), n.end(), [&](int j) { return !group.contains(j); });
        });
        if (attach == atoms.end()) continue;
        const auto head = abbreviationHead(l.text);
        Atom& a = doc.atoms[*attach];
        a.z = head ? head->z : 0, a.charge = head ? head->charge : 0, a.label = l.text, a.pos = l.pos;
        for (int i : atoms)
            if (i != *attach) drop.push_back(i);
    }
    doc.removeAtoms(drop);
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
    std::vector<int> nodeOf;  // CDX node id of each atom, to match label nodes
    for (auto& mol : mols) {
        if (!mol->getNumConformers()) continue;
        if (!mol->getNumAtoms()) continue;
        RDKit::Chirality::wedgeMolBonds(*mol, &mol->getConformer());
        for (const auto* a : mol->atoms()) {
            unsigned id = 0;
            a->getPropIfPresent("CDX_NODE_ID", id);
            nodeOf.push_back(int(id));
        }
        doc.append(fromRDKit(*mol, kScale));  // RDKit scales CDXML to 1.5 Å bonds, like MOL
    }
    if (data.trimmed().startsWith('<')) {
        Document graphics;
        QHash<int, int> bondCount;
        placeLabels(doc, nodeOf, chemDrawGraphics(data, graphics, bondCount), bondCount);
        doc.arrows = graphics.arrows;
        doc.texts.insert(doc.texts.end(), graphics.texts.begin(), graphics.texts.end());
    }
    // Any dummy atom left unbonded and unlabelled is an artefact of the expansion.
    std::vector<int> strays;
    for (int i = 0; i < int(doc.atoms.size()); ++i)
        if (doc.atoms[i].z == 0 && doc.atoms[i].label.isEmpty() && doc.neighbors(i).empty()) strays.push_back(i);
    doc.removeAtoms(strays);
    if (doc.empty()) return std::nullopt;
    return doc;
}

std::vector<Record> readRecords(const QString& path) {
    const QFileInfo info(path);
    const QString base = info.completeBaseName(), ext = info.suffix().toLower();
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return {};
    std::vector<Record> out;
    const QString text = QString::fromUtf8(f.readAll()).remove('\r');
    const QStringList parts = ext == "sdf" ? text.split("$$$$") : text.split('\n');
    for (const QString& raw : parts) {
        const QString part = ext == "sdf" ? raw : raw.trimmed();
        if (part.trimmed().isEmpty() || part.startsWith('#')) continue;
        const QString n = QString("%1-%2").arg(base).arg(out.size() + 1);
        if (ext == "sdf") {
            const QString block = part.startsWith('\n') ? part.mid(1) : part;  // after "$$$$\n"
            const QString title = block.section('\n', 0, 0).trimmed();  // the molfile's name line
            out.push_back({title.isEmpty() ? n : title, fromMolBlock(block.toStdString())});
        } else {
            const QStringList cols = part.split(QRegularExpression("\\s+"));
            const std::string first = cols[0].toStdString();
            out.push_back({cols.size() > 1 ? cols[1] : n, ext == "inchi" ? fromInchi(first) : fromSmiles(first)});
        }
    }
    return out;
}

static QRectF atomBox(const Document& d) {
    QPolygonF pts;
    for (const auto& a : d.atoms) pts << a.pos;
    return pts.boundingRect();
}

// Records side by side, row by row, each centred in a cell as big as the largest.
static std::optional<Document> grid(const std::vector<Record>& records) {
    std::vector<Document> docs;
    for (const auto& r : records)
        if (r.doc && !r.doc->empty()) docs.push_back(*r.doc);
    if (docs.size() <= 1) return docs.empty() ? std::nullopt : std::optional(docs[0]);
    auto box = atomBox;
    QSizeF cell;
    for (const auto& d : docs) cell = cell.expandedTo(box(d).size());
    cell += QSizeF(2 * kBondLength, 2 * kBondLength);
    const int cols = int(std::ceil(std::sqrt(double(docs.size()))));
    Document out;
    for (size_t k = 0; k < docs.size(); ++k) {
        const QPointF centre((k % cols) * cell.width(), (k / cols) * cell.height());
        out.append(docs[k], centre - box(docs[k]).center());
    }
    return out;
}

// ---------------------------------------------------------------- reactions

static std::vector<Document> molecules(const Document& doc) {
    std::vector<Document> out;
    std::vector<int> seen(doc.atoms.size(), 0);
    for (int start = 0; start < int(doc.atoms.size()); ++start) {
        if (seen[start]) continue;
        std::vector<int> stack{start}, keep;
        seen[start] = 1;
        while (!stack.empty()) {
            const int i = stack.back();
            stack.pop_back();
            keep.push_back(i);
            for (int nb : doc.neighbors(i))
                if (!seen[nb]) seen[nb] = 1, stack.push_back(nb);
        }
        Document m;
        m.atoms = doc.atoms;
        m.bonds = doc.bonds;
        std::vector<int> drop;
        std::sort(keep.begin(), keep.end());
        for (int i = 0; i < int(doc.atoms.size()); ++i)
            if (!std::binary_search(keep.begin(), keep.end(), i)) drop.push_back(i);
        m.removeAtoms(drop);
        out.push_back(std::move(m));
    }
    return out;
}

// ponytail: the first straight reaction or equilibrium arrow only; multi-step
// schemes would need a split per arrow.
std::optional<Reaction> reactionOf(const Document& doc) {
    auto arrow = std::find_if(doc.arrows.begin(), doc.arrows.end(), [](const Arrow& a) {
        return a.bend == 0 && (a.kind == ArrowKind::Reaction || a.kind == ArrowKind::Equilibrium);
    });
    if (arrow == doc.arrows.end()) return std::nullopt;
    const QPointF d = arrow->to - arrow->from;
    const double len2 = QPointF::dotProduct(d, d);
    if (len2 <= 0) return std::nullopt;
    Reaction r;
    for (auto& m : molecules(doc)) {
        const double t = QPointF::dotProduct(atomBox(m).center() - arrow->from, d) / len2;
        (t < 0 ? r.reactants : t > 1 ? r.products : r.agents).push_back(std::move(m));
    }
    return r;
}

std::string toReactionSmiles(const Reaction& r) {
    auto side = [](const std::vector<Document>& ms) {
        std::string s;
        for (const auto& m : ms) s += (s.empty() ? "" : ".") + toSmiles(m);
        return s;
    };
    return side(r.reactants) + ">" + side(r.agents) + ">" + side(r.products);
}

std::string toRxn(const Reaction& r) {
    std::string out = "$RXN\n\n  Penzene\n\n" + QString("%1%2").arg(r.reactants.size(), 3).arg(r.products.size(), 3).toStdString() + "\n";
    for (const auto* side : {&r.reactants, &r.products})
        for (const auto& m : *side) out += "$MOL\n" + toMolBlock(m);
    return out;
}

// Reactants + … → (agents above the arrow) → products, left to right.
Document layoutReaction(const Reaction& r) {
    Document out;
    const double gap = kBondLength;
    double x = 0;
    auto place = [&](const std::vector<Document>& ms) {
        for (size_t k = 0; k < ms.size(); ++k) {
            if (k) {  // a "+" between molecules
                out.texts.push_back({{x + gap, kBondLength * 0.25}, "+"});
                x += 3 * gap;
            }
            const QRectF b = atomBox(ms[k]);
            out.append(ms[k], QPointF(x - b.left(), -b.center().y()));
            x += b.width();
        }
    };
    place(r.reactants);
    double above = 0, width = 0;
    for (const auto& m : r.agents) width += atomBox(m).width() + gap, above = std::max(above, atomBox(m).height());
    const double len = std::max(3 * kBondLength, width + gap);
    out.arrows.push_back({{x + gap, 0}, {x + gap + len, 0}});
    double ax = x + gap + (len - width + gap) / 2;
    for (const auto& m : r.agents) {
        const QRectF b = atomBox(m);
        out.append(m, QPointF(ax - b.left(), -gap - above / 2 - b.center().y()));
        ax += b.width() + gap;
    }
    x += 2 * gap + len;
    place(r.products);
    return out;
}

std::optional<Document> fromReactionSmiles(const std::string& smiles) {
    const QStringList sides = QString::fromStdString(smiles).trimmed().split('>');
    if (sides.size() != 3) return std::nullopt;
    Reaction r;
    std::vector<Document>* into[] = {&r.reactants, &r.agents, &r.products};
    for (int k = 0; k < 3; ++k)
        for (const QString& part : sides[k].split('.', Qt::SkipEmptyParts)) {
            auto m = fromSmiles(part.toStdString());
            if (!m) return std::nullopt;
            into[k]->push_back(std::move(*m));
        }
    if (r.reactants.empty() && r.products.empty()) return std::nullopt;
    return layoutReaction(r);
}

std::optional<Document> fromRxn(const std::string& text) {
    const QString s = QString::fromStdString(text).remove('\r');
    if (!s.startsWith("$RXN")) return std::nullopt;
    const QString counts = s.section('\n', 4, 4);
    const int nr = counts.mid(0, 3).trimmed().toInt(), np = counts.mid(3, 3).trimmed().toInt();
    QStringList blocks = s.split("$MOL\n");
    blocks.removeFirst();
    if (nr + np == 0 || int(blocks.size()) < nr + np) return std::nullopt;
    Reaction r;
    for (int k = 0; k < nr + np; ++k) {
        auto m = fromMolBlock(blocks[k].toStdString());
        if (!m) return std::nullopt;
        (k < nr ? r.reactants : r.products).push_back(std::move(*m));
    }
    return layoutReaction(r);
}

std::optional<Document> readFile(const QString& path) {
    const QString ext = QFileInfo(path).suffix().toLower();
    if (ext == "sdf" || ext == "smi" || ext == "inchi") return grid(readRecords(path));
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return std::nullopt;
    const QByteArray data = f.readAll();
    if (ext == "penz") return Document::fromJson(data);
    if (ext == "rxn") return fromRxn(data.toStdString());
    if (ext == "png" || ext == "svg") return Document::fromEmbedded(data);
    if (ext == "cdxml" || ext == "cdx") return fromChemDraw(data);
    return fromMolBlock(data.toStdString());
}

// CDXML in our own coordinates (BondLength = ours, y down, as ChemDraw), so
// ChemDraw and chemDrawGraphics read it back unscaled.
// ponytail: abbreviations are written expanded, free-text labels as generic
// nicknames; ChemDraw's own Fragment/Nickname nodes would keep "OMe" as a label.
QByteArray toCdxml(const Document& in) {
    const Document doc = expanded(in);
    QByteArray out;
    QXmlStreamWriter w(&out);
    w.setAutoFormatting(true);
    w.writeStartDocument();
    w.writeDTD(R"(<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">)");
    int id = 1;
    auto pt = [](QPointF p) { return QString("%1 %2").arg(p.x(), 0, 'f', 2).arg(p.y(), 0, 'f', 2); };
    auto pt3 = [&](QPointF p) { return pt(p) + " 0"; };
    w.writeStartElement("CDXML");
    w.writeAttribute("BondLength", QString::number(kBondLength));
    w.writeAttribute("CreationProgram", "Penzene");
    w.writeStartElement("page");
    w.writeAttribute("id", QString::number(id++));
    if (!doc.atoms.empty()) {
        w.writeStartElement("fragment");
        w.writeAttribute("id", QString::number(id++));
        const int base = id;
        for (size_t i = 0; i < doc.atoms.size(); ++i) {
            const Atom& a = doc.atoms[i];
            // Expansion drops labels it can't draw out; free text (R, X, MgEt) is in the original.
            const QString label = i < in.atoms.size() && a.z == 0 ? in.atoms[i].label : QString();
            w.writeStartElement("n");
            w.writeAttribute("id", QString::number(id++));
            w.writeAttribute("p", pt(a.pos));
            if (!label.isEmpty()) {
                w.writeAttribute("NodeType", "GenericNickname");
                w.writeAttribute("GenericNickname", label);
                w.writeStartElement("t");
                w.writeAttribute("p", pt(a.pos + QPointF(-3, 4)));
                w.writeTextElement("s", label);
                w.writeEndElement();
            } else if (a.z != 6) {
                w.writeAttribute("Element", QString::number(a.z));
            }
            if (a.charge) w.writeAttribute("Charge", QString::number(a.charge));
            w.writeEndElement();
        }
        static const char* display[] = {nullptr, "WedgeBegin", "WedgedHashBegin", "Bold", "Dash", "Wavy", nullptr, "Dash"};
        static const char* side[] = {nullptr, "Left", "Center", "Right"};
        for (const Bond& b : doc.bonds) {
            w.writeStartElement("b");
            w.writeAttribute("id", QString::number(id++));
            w.writeAttribute("B", QString::number(base + b.a));
            w.writeAttribute("E", QString::number(base + b.b));
            // ChemDraw's own hydrogen-bond order, and half orders for partial bonds.
            if (b.stereo == BondStereo::Interaction) w.writeAttribute("Order", "hydrogen");
            else if (b.stereo == BondStereo::Partial) w.writeAttribute("Order", b.order == 2 ? "1.5" : "0.5");
            else if (b.order > 1) w.writeAttribute("Order", QString::number(b.order));
            if (display[int(b.stereo)]) w.writeAttribute("Display", display[int(b.stereo)]);
            if (side[int(b.position)]) w.writeAttribute("DoublePosition", side[int(b.position)]);
            w.writeEndElement();
        }
        w.writeEndElement();
    }
    for (const Text& t : doc.texts) {
        w.writeStartElement("t");
        w.writeAttribute("id", QString::number(id++));
        w.writeAttribute("p", pt(t.pos));
        w.writeStartElement("s");
        w.writeAttribute("size", QString::number(10 * t.scale));  // 10 pt: the ACS label size
        w.writeCharacters(QString(t.text).replace('\n', '\r'));
        w.writeEndElement();
        w.writeEndElement();
    }
    for (const Arrow& a : doc.arrows) {
        w.writeStartElement("arrow");
        w.writeAttribute("id", QString::number(id++));
        w.writeAttribute("Head3D", pt3(a.to));
        w.writeAttribute("Tail3D", pt3(a.from));
        switch (a.kind) {
        case ArrowKind::Reaction: w.writeAttribute("ArrowheadHead", "Full"); break;
        case ArrowKind::Retro:
            w.writeAttribute("ArrowheadHead", "Full");
            w.writeAttribute("ArrowheadType", "Hollow");
            break;
        case ArrowKind::Resonance:
            w.writeAttribute("ArrowheadHead", "Full");
            w.writeAttribute("ArrowheadTail", "Full");
            break;
        case ArrowKind::Equilibrium:
            w.writeAttribute("ArrowheadHead", "HalfLeft");
            w.writeAttribute("ArrowheadTail", "HalfLeft");
            w.writeAttribute("ArrowShaftSpacing", "4");
            break;
        case ArrowKind::Fishhook: w.writeAttribute("ArrowheadHead", "HalfLeft"); break;
        }
        if (a.kind != ArrowKind::Retro) w.writeAttribute("ArrowheadType", "Solid");
        if (std::abs(a.bend) > 1e-6 && a.kind != ArrowKind::Equilibrium) {
            // The circle through both ends and the arc's midpoint (bend off the chord, as read back).
            const QPointF d = a.to - a.from, mid = (a.from + a.to) / 2;
            const double c = std::hypot(d.x(), d.y()) / 2, s = std::abs(a.bend);
            const QPointF n(-d.y() / (2 * c), d.x() / (2 * c)), arcMid = mid - a.bend * n;
            const double r = (c * c + s * s) / (2 * s);
            const QPointF centre = arcMid + (mid - arcMid) / s * r;
            auto angle = [&](QPointF p) { return std::atan2(p.y() - centre.y(), p.x() - centre.x()) * 180 / std::numbers::pi; };
            auto wrap = [](double x) { return std::fmod(std::fmod(x, 360) + 360, 360); };
            double sweep = wrap(angle(a.to) - angle(a.from));
            if (wrap(angle(arcMid) - angle(a.from)) > sweep) sweep -= 360;
            w.writeAttribute("Center3D", pt3(centre));
            w.writeAttribute("MajorAxisEnd3D", pt3(centre + QPointF(r, 0)));
            w.writeAttribute("MinorAxisEnd3D", pt3(centre + QPointF(0, r)));
            w.writeAttribute("AngularSize", QString::number(-sweep, 'f', 2));  // ChemDraw: head = tail turned by -AngularSize
        }
        w.writeEndElement();
    }
    w.writeEndDocument();
    return out;
}

// Binary CDX, through RDKit: molecules only (no arrows or text). Empty where
// RDKit has no ChemDraw writer.
QByteArray toCdx(const Document& doc) {
#ifndef PENZENE_CDX_WRITER
    return {};
#else
    auto mol = toRDKit(doc);
    perceive(*mol);
    try {
        return QByteArray::fromStdString(RDKit::v2::MolToChemDrawBlock(*mol, RDKit::v2::CDXFormat::CDX));
    } catch (...) {
        return {};
    }
#endif
}

std::string toMolBlock(const Document& doc, bool v3000) {
    auto mol = toRDKit(doc);
    perceive(*mol);
    RDKit::Chirality::reapplyMolBlockWedging(*mol);  // keep the user's wedges
    return RDKit::MolToMolBlock(*mol, true, -1, false, v3000);
}

std::string toSmiles(const Document& doc) {
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return "";
    std::string smiles = RDKit::MolToSmiles(*mol);
    // RDKit can call a ring with an odd charged atom ([b-2]1ccccc1) aromatic yet not
    // read that SMILES back; write it in Kekulé form then.
    std::unique_ptr<RWMol> back;
    try {
        back.reset(RDKit::SmilesToMol(smiles));
    } catch (...) {
    }
    if (!back) {  // as drawn, without aromaticity
        mol = toRDKit(doc);
        if (perceive(*mol, false)) smiles = RDKit::MolToSmiles(*mol);
    }
    return smiles;
}

std::optional<Properties> properties(const Document& doc) {
    if (doc.atoms.empty()) return std::nullopt;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return std::nullopt;
    return Properties{RDKit::Descriptors::calcMolFormula(*mol), RDKit::Descriptors::calcAMW(*mol),
                      RDKit::Descriptors::calcExactMW(*mol)};
}

std::optional<Profile> profile(const Document& doc) {
    auto basic = properties(doc);
    if (!basic) return std::nullopt;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return std::nullopt;
    Profile p{*basic};
    double mr = 0;
    RDKit::Descriptors::calcCrippenDescriptors(*mol, p.logP, mr);
    p.tpsa = RDKit::Descriptors::calcTPSA(*mol);
    p.hbd = int(RDKit::Descriptors::calcNumHBD(*mol));
    p.hba = int(RDKit::Descriptors::calcNumHBA(*mol));
    p.rotatable = int(RDKit::Descriptors::calcNumRotatableBonds(*mol));
    p.heavyAtoms = int(RDKit::Descriptors::calcNumHeavyAtoms(*mol));
    // Mass fractions per element, hydrogens included; masses via Atom (PeriodicTable's
    // inline methods don't link on Windows).
    RDKit::RWMol withH(*mol);
    RDKit::MolOps::addHs(withH);
    std::map<std::string, double> mass;
    for (const auto* a : withH.atoms()) mass[a->getSymbol()] += a->getMass();
    double total = 0;
    for (const auto& [s, m] : mass) total += m;
    std::vector<std::string> order;  // Hill: C, H, then alphabetical
    for (const char* s : {"C", "H"})
        if (mass.count(s)) order.push_back(s);
    for (const auto& [s, m] : mass)
        if (s != "C" && s != "H") order.push_back(s);
    for (const auto& s : order) p.elemental.push_back({s, 100 * mass[s] / total});
    p.lipinskiViolations = (p.basic.mw > 500) + (p.logP > 5) + (p.hbd > 5) + (p.hba > 10);
    p.veber = p.rotatable <= 10 && p.tpsa <= 140;
    return p;
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
    // Abbreviations stay single nodes, as in ChemDraw: expanding a ring onto a
    // crowded atom squeezes the depictor's layout (#85).
    auto mol = toRDKit(doc, false);
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
        const int m = int(ids.size());  // defensive: only bonds between this fragment's atoms
        for (auto b : clean.bonds) {
            if (b.a >= m || b.b >= m) continue;
            b.a = ids[b.a], b.b = ids[b.b];
            if (int o = doc.bondBetween(b.a, b.b); o >= 0 && chemicalOrder(doc.bonds[o]) != doc.bonds[o].order)
                continue;  // RDKit's view of a partial bond; the drawn one goes back below
            // RDKit only knows order and wedges: keep what the user chose for display.
            if (int o = doc.bondBetween(b.a, b.b); o >= 0) {
                const Bond& was = doc.bonds[o];
                const bool styled = was.stereo == BondStereo::Bold || was.stereo == BondStereo::Dashed ||
                                    was.stereo == BondStereo::Partial;
                if (styled && b.stereo == BondStereo::None) b.stereo = was.stereo;
                if (b.order == 2 && was.order == 2) b.position = was.position;
            }
            out.bonds.push_back(b);
        }
        // Interactions and partial bonds, which RDKit never saw, as drawn.
        for (const auto& b : doc.bonds)
            if (comp[b.a] == c && chemicalOrder(b) != b.order) out.bonds.push_back(b);
    }
    return out;
}

Document addHydrogens(const Document& doc) {
    Document out = doc;
    const int n = int(doc.atoms.size());
    if (!n) return out;
    auto mol = toRDKit(doc);
    perceive(*mol);
    const int before = int(mol->getNumAtoms());
    try {
        RDKit::MolOps::addHs(*mol, false, true);  // with 2D coordinates
    } catch (...) {
        return out;
    }
    const auto& conf = mol->getConformer();
    for (int i = before; i < int(mol->getNumAtoms()); ++i) {
        const auto* h = mol->getAtomWithIdx(i);
        const int heavy = int((*mol->getAtomNeighbors(h).first));
        if (heavy >= n || !doc.atoms[heavy].label.isEmpty()) continue;  // not onto abbreviations
        const auto& p = conf.getAtomPos(i);
        const int j = out.addAtom({p.x * kScale, -p.y * kScale}, 1);
        out.bonds.push_back({heavy, j});
    }
    return out;
}

Document removeHydrogens(const Document& doc) {
    Document out = doc;
    std::vector<int> drop;
    for (int i = 0; i < int(doc.atoms.size()); ++i) {
        const Atom& a = doc.atoms[i];
        if (a.z != 1 || a.charge || !a.label.isEmpty() || doc.neighbors(i).size() != 1) continue;
        const Bond& b = doc.bonds[doc.bondBetween(i, doc.neighbors(i)[0])];
        if (b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash) continue;  // stereo H stays
        if (doc.atoms[doc.neighbors(i)[0]].z == 1) continue;  // H2 stays
        drop.push_back(i);
    }
    out.removeAtoms(drop);
    return out;
}

std::vector<StereoLabel> stereoLabels(const Document& doc) {
    std::vector<StereoLabel> out;
    if (doc.atoms.empty()) return out;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return out;
    try {
        RDKit::CIPLabeler::assignCIPLabels(*mol);
    } catch (...) {
        return out;  // e.g. the labeller's work limit on huge symmetric molecules
    }
    const int n = int(doc.atoms.size());
    std::string code;
    for (const auto* a : mol->atoms())
        if (int(a->getIdx()) < n && a->getPropIfPresent(RDKit::common_properties::_CIPCode, code))
            out.push_back({int(a->getIdx()), -1, QString::fromStdString(code)});
    for (const auto* b : mol->bonds()) {
        const int i = int(b->getBeginAtomIdx()), j = int(b->getEndAtomIdx());
        if (i < n && j < n && b->getPropIfPresent(RDKit::common_properties::_CIPCode, code))
            out.push_back({-1, doc.bondBetween(i, j), QString::fromStdString(code)});
    }
    return out;
}

std::vector<Problem> checkStructure(const Document& doc) {
    std::vector<Problem> out;
    const int n = int(doc.atoms.size());
    auto name = [&](int i) {
        const Atom& a = doc.atoms[i];
        return a.label.isEmpty() ? QString::fromStdString(symbol(a.z)) : a.label;
    };
    const auto info = atomInfo(doc);
    for (int i = 0; i < n; ++i)
        if (info[i].valenceError)
            out.push_back({QObject::tr("Valence error: %1 has too many bonds").arg(name(i)), {i}});
    for (int i = 0; i < n; ++i)
        if (!doc.atoms[i].label.isEmpty() && !abbreviationHead(doc.atoms[i].label))
            out.push_back({QObject::tr("Unknown label \"%1\": drawn, but treated as an unknown group")
                               .arg(doc.atoms[i].label),
                           {i}});
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (QLineF(doc.atoms[i].pos, doc.atoms[j].pos).length() < 0.3 * kBondLength)
                out.push_back({QObject::tr("Overlapping atoms: %1 and %2").arg(name(i), name(j)), {i, j}});

    auto mol = toRDKit(doc);
    if (n && perceive(*mol)) {
        QSet<int> centres;
        for (const auto& s : RDKit::Chirality::findPotentialStereo(*mol)) {
            if (s.type != RDKit::Chirality::StereoType::Atom_Tetrahedral || int(s.centeredOn) >= n) continue;
            centres.insert(int(s.centeredOn));
            if (s.specified == RDKit::Chirality::StereoSpecified::Unspecified)
                out.push_back({QObject::tr("Stereocentre %1 has no wedge or hash").arg(name(int(s.centeredOn))),
                               {int(s.centeredOn)}});
        }
        for (const Bond& b : doc.bonds)
            if ((b.stereo == BondStereo::Wedge || b.stereo == BondStereo::Hash) && !centres.contains(b.a))
                out.push_back({QObject::tr("Wedge or hash starts at %1, which is not a stereocentre").arg(name(b.a)),
                               {b.a, b.b}});
    }
    return out;
}

std::vector<std::vector<int>> aromaticRings(const Document& doc) {
    std::vector<std::vector<int>> out;
    const int n = int(doc.atoms.size());
    if (!n) return out;
    auto mol = toRDKit(doc);
    if (!perceive(*mol)) return out;
    for (const auto& r : mol->getRingInfo()->atomRings())
        if (std::all_of(r.begin(), r.end(), [&](int i) { return i < n && mol->getAtomWithIdx(i)->getIsAromatic(); }))
            out.push_back(r);
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
