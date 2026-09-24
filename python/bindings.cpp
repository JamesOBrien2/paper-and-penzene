// Python bindings over penzene_core: the same engine and renderer as the app.
#include "Chem.h"
#include "Edit.h"
#include "Render.h"

#include <QBuffer>
#include <QGuiApplication>
#include <QImage>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// Fonts and painting need a QGuiApplication. Reuse the host's if there is one
// (e.g. PySide6); otherwise make an offscreen one that lives as long as the process.
void ensureApp() {
    if (QCoreApplication::instance()) return;
    if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM")) qputenv("QT_QPA_PLATFORM", "offscreen");
    static int argc = 1;
    static char name[] = "penzene";
    static char* argv[] = {name, nullptr};
    new QGuiApplication(argc, argv);  // ponytail: never freed; the process owns it until exit
}

QString qs(const std::string& s) { return QString::fromStdString(s); }

Document fromSmiles(const std::string& smiles) {
    auto doc = chem::fromSmiles(smiles);
    if (!doc) throw nb::value_error(("not a valid SMILES: " + smiles).c_str());
    return *doc;
}

Document readPath(const std::string& path) {
    auto doc = chem::readFile(qs(path));
    if (!doc) throw nb::value_error(("cannot read " + path).c_str());
    return *doc;
}

void save(const Document& doc, const std::string& path) {
    const bool penz = qs(path).endsWith(".penz", Qt::CaseInsensitive);
    const QByteArray data = penz ? doc.toJson() : QByteArray::fromStdString(chem::toMolBlock(doc));
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f || std::fwrite(data.constData(), 1, data.size(), f) != size_t(data.size()))
        throw nb::value_error(("cannot write " + path).c_str());
    std::fclose(f);
}

std::string symbol(const Atom& a) { return a.label.isEmpty() ? chem::symbol(a.z) : a.label.toStdString(); }

}  // namespace

NB_MODULE(_penzene, m) {
    ensureApp();
    m.doc() = "Penzene: ChemDraw-style 2D structures, drawn by the same engine as the app.";
    m.attr("__version__") = PENZENE_VERSION;

    nb::class_<Atom>(m, "Atom")
        .def_prop_ro("symbol", &symbol, "Element symbol, or the abbreviation (Boc, OMe…) if it has one")
        .def_prop_ro("x", [](const Atom& a) { return a.pos.x(); })
        .def_prop_ro("y", [](const Atom& a) { return a.pos.y(); }, "Points, y down")
        .def_ro("charge", &Atom::charge)
        .def("__repr__", [](const Atom& a) { return "<Atom " + symbol(a) + ">"; });

    nb::class_<Bond>(m, "Bond")
        .def_ro("a", &Bond::a)
        .def_ro("b", &Bond::b)
        .def_ro("order", &Bond::order)
        .def("__repr__", [](const Bond& b) {
            return "<Bond " + std::to_string(b.a) + "-" + std::to_string(b.b) + " order " + std::to_string(b.order) + ">";
        });

    nb::class_<Document>(m, "Document")
        .def(nb::init<>())
        .def_ro("atoms", &Document::atoms)
        .def_ro("bonds", &Document::bonds)
        .def_prop_rw(
            "style", [](const Document& d) { return drawingStyle(d.style).name.toStdString(); },
            [](Document& d, const std::string& name) {
                const DrawingStyle& s = drawingStyle(qs(name));
                if (s.name != qs(name)) throw nb::value_error(("unknown drawing style: " + name).c_str());
                d.style = s.name == drawingStyles()[0].name ? QString() : s.name;
            },
            "Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC'")
        .def(
            "add_atom",
            [](Document& d, const std::string& label, double x, double y) {
                int i = d.addAtom({x, y});
                if (!edit::applyLabel(d, i, qs(label))) {
                    d.atoms.pop_back();
                    throw nb::value_error(("not an element, abbreviation or SMILES: " + label).c_str());
                }
                return i;
            },
            "label"_a = "C", "x"_a = 0.0, "y"_a = 0.0,
            "Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index.")
        .def(
            "add_bond",
            [](Document& d, int a, int b, int order) {
                const int n = int(d.atoms.size());
                if (a < 0 || b < 0 || a >= n || b >= n || a == b) throw nb::index_error("no such atoms");
                if (order < 1 || order > 3) throw nb::value_error("order must be 1, 2 or 3");
                edit::link(d, a, b, order);
            },
            "a"_a, "b"_a, "order"_a = 1)
        .def(
            "hotkeys",
            [](Document& d, int atom, const std::string& keys, std::optional<int> bond) {
                edit::Hotspot h{bond ? -1 : atom, bond.value_or(-1)};
                if (h.atom >= int(d.atoms.size()) || h.bond >= int(d.bonds.size()))
                    throw nb::index_error("no such atom or bond");
                for (QChar k : qs(keys)) {
                    h = edit::hotkey(d, h, QString(k));
                    if (!h.valid()) throw nb::value_error(("not a hotkey here: " + QString(k).toStdString()).c_str());
                }
                return std::make_tuple(h.atom, h.bond);
            },
            "atom"_a, "keys"_a, "bond"_a = nb::none(),
            "Type ChemDraw hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). "
            "Returns the final (atom, bond) hotspot.")
        .def(
            "clean", [](Document& d, std::vector<int> atoms) { d = chem::clean2D(d, atoms); }, "atoms"_a = std::vector<int>{},
            "Lay out afresh with RDKit; with atoms, only the molecules containing them.")
        .def("to_smiles", [](const Document& d) { return chem::toSmiles(d); })
        .def("to_molblock", [](const Document& d) { return chem::toMolBlock(d); })
        .def("to_inchi", [](const Document& d) { return chem::toInchi(d); })
        .def("to_inchikey", [](const Document& d) { return chem::toInchiKey(d); })
        .def("to_json", [](const Document& d) { return d.toJson().toStdString(); })
        .def("to_svg", [](const Document& d) { return renderSvg(d).toStdString(); })
        .def(
            "to_png",
            [](const Document& d, double dpi) {
                QByteArray png;
                QBuffer buf(&png);
                buf.open(QIODevice::WriteOnly);
                renderImage(d, dpi).save(&buf, "PNG");
                return nb::bytes(png.constData(), png.size());
            },
            "dpi"_a = 300)
        .def_prop_ro("formula", [](const Document& d) { auto p = chem::properties(d); return p ? p->formula : ""; })
        .def_prop_ro("mw", [](const Document& d) { auto p = chem::properties(d); return p ? p->mw : 0.0; })
        .def_prop_ro("exact_mass", [](const Document& d) { auto p = chem::properties(d); return p ? p->exactMass : 0.0; })
        .def("save", &save, "path"_a, "Write .penz (full fidelity) or MOL for any other extension.")
        .def(
            "export",
            [](const Document& d, const std::string& path) {
                if (!exportDocument(d, qs(path))) throw nb::value_error(("cannot export " + path).c_str());
            },
            "path"_a, "Write .svg, .png or .pdf, exactly as the app exports.")
        .def("_repr_svg_", [](const Document& d) { return renderSvg(d).toStdString(); })
        .def("__repr__", [](const Document& d) {
            auto p = chem::properties(d);
            return "<penzene.Document " + (p ? p->formula : std::to_string(d.atoms.size()) + " atoms") + ">";
        });

    m.def("from_smiles", &fromSmiles, "smiles"_a, "A structure from SMILES, laid out in 2D.");
    m.def("from_json", [](const std::string& json) {
        auto d = Document::fromJson(QByteArray::fromStdString(json));
        if (!d) throw nb::value_error("not a .penz document");
        return *d;
    });
    m.def("read", &readPath, "path"_a, "Open .penz, .mol/.sdf (first record) or ChemDraw .cdxml.");
    m.def("drawing_styles", [] {
        std::vector<std::string> out;
        for (const auto& s : drawingStyles()) out.push_back(s.name.toStdString());
        return out;
    });
}
