#include "Canvas.h"
#include "Chem.h"
#include "MainWindow.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QIcon>
#include <QRegularExpression>
#include <QSet>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

using chem::Record;

// One input: a SMILES string, a multi-record file (SDF, .smi, .inchi) or any
// single file Penzene opens.
std::vector<Record> records(const QString& in) {
    const QFileInfo info(in);
    if (!info.exists()) return {{"structure", chem::fromSmiles(in.toStdString())}};
    const QString ext = info.suffix().toLower();
    if (ext == "smi" || ext == "sdf" || ext == "inchi") return chem::readRecords(in);
    return {{info.completeBaseName(), chem::readFile(in)}};
}

// penzene --render IN [IN...] (OUT | --out DIR) [--format svg|png|pdf] [--drawing-style NAME] [--clean]
int render(const QStringList& args) {
    QCommandLineParser p;
    p.addOption({"render", "Render structures without opening a window."});
    p.addOption({{"o", "out"}, "Directory for the rendered files.", "dir"});
    p.addOption({"format", "svg, png or pdf (default svg) when writing to --out.", "ext", "svg"});
    // Not --style: QApplication claims that for widget styles.
    p.addOption({"drawing-style", "ACS 1996 (default), JDP or RSC.", "name"});
    p.addOption({"clean", "Lay out each structure afresh with RDKit."});
    p.addPositionalArgument("inputs", "SMILES, .smi, .sdf, .inchi, .mol, .penz or .cdxml");
    if (!p.parse(args)) {
        std::fprintf(stderr, "penzene: %s\n", qPrintable(p.errorText()));
        return 2;
    }
    QStringList inputs = p.positionalArguments();
    QString single;  // legacy form: penzene --render IN OUT.ext
    if (!p.isSet("out") && inputs.size() == 2) single = inputs.takeLast();
    if (inputs.isEmpty() || (!p.isSet("out") && single.isEmpty())) {
        std::fprintf(stderr, "usage: penzene --render IN... (OUT.svg|png|pdf | --out DIR) "
                             "[--format svg|png|pdf] [--drawing-style NAME] [--clean]\n");
        return 2;
    }
    if (p.isSet("out")) QDir().mkpath(p.value("out"));
    int failed = 0;
    QSet<QString> used;
    for (const QString& in : inputs)
        for (auto& [name, doc] : records(in)) {
            QString path = single;
            if (path.isEmpty()) {
                const QString safe = QString(name).replace(QRegularExpression("[^A-Za-z0-9._-]+"), "_");
                QString stem = safe;
                for (int k = 2; used.contains(stem); ++k) stem = safe + QString("-%1").arg(k);
                used.insert(stem);
                path = QDir(p.value("out")).filePath(stem + "." + p.value("format"));
            }
            if (doc && p.isSet("clean")) doc = chem::clean2D(*doc);
            if (doc && p.isSet("drawing-style")) doc->style = drawingStyle(p.value("drawing-style")).name;
            if (!doc || !exportDocument(*doc, path)) {
                std::fprintf(stderr, "penzene: could not render %s\n", qPrintable(name));
                ++failed;
            } else {
                std::printf("%s\n", qPrintable(path));  // one line per file, for scripts
            }
        }
    return failed ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--version")) {
            std::printf("penzene %s\n", PENZENE_BUILD);
            return 0;
        } else if (!std::strcmp(argv[i], "--help") || !std::strcmp(argv[i], "-h")) {
            std::printf("usage: penzene [FILE]\n"
                        "       penzene --render IN... (OUT.svg|png|pdf | --out DIR) "
                        "[--format svg|png|pdf] [--drawing-style NAME] [--clean]\n"
                        "       penzene --version\n");
            return 0;
        }
    // In an AppImage, conda's Qt doesn't find the bundled plugins on its own.
    if (const char* appdir = std::getenv("APPDIR"); appdir && !std::getenv("QT_PLUGIN_PATH"))
        qputenv("QT_PLUGIN_PATH", QByteArray(appdir) + "/usr/plugins");
    QApplication app(argc, argv);
    if (argc > 1 && !std::strcmp(argv[1], "--render")) return render(app.arguments());
    QApplication::setApplicationName("Penzene");
    QApplication::setOrganizationName("Penzene");
    QApplication::setWindowIcon(QIcon(":/logo.svg"));
    MainWindow w;
    if (argc == 2) w.openFile(QString::fromLocal8Bit(argv[1]));
    w.show();
    w.offerRecovery();  // after a crash, the last unsaved drawing
    w.maybeCheckForUpdates();  // only if turned on in Preferences
    return app.exec();
}
