#include "Canvas.h"
#include "Chem.h"
#include "MainWindow.h"

#include <QApplication>
#include <QFile>
#include <QIcon>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--version")) {
            std::printf("penzene %s\n", PENZENE_VERSION);
            return 0;
        }
    // In an AppImage, conda's Qt doesn't find the bundled plugins on its own.
    if (const char* appdir = std::getenv("APPDIR"); appdir && !std::getenv("QT_PLUGIN_PATH"))
        qputenv("QT_PLUGIN_PATH", QByteArray(appdir) + "/usr/plugins");
    QApplication app(argc, argv);
    // Headless: penzene --render <SMILES|file.penz> <out.svg|png|pdf>
    if (argc == 4 && !std::strcmp(argv[1], "--render")) {
        QFile in(QString::fromLocal8Bit(argv[2]));
        auto doc = in.open(QIODevice::ReadOnly) ? Document::fromJson(in.readAll()) : chem::fromSmiles(argv[2]);
        if (!doc || !exportDocument(*doc, argv[3])) {
            std::fprintf(stderr, "penzene: could not render %s\n", argv[2]);
            return 1;
        }
        return 0;
    }
    QApplication::setApplicationName("Penzene");
    QApplication::setWindowIcon(QIcon(":/logo.svg"));
    MainWindow w;
    if (argc == 2) w.openFile(QString::fromLocal8Bit(argv[1]));
    w.show();
    return app.exec();
}
