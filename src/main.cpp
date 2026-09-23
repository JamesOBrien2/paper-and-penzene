#include "Chem.h"

#include <QApplication>
#include <QMainWindow>
#include <cstdio>
#include <cstring>

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--version")) {
            std::printf("penzene %s\n", PENZENE_VERSION);
            return 0;
        }
    QApplication app(argc, argv);
    QApplication::setApplicationName("Paper & Penzene");
    QMainWindow w;
    w.resize(1000, 700);
    w.show();
    return app.exec();
}
