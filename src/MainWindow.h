#pragma once
#include <QMainWindow>

class Canvas;
class QUndoStack;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();

private:
    void buildTools();
    void buildMenus();

    QUndoStack* undo_;
    Canvas* canvas_;
};
