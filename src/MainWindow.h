#pragma once
#include <QMainWindow>

class Canvas;
class QLabel;
class QUndoStack;
struct Document;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();
    bool openFile(const QString& path);

protected:
    void closeEvent(QCloseEvent* e) override;

private:
    void buildTools();
    void buildMenus();
    bool saveTo(const QString& path);
    bool save();
    bool saveAs();
    bool maybeSave();
    void exportImage();
    void importSmiles();
    void copy();
    void paste();
    void updateTitle();
    void updateInfo();

    QUndoStack* undo_;
    Canvas* canvas_;
    QString path_;
    QLabel* info_;
};
