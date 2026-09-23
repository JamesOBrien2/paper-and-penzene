#pragma once
#include <QMainWindow>
#include <functional>
#include <vector>

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
    void applyTheme(const QString& name);

    QUndoStack* undo_;
    Canvas* canvas_;
    QString path_;
    QLabel* info_;
    std::vector<std::pair<QAction*, std::function<QIcon()>>> icons_;
};
