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
    ~MainWindow() override;
    bool openFile(const QString& path);
    // Crash recovery: offer the autosaved document left behind by a crash, if any.
    void offerRecovery();
    void autosave();  // writes unsaved changes to autosavePath() (every minute)
    static QString autosavePath();
    QStringList recentFiles() const;
    void showPreferences();  // theme, default style, export resolution and background
    QWidget* checkStructure();  // lists problems; clicking one selects its atoms

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
    void updateProfile();  // the Properties panel, while it's visible
    void applyTheme(const QString& name);
    void remember(const QString& path);  // most recent first, at most 10

    QUndoStack* undo_;
    Canvas* canvas_;
    QString path_;
    QLabel* info_;
    class QActionGroup* themeGroup_ = nullptr;
    class QDockWidget* profileDock_;
    QLabel* profile_;
    QString profileText_;  // plain-text copy of the panel, for the Copy button
    std::vector<std::pair<QAction*, std::function<QIcon()>>> icons_;
};
