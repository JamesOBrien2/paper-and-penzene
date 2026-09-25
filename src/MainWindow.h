#pragma once
#include <QMainWindow>
#include <functional>
#include <vector>

#ifdef Q_OS_MACOS
#include <QUtiMimeConverter>
// ChemDraw's pasteboard type (binary CDX) as chemical/x-cdx.
struct ChemDrawPasteboard : QUtiMimeConverter {
    ChemDrawPasteboard();
    QString mimeForUti(const QString& uti) const override;
    QString utiForMime(const QString& mime) const override;
    QVariant convertToMime(const QString&, const QList<QByteArray>& data, const QString&) const override;
    QList<QByteArray> convertFromMime(const QString&, const QVariant& data, const QString&) const override;
};
#endif

class Canvas;
class QLabel;
class QPrinter;
class QUndoStack;
struct Document;

// The drawing at its export size, centred on the page; shrunk to fit if it's bigger.
bool printDocument(QPrinter& printer, const Document& doc);

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow() override;
    bool openFile(const QString& path);
    // Crash recovery: offer the autosaved document left behind by a crash, if any.
    void offerRecovery();
    // A newer release on GitHub? By hand (Help menu) it always answers; quietly (the weekly
    // check, off unless turned on in Preferences) it only speaks up if there is one.
    void checkForUpdates(bool quietly);
    void maybeCheckForUpdates();  // at startup: weekly, if turned on
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
    bool saveTo(const QString& path, bool v3000 = false);
    bool save();
    bool saveAs();
    bool maybeSave();
    void exportImage();
    void importSmiles();
    void importName();
    void print();
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
    class QDockWidget* templateDock_;
    class QTreeWidget* templates_;
    void fillTemplates();
    void insertTemplate(class QTreeWidgetItem* item);
    void saveTemplate();
    QLabel* profile_;
    QString profileText_;  // plain-text copy of the panel, for the Copy button
    std::vector<std::pair<QAction*, std::function<QIcon()>>> icons_;
};
