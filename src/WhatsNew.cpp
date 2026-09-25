#include "WhatsNew.h"
#include "Online.h"

#include <QApplication>
#include <QDialog>
#include <QFile>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QSvgRenderer>
#include <QVBoxLayout>

namespace {

struct Tokens {
    QString page, card, border, text, secondary, accent, accentBg;
};

// The lab notebook look: warm paper, one teal accent; its dark twin when the app is dark.
Tokens tokens(const QWidget* w) {
    if (w->palette().color(QPalette::Window).lightness() < 128)
        return {"#22211F", "#2C2C2A", "#444441", "#F1EFE8", "#B4B2A9", "#5DCAA5", "#0B3B30"};
    return {"#FBF8F1", "#FFFFFF", "#E4E1D6", "#2C2C2A", "#5F5E5A", "#0F6E56", "#E1F5EE"};
}

// An SVG from the resources, drawn in `color` (the icons use currentColor).
QPixmap svg(const QString& path, QColor color, int size, qreal ratio) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) return {};
    QByteArray data = f.readAll();
    data.replace("currentColor", color.name().toUtf8());
    QSvgRenderer renderer(data);
    QPixmap pm(QSize(size, size) * ratio);
    pm.setDevicePixelRatio(ratio);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    renderer.render(&p, QRectF(0, 0, size, size));
    return pm;
}

}  // namespace

void showWhatsNewDialog(QWidget* parent, const QString& changelog, const QString& version) {
    const auto notes = online::parseReleaseNotes(online::releaseNotes(changelog, version));
    QDialog dialog(parent);
    dialog.setObjectName("whatsNew");
    dialog.setWindowTitle(QObject::tr("What's New"));
    const Tokens t = tokens(&dialog);
    const qreal ratio = dialog.devicePixelRatioF();
    dialog.setStyleSheet(QString(R"(
        QDialog#whatsNew { background: %1; }
        QLabel { color: %4; background: transparent; }
        QLabel#title { font-size: 20px; font-weight: 600; }
        QLabel#subtitle, QLabel#detail, QLabel#others { color: %5; }
        QLabel#section { color: %6; font-weight: 600; }
        QFrame#highlight { background: %2; border: 1px solid %3; border-radius: 14px; }
        QLabel#badge { background: %7; border-radius: 18px; }
        QPushButton#continue { background: %6; color: %1; border: none; border-radius: 10px; padding: 7px 22px; font-weight: 600; }
        QPushButton#continue:hover { background: %4; }
    )").arg(t.page, t.card, t.border, t.text, t.secondary, t.accent, t.accentBg));

    auto* layout = new QVBoxLayout(&dialog);
    layout->setContentsMargins(28, 24, 28, 20);
    layout->setSpacing(14);

    auto* logo = new QLabel;
    logo->setPixmap(svg(":/logo.svg", QColor(t.accent), 56, ratio));
    logo->setAlignment(Qt::AlignCenter);
    auto* title = new QLabel(QObject::tr("What's new in Penzene %1").arg(version));
    title->setObjectName("title");
    title->setAlignment(Qt::AlignCenter);
    layout->addWidget(logo);
    layout->addWidget(title);

    if (!notes.highlights.empty()) {
        auto* grid = new QGridLayout;
        grid->setSpacing(10);
        const int n = int(notes.highlights.size());
        for (int k = 0; k < n; ++k) {
            const auto& h = notes.highlights[k];
            auto* card = new QFrame;
            card->setObjectName("highlight");
            auto* row = new QHBoxLayout(card);
            row->setContentsMargins(12, 12, 12, 12);
            row->setSpacing(12);
            auto* badge = new QLabel;
            badge->setObjectName("badge");
            badge->setFixedSize(36, 36);
            badge->setAlignment(Qt::AlignCenter);
            QPixmap icon = svg(":/whatsnew/" + h.icon + ".svg", QColor(t.accent), 20, ratio);
            if (icon.isNull()) icon = svg(":/whatsnew/sparkles.svg", QColor(t.accent), 20, ratio);
            badge->setPixmap(icon);
            row->addWidget(badge, 0, Qt::AlignTop);
            auto* text = new QVBoxLayout;
            text->setSpacing(2);
            auto* name = new QLabel(h.title);
            name->setStyleSheet("font-weight: 600;");
            auto* detail = new QLabel(h.detail);
            detail->setObjectName("detail");
            detail->setWordWrap(true);
            text->addWidget(name);
            text->addWidget(detail);
            text->addStretch();
            row->addLayout(text, 1);
            // Two columns; an odd one out takes the whole last row.
            const bool last = k == n - 1 && n % 2;
            grid->addWidget(card, k / 2, k % 2, 1, last ? 2 : 1);
        }
        grid->setColumnStretch(0, 1), grid->setColumnStretch(1, 1);
        layout->addLayout(grid);
    }

    if (!notes.others.isEmpty()) {
        auto* section = new QLabel(QObject::tr("Also in this release"));
        section->setObjectName("section");
        auto* others = new QLabel("• " + notes.others.join("<br>• "));
        others->setObjectName("others");
        others->setTextFormat(Qt::RichText);
        others->setWordWrap(true);
        layout->addWidget(section);
        layout->addWidget(others);
    }
    if (notes.highlights.empty() && notes.others.isEmpty()) layout->addWidget(new QLabel(QObject::tr("No notes for this version.")));

    auto* footer = new QHBoxLayout;
    auto* all = new QLabel(QString("<a style='color:%1' href='https://github.com/JamesOBrien2/penzene/blob/main/CHANGELOG.md'>%2</a>")
                               .arg(t.accent, QObject::tr("All releases")));
    all->setOpenExternalLinks(true);
    auto* ok = new QPushButton(QObject::tr("Continue"));
    ok->setObjectName("continue");
    ok->setDefault(true);
    QObject::connect(ok, &QPushButton::clicked, &dialog, &QDialog::accept);
    footer->addWidget(all);
    footer->addStretch();
    footer->addWidget(ok);
    layout->addSpacing(4);
    layout->addLayout(footer);
    dialog.setFixedWidth(620);
    dialog.exec();
}
