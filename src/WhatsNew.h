#pragma once
class QString;
class QWidget;

// What's New in Penzene <version>: the logo, the release's highlights as cards with an
// icon, the smaller changes as text (from that version's section of CHANGELOG.md).
void showWhatsNewDialog(QWidget* parent, const QString& changelog, const QString& version);
