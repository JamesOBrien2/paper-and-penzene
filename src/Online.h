#pragma once

#include <QByteArray>
#include <QString>
#include <QUrl>

// The app's only network use, always on the user's say-so: PubChem lookups, and
// the check for a newer release (by hand, or weekly if they turn it on).
namespace online {
// The latest release's version tag (e.g. "v0.9.0") and page, from GitHub's API JSON.
struct Release {
    QString tag, url;
};
QUrl latestReleaseUrl();
Release parseRelease(const QByteArray& json);
bool isNewer(const QString& tag, const QString& current);  // "v0.9.0" vs "0.8.0"
// One release's section of CHANGELOG.md (its "## 0.9.0 (date)" heading to the next); empty if none.
QString releaseNotes(const QString& changelog, const QString& version);
}  // namespace online

// Optional online lookups against PubChem, only when the user asks for one.
// ponytail: only compounds PubChem already knows resolve; a local namer (OPSIN
// needs a JVM, STOUT an ML runtime) would lift that, at a heavy dependency cost.
namespace pubchem {

QUrl nameToSmilesUrl(const QString& name);
// POSTed with smilesToNameForm: SMILES can hold '/' and '#', which a URL path can't.
QUrl smilesToNameUrl();
QByteArray smilesToNameForm(const QString& smiles);
// The first record's `key` from a PUG REST property table; empty when absent.
QString property(const QByteArray& json, const QString& key);
// Blocking GET, or POST of `form` (with a timeout), of `url`'s `key`; empty with `error` set on failure.
QString fetch(const QUrl& url, const QString& key, QString* error, const QByteArray& form = {});

}  // namespace pubchem
