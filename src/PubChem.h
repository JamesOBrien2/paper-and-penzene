#pragma once

#include <QByteArray>
#include <QString>
#include <QUrl>

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
