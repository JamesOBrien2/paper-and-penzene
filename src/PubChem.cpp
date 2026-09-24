#include "PubChem.h"

#include <QEventLoop>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QTimer>

namespace pubchem {

static const QString kBase = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/";

QUrl nameToSmilesUrl(const QString& name) {
    return QUrl(kBase + "name/" + QUrl::toPercentEncoding(name.trimmed()) + "/property/SMILES/JSON");
}

QUrl smilesToNameUrl() { return QUrl(kBase + "smiles/property/IUPACName/JSON"); }

QByteArray smilesToNameForm(const QString& smiles) { return "smiles=" + QUrl::toPercentEncoding(smiles); }

QString property(const QByteArray& json, const QString& key) {
    const auto rows = QJsonDocument::fromJson(json).object()["PropertyTable"].toObject()["Properties"].toArray();
    return rows.isEmpty() ? QString() : rows[0].toObject()[key].toString();
}

QString fetch(const QUrl& url, const QString& key, QString* error, const QByteArray& form) {
    QNetworkAccessManager net;
    QNetworkRequest request(url);
    request.setTransferTimeout(15000);
    request.setHeader(QNetworkRequest::ContentTypeHeader, "application/x-www-form-urlencoded");
    QNetworkReply* reply = form.isEmpty() ? net.get(request) : net.post(request, form);
    QEventLoop wait;
    QObject::connect(reply, &QNetworkReply::finished, &wait, &QEventLoop::quit);
    wait.exec(QEventLoop::ExcludeUserInputEvents);  // no closing the window mid-request
    const QByteArray body = reply->readAll();
    const QString value = property(body, key);
    if (value.isEmpty())
        *error = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt() == 404
                     ? QObject::tr("PubChem has no match.")
                     : reply->errorString();
    reply->deleteLater();
    return value;
}

}  // namespace pubchem
