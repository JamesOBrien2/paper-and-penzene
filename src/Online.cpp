#include "Online.h"

#include <QEventLoop>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QRegularExpression>
#include <QTimer>
#include <QVersionNumber>

namespace online {

QUrl latestReleaseUrl() { return QUrl("https://api.github.com/repos/JamesOBrien2/penzene/releases/latest"); }

Release parseRelease(const QByteArray& json) {
    const auto o = QJsonDocument::fromJson(json).object();
    return {o["tag_name"].toString(), o["html_url"].toString()};
}

bool isNewer(const QString& tag, const QString& current) {
    const auto latest = QVersionNumber::fromString(QString(tag).remove(QRegularExpression("^v")));
    return !latest.isNull() && latest > QVersionNumber::fromString(current);
}

QString releaseNotes(const QString& changelog, const QString& version) {
    const QRegularExpression heading("^## " + QRegularExpression::escape(version) + "(?=[ \\t]|$)[^\n]*\n",
                                     QRegularExpression::MultilineOption);
    const auto m = heading.match(changelog);
    if (!m.hasMatch()) return {};
    const qsizetype next = changelog.indexOf("\n## ", m.capturedEnd());
    return changelog.mid(m.capturedEnd(), next < 0 ? -1 : next - m.capturedEnd()).trimmed();
}

ReleaseNotes parseReleaseNotes(const QString& section) {
    static const QRegularExpression highlight(R"(^\*\*(.+?)\*\*:?\s*(.*?)\s*(?:<!--\s*icon:\s*([\w-]+)\s*-->)?\s*$)");
    static const QRegularExpression comment(R"(\s*<!--.*?-->)");
    ReleaseNotes notes;
    for (QString line : section.split('\n')) {
        line = line.trimmed();
        if (!line.startsWith("- ")) continue;
        line = line.mid(2);
        if (const auto m = highlight.match(line); m.hasMatch())
            notes.highlights.push_back({m.captured(3).isEmpty() ? "sparkles" : m.captured(3), m.captured(1), m.captured(2)});
        else
            notes.others << QString(line).remove(comment).remove("**");
    }
    return notes;
}

}  // namespace online

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
