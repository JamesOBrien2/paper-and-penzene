// Binary CDX ⇄ CDXML through Revvity's ChemDraw file library (third_party/chemdraw), so
// every platform reads and writes .cdx the same way: CDXML is what the rest of Penzene speaks.
#include "Chem.h"

#include "chemdraw/CDXMLParser.h"
#include "chemdraw/CDXStdObjects.h"

#include <sstream>

namespace chem {

QByteArray cdxToCdxml(const QByteArray& cdx) {
    try {
        std::istringstream in(std::string(cdx.constData(), size_t(cdx.size())), std::ios::binary);
        CDXistream input(in);
        std::unique_ptr<CDXDocument> doc(CDXReadDocFromStorage(input, true));
        if (!doc) return {};
        std::ostringstream out;
        out << kCDXML_HeaderString;
        XMLDataSink sink(out);
        doc->XMLWrite(sink);
        return QByteArray::fromStdString(out.str());
    } catch (...) {
        return {};
    }
}

QByteArray cdxmlToCdx(const QByteArray& cdxml) {
    try {
        CDXMLParser parser;
        if (parser.XML_Parse(cdxml.constData(), int(cdxml.size()), true) != XML_STATUS_OK) return {};
        std::unique_ptr<CDXDocument> doc = parser.ReleaseDocument();
        if (!doc) return {};
        std::ostringstream out(std::ios::binary);
        CDXostream sink(out);
        CDXWriteDocToStorage(doc.get(), sink);
        return QByteArray::fromStdString(out.str());
    } catch (...) {
        return {};
    }
}

}  // namespace chem
