// Penzene: see expatpp.h.
#include "expatpp.h"

expatpp::expatpp() : m_parser(XML_ParserCreate(nullptr)) {
    XML_SetUserData(m_parser, this);
    XML_SetElementHandler(m_parser, onStart, onEnd);
    XML_SetCharacterDataHandler(m_parser, onText);
    XML_SetNotStandaloneHandler(m_parser, onNotStandalone);
}

expatpp::~expatpp() { XML_ParserFree(m_parser); }

XML_Status expatpp::XML_Parse(const char* data, int length, int isFinal) {
    return ::XML_Parse(m_parser, data, length, isFinal);
}

XML_Error expatpp::XML_GetErrorCode() const { return ::XML_GetErrorCode(m_parser); }

XML_Size expatpp::XML_GetCurrentLineNumber() const { return ::XML_GetCurrentLineNumber(m_parser); }

void XMLCALL expatpp::onStart(void* self, const XML_Char* name, const XML_Char** atts) {
    static_cast<expatpp*>(self)->startElement(name, atts);
}
void XMLCALL expatpp::onEnd(void* self, const XML_Char* name) { static_cast<expatpp*>(self)->endElement(name); }
void XMLCALL expatpp::onText(void* self, const XML_Char* s, int len) { static_cast<expatpp*>(self)->charData(s, len); }
int XMLCALL expatpp::onNotStandalone(void* self) { return static_cast<expatpp*>(self)->notStandaloneHandler(); }
