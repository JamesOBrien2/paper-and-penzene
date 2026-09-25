// Penzene: a minimal stand-in for the expatpp C++ wrapper's interface, written for
// Penzene over expat (MIT) because expatpp itself is MPL 1.0, which the GPL can't take.
// Only what this library uses: parse a buffer, and virtual element and text callbacks.
#pragma once
#include <expat.h>

class expatpp {
public:
    expatpp();
    virtual ~expatpp();
    expatpp(const expatpp&) = delete;
    expatpp& operator=(const expatpp&) = delete;

    XML_Status XML_Parse(const char* data, int length, int isFinal);
    XML_Error XML_GetErrorCode() const;
    XML_Size XML_GetCurrentLineNumber() const;
    operator XML_Parser() const { return m_parser; }

protected:
    virtual void startElement(const XML_Char*, const XML_Char**) {}
    virtual void endElement(const XML_Char*) {}
    virtual void charData(const XML_Char*, int) {}
    virtual int notStandaloneHandler() { return 1; }

private:
    XML_Parser m_parser;
    static void XMLCALL onStart(void* self, const XML_Char* name, const XML_Char** atts);
    static void XMLCALL onEnd(void* self, const XML_Char* name);
    static void XMLCALL onText(void* self, const XML_Char* s, int len);
    static int XMLCALL onNotStandalone(void* self);
};
