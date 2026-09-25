#pragma once
#include "Document.h"

#include <QString>
#include <utility>
#include <vector>

// The template library: built-in structures by category, and the user's own
// (saved selections, one .penz each in the app's data folder).
struct Template {
    QString category, name, smiles;
    Document (*draw)() = nullptr;  // a drawing (projections) instead of laying out `smiles`
};
Document templateDocument(const Template& t);  // empty if it doesn't parse
const std::vector<Template>& builtinTemplates();
std::vector<std::pair<QString, Document>> userTemplates();  // by name
bool saveUserTemplate(const QString& name, const Document& doc);  // replaces one of the same name
bool removeUserTemplate(const QString& name);
std::vector<std::pair<QString, Document>> exampleDocuments();  // for the welcome screen
