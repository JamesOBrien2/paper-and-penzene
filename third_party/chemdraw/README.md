# ChemDraw file format library

The CDX/CDXML object model and readers/writers released by Revvity (CambridgeSoft) under the
BSD 3-Clause license (`chemdraw/license.txt`), from https://github.com/Glysade/chemdraw v1.0.14,
the version RDKit 2026.03 uses. Penzene builds it as a private static library to read and write
binary CDX on every platform.

Changes for Penzene:
- `CoreChemistryAPI.h`: a static build exports nothing on any platform (marked in the file).
- `expatpp.h` / `expatpp.cpp`: Penzene's own minimal replacement for the expatpp wrapper over
  expat (MIT). The original expatpp is MPL 1.0, which isn't GPL-compatible, so it isn't included.
- `XMLDoc.cpp` (unused) is left out.
- `UTF8Iterator.h`: fixed (its length and position were wrong, garbling text on writing).
- `CDXUnicode.cpp` includes CoreFoundation on macOS; `XMLParser.cpp` calls the replacement
  expatpp's constructor. Both are marked "Modified for Penzene".
