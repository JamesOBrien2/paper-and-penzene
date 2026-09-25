#!/usr/bin/env bash
# Fails if a package carries RDKit's ChemDraw library (MPL-1.0 expatpp inside), as a file or as
# a library that another one loads. Penzene's own copy is linked in statically and has no file.
#   bash cmake/check-no-chemdraw.sh dist/penzene.app
set -euo pipefail
found=$( { find "$@" -iname '*chemdraw*' -print; grep -rlaE 'libRDKit(RD)?ChemDraw' "$@" || true; } | sort -u)
if [ -n "$found" ]; then
    echo "RDKit's ChemDraw library is in the package:"
    echo "$found"
    exit 1
fi
echo "no RDKit ChemDraw library in: $*"
