#!/usr/bin/env bash
# Builds the RDKit that release packages ship: the version in pixi.lock, without RDKit's own
# ChemDraw library (it bundles expatpp, MPL 1.0, which isn't GPL-compatible; Penzene reads and
# writes CDX with its own copy in third_party/chemdraw). Same options as the Flatpak.
#
#   pixi run bash cmake/build-rdkit.sh rdkit-prefix
#   export RDKit_ROOT=$PWD/rdkit-prefix      # then configure Penzene as usual
set -euo pipefail
export MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-11.0}  # as the wheels; ignored elsewhere
prefix=$(mkdir -p "$1" && cd "$1" && pwd)
work=${RUNNER_TEMP:-${TMPDIR:-/tmp}}/rdkit-build
rm -rf "$work" && mkdir -p "$work" && cd "$work"

fetch() {  # url sha256 file
    curl -fsSL --retry 3 -o "$3" "$1"
    if command -v sha256sum >/dev/null; then echo "$2  $3" | sha256sum -c -; else echo "$2  $3" | shasum -a 256 -c -; fi
}
fetch https://github.com/rdkit/rdkit/archive/refs/tags/Release_2026_03_6.tar.gz \
      d4d20b3b140237084694518aab34fdba6929d44bd7f720bce69329516abef663 rdkit.tar.gz
fetch https://github.com/IUPAC-InChI/InChI/releases/download/v1.07.3/INCHI-1-SRC.zip \
      b42d828b5d645bd60bc43df7e0516215808d92e5a46c28e12b1f4f75dfaae333 inchi.zip
fetch https://github.com/rareylab/RingDecomposerLib/archive/v1.1.3_rdkit.tar.gz \
      944b5816712a48bbf88aa25d4300ce11871ddf6e971218eac08f90ed2192f715 rdl.tar.gz
fetch https://github.com/aantron/better-enums/archive/c35576bed0295689540b39873126129adfa0b4c8.tar.gz \
      9b78dcef7f88d1345b6f25335bfbcba5f024b08990c7d6dc605b833f4128b8dd better-enums.tar.gz

mkdir rdkit better-enums
tar xzf rdkit.tar.gz --strip-components=1 -C rdkit
tar xzf better-enums.tar.gz --strip-components=1 -C better-enums
mkdir -p rdkit/External/INCHI-API/src rdkit/External/RingFamilies/RingDecomposerLib
python -c "import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" inchi.zip inchi
cp -R inchi/*/. rdkit/External/INCHI-API/src/
tar xzf rdl.tar.gz --strip-components=1 -C rdkit/External/RingFamilies/RingDecomposerLib

# Boost and Eigen from the pixi environment, so the libraries match the rest of the build.
cmake -G Ninja -S rdkit -B build -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_PREFIX_PATH="${CONDA_PREFIX:?run inside pixi}" \
    -DRDK_INSTALL_INTREE=OFF -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_CONTRIB=OFF -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_CHEMDRAW_SUPPORT=OFF \
    -DRDK_BUILD_COORDGEN_SUPPORT=OFF -DRDK_BUILD_MAEPARSER_SUPPORT=OFF -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_BUILD_CAIRO_SUPPORT=OFF -DRDK_BUILD_YAEHMOP_SUPPORT=OFF -DRDK_BUILD_FREESASA_SUPPORT=OFF \
    -DRDK_BUILD_AVALON_SUPPORT=OFF -DRDK_BUILD_PUBCHEMSHAPE_SUPPORT=OFF \
    -DRDK_USE_BOOST_SERIALIZATION=OFF -DRDK_USE_BOOST_IOSTREAMS=OFF \
    -DFETCHCONTENT_SOURCE_DIR_BETTER_ENUMS="$work/better-enums" -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build build
cmake --install build
echo "RDKit (no ChemDraw library) in $prefix"
