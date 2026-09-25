cmake_minimum_required(VERSION 3.24)
project(penzene VERSION 0.9.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# A wheel build (scikit-build-core sets SKBUILD) wants only the Python module.
if(SKBUILD)  # "2" under scikit-build-core, which option() won't read as ON
    set(app_default OFF)
    set(python_default ON)
else()
    set(app_default ON)
    set(python_default OFF)
endif()
option(PENZENE_APP "Build the desktop app and its tests" ${app_default})
option(PENZENE_PYTHON "Build the penzene Python module" ${python_default})
option(PENZENE_TESTS "Build the tests (needs Catch2)" ${app_default})

find_package(Qt6 REQUIRED COMPONENTS Gui Svg)
find_package(RDKit REQUIRED CONFIG)
message(STATUS "RDKit: ${RDKit_DIR}")
if(DEFINED ENV{RDKit_ROOT})  # release packages: the RDKit without its ChemDraw library (#233), never another
    cmake_path(SET rdkit_root NORMALIZE "$ENV{RDKit_ROOT}")
    cmake_path(IS_PREFIX rdkit_root "${RDKit_DIR}" NORMALIZE rdkit_from_root)
    if(NOT rdkit_from_root)
        message(FATAL_ERROR "RDKit_ROOT is ${rdkit_root}, but RDKit was found at ${RDKit_DIR}")
    endif()
endif()

# The engine: document, chemistry, editing and rendering. QtGui only, no
# widgets, so it can back a command line or Python module as well as the app.
add_library(penzene_core STATIC src/Document.cpp src/Chem.cpp src/Edit.cpp src/Render.cpp src/ChemDrawIO.cpp)
target_include_directories(penzene_core PUBLIC src)
target_compile_definitions(penzene_core PUBLIC PENZENE_VERSION="${PROJECT_VERSION}")
# Revvity's ChemDraw file library (BSD-3, third_party/chemdraw): binary CDX in and out on
# every platform. Private and static, so it can't clash with RDKit's own copy.
# The active conda/pixi environment's expat, not an older system one (Linux runners have both).
if(DEFINED ENV{CONDA_PREFIX})
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()
find_package(EXPAT REQUIRED)
find_package(Boost CONFIG REQUIRED)
file(GLOB chemdraw_sources third_party/chemdraw/chemdraw/*.cpp)
add_library(penzene_chemdraw STATIC ${chemdraw_sources})
target_include_directories(penzene_chemdraw PUBLIC third_party/chemdraw)
# UTF8_STD_STRING: text is UTF-8 throughout (as in Penzene), written with ChemDraw's UTF-8 text property.
target_compile_definitions(penzene_chemdraw PUBLIC CHEMDRAW_STATIC_LINK PLATFORM_LITTLEENDIAN UTF8_STD_STRING=1)
target_link_libraries(penzene_chemdraw PUBLIC EXPAT::EXPAT Boost::headers)
set_target_properties(penzene_chemdraw PROPERTIES CXX_VISIBILITY_PRESET hidden POSITION_INDEPENDENT_CODE ON AUTOMOC OFF)
# Its platform switches, set here for every file that includes its headers: Apple's own
# TargetConditionals.h isn't reliably included first, and nothing else sets the others.
if(APPLE)
    target_compile_definitions(penzene_chemdraw PUBLIC TARGET_OS_MAC=1)
    target_link_libraries(penzene_chemdraw PRIVATE "-framework CoreFoundation")
elseif(WIN32)
    target_compile_definitions(penzene_chemdraw PUBLIC TARGET_OS_WIN32=1 TARGET_API_LIB _WINDOWS)
else()  # Linux converts text encodings with Boost.Locale
    target_compile_definitions(penzene_chemdraw PUBLIC TARGET_OS_LINUX=1)
    find_package(Boost CONFIG REQUIRED COMPONENTS locale)
    target_link_libraries(penzene_chemdraw PRIVATE Boost::locale)
endif()
if(MSVC)
    target_compile_options(penzene_chemdraw PRIVATE /W0)
    target_compile_definitions(penzene_chemdraw PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
    target_compile_options(penzene_chemdraw PRIVATE -w)
endif()

target_link_libraries(penzene_core PUBLIC
    Qt6::Gui Qt6::Svg
    RDKit::SmilesParse RDKit::FileParsers RDKit::Depictor RDKit::GraphMol
    RDKit::Descriptors RDKit::RDInchiLib RDKit::CIPLabeler RDKit::DistGeomHelpers)
# Whole: RDKit's shared ChemDraw library (macOS, Linux) defines the same functions, and a
# plain archive link would take RDKit's copies instead of these (which carry Penzene's fixes).
target_link_libraries(penzene_core PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,penzene_chemdraw>")

if(PENZENE_APP)
find_package(Qt6 REQUIRED COMPONENTS Widgets Network PrintSupport Test)
# The app's widgets; tests link this too, so they exercise the same code.
add_library(penzene_ui STATIC src/Canvas.cpp src/MainWindow.cpp src/Online.cpp src/Templates.cpp src/WhatsNew.cpp)
target_link_libraries(penzene_ui PUBLIC penzene_core Qt6::Widgets Qt6::Network Qt6::PrintSupport)
qt_add_resources(penzene_ui "changelog" PREFIX "/" FILES CHANGELOG.md)  # What's New
# The logo (window icon, What's New) and What's New's icons, with the widgets so tests see them too.
file(GLOB whatsnew_icons RELATIVE ${CMAKE_SOURCE_DIR}/resources ${CMAKE_SOURCE_DIR}/resources/whatsnew/*.svg)
list(TRANSFORM whatsnew_icons PREPEND resources/)
qt_add_resources(penzene_ui "art" PREFIX "/" BASE resources FILES resources/logo.svg ${whatsnew_icons})
# Translations: each translations/penzene_<lang>.ts is compiled and embedded as :/i18n/penzene_<lang>.qm.
find_package(Qt6 REQUIRED COMPONENTS LinguistTools)
file(GLOB penzene_ts CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/translations/penzene_*.ts)
if(penzene_ts)
    qt_add_lrelease(TS_FILES ${penzene_ts} LRELEASE_TARGET penzene_lrelease NO_GLOBAL_TARGET
                    QM_OUTPUT_DIRECTORY i18n QM_FILES_OUTPUT_VARIABLE penzene_qm)
    qt_add_resources(penzene_ui "i18n" PREFIX "/i18n" BASE ${CMAKE_CURRENT_BINARY_DIR}/i18n FILES ${penzene_qm})
endif()
# The build shown in the title bar: git describe tells builds between releases
# apart (0.5.0-12-gabc1234; 0.6.0-dev-gabc1234 once bumped, before the tag);
# a tagged or git-less build is just the version.
set(PENZENE_BUILD ${PROJECT_VERSION})
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --match "v*"
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE described
                    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE failed)
    if(NOT failed)
        string(REGEX REPLACE "^v" "" described "${described}")
        if(described MATCHES "^${PROJECT_VERSION}(-|$)")
            set(PENZENE_BUILD "${described}")
        else()  # bumped but not tagged yet: this version, in development
            string(REGEX MATCH "g[0-9a-f]+(-dirty)?$" commit "${described}")
            set(PENZENE_BUILD "${PROJECT_VERSION}-dev")
            if(commit)
                set(PENZENE_BUILD "${PROJECT_VERSION}-dev-${commit}")
            endif()
        endif()
    endif()
endif()
target_compile_definitions(penzene_ui PUBLIC PENZENE_BUILD="${PENZENE_BUILD}")

set(PENZENE_ICNS resources/penzene.icns)
set_source_files_properties(${PENZENE_ICNS} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
add_executable(penzene WIN32 MACOSX_BUNDLE src/main.cpp ${PENZENE_ICNS} resources/penzene.rc)
target_link_libraries(penzene PRIVATE penzene_ui)
set_target_properties(penzene PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER io.github.jamesobrien2.penzene
    MACOSX_BUNDLE_BUNDLE_NAME "Penzene"
    MACOSX_BUNDLE_ICON_FILE penzene.icns
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION})
if(UNIX AND NOT APPLE)  # Linux desktop integration: Flatpak and distribution packages
    include(GNUInstallDirs)
    set(app_id io.github.jamesobrien2.penzene)
    install(TARGETS penzene RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES packaging/linux/${app_id}.desktop DESTINATION ${CMAKE_INSTALL_DATADIR}/applications)
    install(FILES packaging/linux/${app_id}.metainfo.xml DESTINATION ${CMAKE_INSTALL_DATADIR}/metainfo)
    install(FILES resources/logo.svg DESTINATION ${CMAKE_INSTALL_DATADIR}/icons/hicolor/scalable/apps RENAME ${app_id}.svg)
endif()
if(PENZENE_TESTS)
find_package(Catch2 3 REQUIRED)
enable_testing()
add_executable(penzene_tests tests/test_core.cpp tests/test_canvas.cpp)
# GuiPrivate: tests inject input the way the window system does (QWindowSystemInterface).
find_package(Qt6 REQUIRED COMPONENTS GuiPrivate)
target_link_libraries(penzene_tests PRIVATE penzene_ui Qt6::Test Qt6::GuiPrivate Catch2::Catch2WithMain)
# A tiny made-up language ("xx"), to test that a translation is found and applied.
qt_add_lrelease(TS_FILES tests/data/penzene_xx.ts LRELEASE_TARGET penzene_tests_lrelease NO_GLOBAL_TARGET
                QM_OUTPUT_DIRECTORY test-i18n QM_FILES_OUTPUT_VARIABLE test_qm)
qt_add_resources(penzene_tests "test_i18n" PREFIX "/i18n" BASE ${CMAKE_CURRENT_BINARY_DIR}/test-i18n FILES ${test_qm})
add_test(NAME core COMMAND penzene_tests)
add_test(NAME version COMMAND penzene --version)
# --help prints usage; it must not open a window that tries to read "--help" as a file.
add_test(NAME help COMMAND penzene --help)
set_tests_properties(help PROPERTIES TIMEOUT 20 ENVIRONMENT QT_QPA_PLATFORM=offscreen
                     PASS_REGULAR_EXPRESSION "usage: penzene")
# Batch CLI: a MOL file and a SMILES into a directory, RSC style, freshly laid out.
add_test(NAME render COMMAND penzene --render ${CMAKE_SOURCE_DIR}/tests/data/aspirin.mol CCO
         --out ${CMAKE_BINARY_DIR}/render-test --format png --drawing-style RSC --clean)
set_tests_properties(render PROPERTIES ENVIRONMENT QT_QPA_PLATFORM=offscreen
                     PASS_REGULAR_EXPRESSION "aspirin\\.png.*structure\\.png")
target_compile_definitions(penzene_tests PRIVATE PENZENE_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
endif()
endif()

if(PENZENE_PYTHON)
    # The module is a shared object linking the static core, so the core must be PIC.
    set_target_properties(penzene_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
    find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
    execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
    find_package(nanobind CONFIG REQUIRED)
    # Stable ABI: one build serves every Python from 3.12 on.
    nanobind_add_module(_penzene STABLE_ABI NB_STATIC python/bindings.cpp)
    target_link_libraries(_penzene PRIVATE penzene_core)
    target_compile_definitions(_penzene PRIVATE QT_NO_KEYWORDS)  # Qt's "slots" macro breaks Python.h
    if(SKBUILD)
        # Keep the conda library paths in the module so auditwheel/delocate/delvewheel
        # can find and bundle RDKit and Qt, plus Qt's offscreen platform plugin.
        set_target_properties(_penzene PROPERTIES INSTALL_RPATH_USE_LINK_PATH ON)
        install(TARGETS _penzene LIBRARY DESTINATION penzene)
        install(FILES $<TARGET_FILE:Qt6::QOffscreenIntegrationPlugin> DESTINATION penzene/plugins/platforms)
    else()
        # Importable from build/python during development and in ctest.
        set_target_properties(_penzene PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/penzene)
        configure_file(python/penzene/__init__.py ${CMAKE_BINARY_DIR}/python/penzene/__init__.py COPYONLY)
        if(PENZENE_APP AND PENZENE_TESTS)
            add_test(NAME python COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_python.py)
            set_tests_properties(python PROPERTIES ENVIRONMENT
                "PYTHONPATH=${CMAKE_BINARY_DIR}/python;PENZENE_TEST_DATA=${CMAKE_SOURCE_DIR}/tests/data")
            # The docs' generated pages (keys, Python API) match the code.
            add_test(NAME docs COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/docs/generate.py --check)
            set_tests_properties(docs PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python;QT_QPA_PLATFORM=offscreen")
            # Release notes come from CHANGELOG.md; a version without notes can't be released.
            add_test(NAME release-notes COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/release-notes.py v0.8.0)
            set_tests_properties(release-notes PROPERTIES PASS_REGULAR_EXPRESSION "Template library" FAIL_REGULAR_EXPRESSION "<!--")
            add_test(NAME release-notes-missing COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/release-notes.py v0.0.9)
            set_tests_properties(release-notes-missing PROPERTIES WILL_FAIL ON)
        endif()
    endif()
endif()
