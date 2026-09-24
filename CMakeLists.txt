cmake_minimum_required(VERSION 3.21)
project(penzene VERSION 0.5.0 LANGUAGES CXX)
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

find_package(Qt6 REQUIRED COMPONENTS Gui Svg)
find_package(RDKit REQUIRED CONFIG)

# The engine: document, chemistry, editing and rendering. QtGui only, no
# widgets, so it can back a command line or Python module as well as the app.
add_library(penzene_core STATIC src/Document.cpp src/Chem.cpp src/Edit.cpp src/Render.cpp)
target_include_directories(penzene_core PUBLIC src)
target_compile_definitions(penzene_core PUBLIC PENZENE_VERSION="${PROJECT_VERSION}")
target_link_libraries(penzene_core PUBLIC
    Qt6::Gui Qt6::Svg
    RDKit::SmilesParse RDKit::FileParsers RDKit::Depictor RDKit::GraphMol
    RDKit::Descriptors RDKit::RDInchiLib RDKit::CIPLabeler)

if(PENZENE_APP)
find_package(Qt6 REQUIRED COMPONENTS Widgets Network Test)
# The app's widgets; tests link this too, so they exercise the same code.
add_library(penzene_ui STATIC src/Canvas.cpp src/MainWindow.cpp src/PubChem.cpp)
target_link_libraries(penzene_ui PUBLIC penzene_core Qt6::Widgets Qt6::Network)

set(PENZENE_ICNS resources/penzene.icns)
set_source_files_properties(${PENZENE_ICNS} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
add_executable(penzene WIN32 MACOSX_BUNDLE src/main.cpp ${PENZENE_ICNS} resources/penzene.rc)
qt_add_resources(penzene "icons" PREFIX "/" BASE resources FILES resources/logo.svg)
target_link_libraries(penzene PRIVATE penzene_ui)
set_target_properties(penzene PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER io.github.jamesobrien2.penzene
    MACOSX_BUNDLE_BUNDLE_NAME "Penzene"
    MACOSX_BUNDLE_ICON_FILE penzene.icns
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION})
find_package(Catch2 3 REQUIRED)
enable_testing()
add_executable(penzene_tests tests/test_core.cpp tests/test_canvas.cpp)
target_link_libraries(penzene_tests PRIVATE penzene_ui Qt6::Test Catch2::Catch2WithMain)
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
        if(PENZENE_APP)
            add_test(NAME python COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_python.py)
            set_tests_properties(python PROPERTIES ENVIRONMENT
                "PYTHONPATH=${CMAKE_BINARY_DIR}/python;PENZENE_TEST_DATA=${CMAKE_SOURCE_DIR}/tests/data")
        endif()
    endif()
endif()
