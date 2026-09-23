cmake_minimum_required(VERSION 3.21)
project(penzene VERSION 0.3.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Qt6 REQUIRED COMPONENTS Widgets Svg)
find_package(RDKit REQUIRED CONFIG)

# The engine: document, chemistry, editing and rendering. QtGui only, no
# widgets, so it can back a command line or Python module as well as the app.
add_library(penzene_core STATIC src/Document.cpp src/Chem.cpp src/Edit.cpp src/Render.cpp)
target_include_directories(penzene_core PUBLIC src)
target_compile_definitions(penzene_core PUBLIC PENZENE_VERSION="${PROJECT_VERSION}")
target_link_libraries(penzene_core PUBLIC
    Qt6::Gui Qt6::Svg
    RDKit::SmilesParse RDKit::FileParsers RDKit::Depictor RDKit::GraphMol
    RDKit::Descriptors RDKit::RDInchiLib)

# The app's widgets; tests link this too, so they exercise the same code.
add_library(penzene_ui STATIC src/Canvas.cpp src/MainWindow.cpp)
target_link_libraries(penzene_ui PUBLIC penzene_core Qt6::Widgets)

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
find_package(Qt6 REQUIRED COMPONENTS Test)
enable_testing()
add_executable(penzene_tests tests/test_core.cpp tests/test_canvas.cpp)
target_link_libraries(penzene_tests PRIVATE penzene_ui Qt6::Test Catch2::Catch2WithMain)
add_test(NAME core COMMAND penzene_tests)
add_test(NAME version COMMAND penzene --version)
# Batch CLI: a MOL file and a SMILES into a directory, RSC style, freshly laid out.
add_test(NAME render COMMAND penzene --render ${CMAKE_SOURCE_DIR}/tests/data/caffeine.mol CCO
         --out ${CMAKE_BINARY_DIR}/render-test --format png --drawing-style RSC --clean)
set_tests_properties(render PROPERTIES ENVIRONMENT QT_QPA_PLATFORM=offscreen
                     PASS_REGULAR_EXPRESSION "caffeine\\.png.*structure\\.png")
target_compile_definitions(penzene_tests PRIVATE PENZENE_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
