cmake_minimum_required(VERSION 3.21)
project(penzene VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Qt6 REQUIRED COMPONENTS Widgets Svg)
find_package(RDKit REQUIRED CONFIG)

# Everything except main.cpp, so tests link the same code as the app.
add_library(penzene_core STATIC
    src/Document.cpp src/Chem.cpp src/Canvas.cpp src/MainWindow.cpp)
target_include_directories(penzene_core PUBLIC src)
target_compile_definitions(penzene_core PUBLIC PENZENE_VERSION="${PROJECT_VERSION}")
target_link_libraries(penzene_core PUBLIC
    Qt6::Widgets Qt6::Svg
    RDKit::SmilesParse RDKit::FileParsers RDKit::Depictor RDKit::GraphMol
    RDKit::Descriptors RDKit::RDInchiLib RDKit::RDChemDrawLib)

set(PENZENE_ICNS resources/penzene.icns)
set_source_files_properties(${PENZENE_ICNS} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
add_executable(penzene WIN32 MACOSX_BUNDLE src/main.cpp ${PENZENE_ICNS} resources/penzene.rc)
qt_add_resources(penzene "icons" PREFIX "/" BASE resources FILES resources/logo.svg)
target_link_libraries(penzene PRIVATE penzene_core)
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
target_link_libraries(penzene_tests PRIVATE penzene_core Qt6::Test Catch2::Catch2WithMain)
add_test(NAME core COMMAND penzene_tests)
add_test(NAME version COMMAND penzene --version)
target_compile_definitions(penzene_tests PRIVATE PENZENE_TEST_DATA="${CMAKE_SOURCE_DIR}/tests/data")
