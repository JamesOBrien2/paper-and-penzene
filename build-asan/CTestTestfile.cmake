# CMake generated Testfile for 
# Source directory: /Users/user/PhD/Github/penzene-work
# Build directory: /Users/user/PhD/Github/penzene-work/build-asan
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test("core" "/Users/user/PhD/Github/penzene-work/build-asan/bin/penzene_tests")
set_tests_properties("core" PROPERTIES  _BACKTRACE_TRIPLES "/Users/user/PhD/Github/penzene-work/CMakeLists.txt;54;add_test;/Users/user/PhD/Github/penzene-work/CMakeLists.txt;0;")
add_test("version" "/Users/user/PhD/Github/penzene-work/build-asan/bin/penzene.app/Contents/MacOS/penzene" "--version")
set_tests_properties("version" PROPERTIES  _BACKTRACE_TRIPLES "/Users/user/PhD/Github/penzene-work/CMakeLists.txt;55;add_test;/Users/user/PhD/Github/penzene-work/CMakeLists.txt;0;")
add_test("help" "/Users/user/PhD/Github/penzene-work/build-asan/bin/penzene.app/Contents/MacOS/penzene" "--help")
set_tests_properties("help" PROPERTIES  ENVIRONMENT "QT_QPA_PLATFORM=offscreen" PASS_REGULAR_EXPRESSION "usage: penzene" TIMEOUT "20" _BACKTRACE_TRIPLES "/Users/user/PhD/Github/penzene-work/CMakeLists.txt;57;add_test;/Users/user/PhD/Github/penzene-work/CMakeLists.txt;0;")
add_test("render" "/Users/user/PhD/Github/penzene-work/build-asan/bin/penzene.app/Contents/MacOS/penzene" "--render" "/Users/user/PhD/Github/penzene-work/tests/data/aspirin.mol" "CCO" "--out" "/Users/user/PhD/Github/penzene-work/build-asan/render-test" "--format" "png" "--drawing-style" "RSC" "--clean")
set_tests_properties("render" PROPERTIES  ENVIRONMENT "QT_QPA_PLATFORM=offscreen" PASS_REGULAR_EXPRESSION [[aspirin\.png.*structure\.png]] _BACKTRACE_TRIPLES "/Users/user/PhD/Github/penzene-work/CMakeLists.txt;61;add_test;/Users/user/PhD/Github/penzene-work/CMakeLists.txt;0;")
