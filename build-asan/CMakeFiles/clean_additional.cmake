# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles/penzene_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/penzene_autogen.dir/ParseCache.txt"
  "CMakeFiles/penzene_core_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/penzene_core_autogen.dir/ParseCache.txt"
  "CMakeFiles/penzene_tests_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/penzene_tests_autogen.dir/ParseCache.txt"
  "CMakeFiles/penzene_ui_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/penzene_ui_autogen.dir/ParseCache.txt"
  "penzene_autogen"
  "penzene_core_autogen"
  "penzene_tests_autogen"
  "penzene_ui_autogen"
  )
endif()
