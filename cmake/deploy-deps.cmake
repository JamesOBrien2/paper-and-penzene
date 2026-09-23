# Copies every non-system DLL that the exe and its Qt plugins need next to the exe.
# windeployqt only handles Qt itself; RDKit, Boost, ICU etc. come from here.
# Usage: cmake -DDEST=<dir with penzene.exe> -DSEARCH=<env>/Library/bin -P deploy-deps.cmake
set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "windows+pe")
set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "dumpbin")
file(GLOB_RECURSE plugins "${DEST}/*.dll")
file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES "${DEST}/penzene.exe"
    LIBRARIES ${plugins}
    DIRECTORIES "${SEARCH}"
    RESOLVED_DEPENDENCIES_VAR deps
    UNRESOLVED_DEPENDENCIES_VAR missing
    PRE_EXCLUDE_REGEXES "^api-ms-" "^ext-ms-"
    POST_EXCLUDE_REGEXES "[Ss]ystem32")
foreach(dep ${deps})
    get_filename_component(name "${dep}" NAME)
    if(NOT EXISTS "${DEST}/${name}")
        file(COPY "${dep}" DESTINATION "${DEST}")
    endif()
endforeach()
if(missing)
    message(STATUS "Unresolved (expected to be system DLLs): ${missing}")
endif()
