cmake_minimum_required(VERSION 3.8)
project(biosequences)

add_subdirectory("tests/googletest")

ENABLE_TESTING()
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
OPTION(BUILD_TESTING "Build tests." ON)

# Get the remote url
execute_process(
        COMMAND git remote get-url origin
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_ORIGIN_URL
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

## Include directories
include_directories("tests/googletest/googletest/include/")
include_directories("tests/googletest/googletest/src/")

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES Biosequence.cpp Biosequence.h biosymbols.h alphabets.h)
add_library(biosequences SHARED ${SOURCE_FILES})

IF(BUILD_TESTING)
    add_executable(runBasicTests tests/biosequence_checks.cpp)
    set_target_properties(runBasicTests PROPERTIES
            COMPILE_FLAGS "-g -O0 -ftest-coverage -fprofile-arcs"
            LINK_FLAGS    "--coverage"
            )
    target_link_libraries(runBasicTests gtest gtest_main biosequences)
    add_test(basicTests runBasicTests)
ENDIF()