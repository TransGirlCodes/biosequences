cmake_minimum_required(VERSION 3.8)
project(biosequences)

add_subdirectory("tests/googletest")

ENABLE_TESTING()
OPTION(BUILD_TESTING "Build tests." ON)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES Biosequence.cpp Biosequence.h biosymbols.h alphabets.h)
add_library(biosequences SHARED ${SOURCE_FILES})

IF(BUILD_TESTING)
    add_executable(runBioseqTests tests/biosequence_checks.cpp)
    set_target_properties(runBioseqTests PROPERTIES
            COMPILE_FLAGS "-g -O0 -ftest-coverage -fprofile-arcs"
            LINK_FLAGS "--coverage")
    target_link_libraries(runBioseqTests biosequences gtest gtest_main)
    add_test(sequenceTests runBioseqTests)
ENDIF()