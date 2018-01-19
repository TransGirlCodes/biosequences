cmake_minimum_required(VERSION 3.8)
project(biosequences)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES Biosequence.cpp Biosequence.h biosymbols.h)
add_library(biosequences SHARED ${SOURCE_FILES})

add_executable(biosymbols_check ../test/biosymbols_checks.cpp)
add_executable(bioseq_check ../test/biosequence_checks.cpp)
target_link_libraries(bioseq_check biosequences)