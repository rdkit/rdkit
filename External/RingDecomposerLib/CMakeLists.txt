cmake_minimum_required (VERSION 3.1)
project (RingDecomposerLib VERSION 1.0)

enable_testing()

# external options -- build our RDKit testing tool
option(BUILD_RDKIT_BENCHMARK "Build the Benchmark with RDKit" OFF)

# external options -- build python wrapper
option(BUILD_PYTHON_WRAPPER "Build python wrapper" OFF)

# external options -- python executable
set(PYTHON_EXECUTABLE "" CACHE STRING "Python executable for building of the extension")

if(PYTHON_EXECUTABLE STREQUAL "")
  find_program(PYTHON "python")
else()
  set(PYTHON "${PYTHON_EXECUTABLE}")
endif()

# external options -- python flags
set(PYTHON_FLAGS "" CACHE STRING "Python flags")

# external options -- rdkit include directories
set(RDKIT_INCLUDE_DIR "" CACHE STRING "Rdkit include dir")

# external options -- rdkit include directories
set(RDKIT_LIB_DIR "" CACHE STRING "Rdkit librart dir")

set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(LIB_SOURCE_DIR "${SOURCE_DIR}/RingDecomposerLib")
set(PYTHON_SOURCE_DIR "${SOURCE_DIR}/python")

set(SHARED_LIBRARY "RingDecomposerLib")
set(STATIC_LIBRARY "RingDecomposerLibStatic")

set(SHARED_TEST "Test")
set(STATIC_TEST "Test_static")

set(SHARED_MINIMALEXAMPLE "MinimalExample")
set(STATIC_MINIMALEXAMPLE "MinimalExample_static")

set(PYTHON_TEST "${PYTHON_SOURCE_DIR}/test/test.py")
set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")

SET(CMAKE_INSTALL_RPATH "\$ORIGIN/lib")

include_directories(${LIB_SOURCE_DIR})
add_subdirectory("${LIB_SOURCE_DIR}")
add_subdirectory("${SOURCE_DIR}/Test")
add_subdirectory("${SOURCE_DIR}/MinimalExample")
add_subdirectory("${SOURCE_DIR}/RDKitBenchmark")
add_subdirectory("${PYTHON_SOURCE_DIR}")
add_subdirectory("${TEST_DIR}")
