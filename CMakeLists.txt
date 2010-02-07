cmake_minimum_required(VERSION 2.6)
project (RDKit)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

enable_testing()

set(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Code/cmake/Modules/")

option(RDK_BUILD_SWIG_WRAPPERS "build the experimental SWIG wrappers" OFF )
option(RDK_BUILD_PYTHON_WRAPPERS "build the standard python wrappers" ON )

#------- 
# Config variables:
set(RDKit_CodeDir "${CMAKE_SOURCE_DIR}/Code")
set(RDKit_ExternalDir "${CMAKE_SOURCE_DIR}/External")
set(RDKit_BinDir "${CMAKE_SOURCE_DIR}/bin")
set(RDKit_LibDir "${CMAKE_SOURCE_DIR}/lib")
set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)

# extra boost versions
if(MSVC)
SET(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0" "1.40" "1.40.0")
endif(MSVC)

if(RDK_BUILD_PYTHON_WRAPPERS)
  #------- 
  # pull in python:
  find_package(PythonLibs)
  include_directories(${PYTHON_INCLUDE_DIR})
  link_directories(${PYTHON_LIBRARIES})
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_NUMPY_INCLUDE_PATH})
  find_package(Boost 1.39.0 COMPONENTS python REQUIRED)
else(RDK_BUILD_PYTHON_WRAPPERS)
  find_package(Boost 1.39.0 REQUIRED)
endif(RDK_BUILD_PYTHON_WRAPPERS)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${RDKit_CodeDir})

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing -fPIC")
endif()


add_subdirectory(External)
add_subdirectory(Code)
add_subdirectory(Projects)
add_subdirectory(rdkit)
