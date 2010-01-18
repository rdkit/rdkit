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

#------- 
# Config variables:
set(RDKit_CodeDir "${CMAKE_SOURCE_DIR}/Code")
set(RDKit_ExternalDir "${CMAKE_SOURCE_DIR}/External")
set(RDKit_LibDir "${CMAKE_SOURCE_DIR}/bin")
set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)

#------- 
# pull in python:
# start with a bit of hackery to allow the user to provide their own
# path to python:
if(PYTHON_LIBRARIES)
  set(oPYTHON_LIBRARIES ${PYTHON_LIBRARIES})
endif(PYTHON_LIBRARIES)
if(PYTHON_INCLUDE_DIR)
  set(oPYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR})
endif(PYTHON_INCLUDE_DIR)
find_package(PythonLibs)
if(oPYTHON_LIBRARIES)
  set(PYTHON_LIBRARIES ${oPYTHON_LIBRARIES})
endif(oPYTHON_LIBRARIES)
if(oPYTHON_INCLUDE_DIR)
  set(PYTHON_INCLUDE_DIR ${oPYTHON_INCLUDE_DIR})
endif(oPYTHON_INCLUDE_DIR)


if(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIR)
  set(PYTHON_FOUND "NO")
else(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIR)
  set(PYTHON_FOUND "YES")
endif(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIR)
if(PYTHON_FOUND)
  MESSAGE(STATUS "Found Python libraries in ${PYTHON_INCLUDE_DIR} as ${PYTHON_LIBRARIES}")
else(PYTHON_FOUND)
  MESSAGE(FATAL_ERROR "Python libraries not found")
endif(PYTHON_FOUND)


include_directories(${PYTHON_INCLUDE_DIR})
link_directories(${PYTHON_LIBRARIES})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_PATH})

# pull in boost:
if(MSVC)
SET(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0" "1.40" "1.40.0")
endif(MSVC)
find_package(Boost 1.39.0 COMPONENTS python thread REQUIRED)
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
