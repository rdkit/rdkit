cmake_minimum_required(VERSION 2.6)
project (RDKit)

#------- 
# Config variables:
set(RDKit_LibDir $ENV{RDBASE}/bin)
# pull in python:
set(RDKit_PythonDir $ENV{RDBASE}/rdkit)


include(RDKitUtils.cmake)

# pull in python:
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
link_directories(${PYTHON_LIBRARIES})

# pull in boost:
find_package(Boost COMPONENTS python thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

enable_testing()

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing")
endif()

include_directories(.)
add_subdirectory(RDGeneral)
add_subdirectory(RDBoost)
add_subdirectory(DataStructs)
add_subdirectory(Geometry)
add_subdirectory(GraphMol)

