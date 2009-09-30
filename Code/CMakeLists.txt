cmake_minimum_required(VERSION 2.6)
project (RDKit)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#------- 
# Config variables:
set(RDKit_LibDir $ENV{RDBASE}/bin)
set(RDKit_PythonDir $ENV{RDBASE}/rdkit)

include(RDKitUtils.cmake)

# pull in python:
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
link_directories(${PYTHON_LIBRARIES})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_PATH})

# pull in boost:
if(MSVC)
# I cannot, for the life of me, get this working without hard coding BOOST_ROOT
# here
set(BOOST_ROOT "c:\\boost\\boost_1_40_0")
endif(MSVC)
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
add_subdirectory(Geometry)
add_subdirectory(GraphMol)
add_subdirectory(Query)

