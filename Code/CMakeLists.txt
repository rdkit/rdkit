cmake_minimum_required(VERSION 2.6)
project (RDKit)

# pull in boost:
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

enable_testing()

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing")
endif()

include_directories(.)
add_subdirectory(RDGeneral)
add_subdirectory(DataStructs)
add_subdirectory(Geometry)
add_subdirectory(GraphMol)

