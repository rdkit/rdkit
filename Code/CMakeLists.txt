
#------- 
# Config variables:
set(RDKit_LibDir $ENV{RDBASE}/bin)
set(RDKit_PythonDir $ENV{RDBASE}/rdkit)

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)

# pull in python:
# start with a bit of hackery to allow the user to provide their own
# path to python:
if(PYTHON_LIBRARIES)
  set(oPYTHON_LIBRARIES ${PYTHON_LIBRARIES})
endif(PYTHON_LIBRARIES)
if(PYTHON_INCLUDE_DIRS)
  set(oPYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
endif(PYTHON_INCLUDE_DIRS)
find_package(PythonLibs)
if(oPYTHON_LIBRARIES)
  set(PYTHON_LIBRARIES ${oPYTHON_LIBRARIES})
endif(oPYTHON_LIBRARIES)
if(oPYTHON_INCLUDE_DIRS)
  set(PYTHON_INCLUDE_DIRS ${oPYTHON_INCLUDE_DIRS})
endif(oPYTHON_INCLUDE_DIRS)


if(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIRS)
  set(PYTHON_FOUND "NO")
else(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIRS)
  set(PYTHON_FOUND "YES")
endif(NOT PYTHON_LIBRARIES AND NOT PYTHON_INCLUDE_DIRS)
if(PYTHON_FOUND)
  MESSAGE(STATUS "Found Python libraries in ${PYTHON_INCLUDE_DIRS} as ${PYTHON_LIBRARIES}")
else(PYTHO_NFOUND)
  MESSAGE(FATAL_ERROR "Python libraries not found")
endif(PYTHON_FOUND)


include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_PATH})

# pull in boost:
if(MSVC)
SET(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0" "1.40" "1.40.0")
endif(MSVC)
find_package(Boost COMPONENTS python thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(RDGeneral)
if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(RDBoost)
endif(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Geometry)
add_subdirectory(Numerics)
add_subdirectory(GraphMol)
add_subdirectory(Query)

