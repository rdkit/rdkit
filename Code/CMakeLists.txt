set(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

#------- 
# Config variables:
set(RDKit_LibDir "${CMAKE_SOURCE_DIR}/bin")
set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)

# pull in python:
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
link_directories(${PYTHON_LIBRARIES})
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_PATH})

# pull in boost:
find_package(Boost COMPONENTS python thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# pull BLAS and LAPACK (Fortran required by LAPACK)
enable_language(Fortran)
find_package(BLAS)
find_package(LAPACK)

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(RDGeneral)
add_subdirectory(RDBoost)
add_subdirectory(DataStructs)
add_subdirectory(Geometry)

add_subdirectory(Numerics)
add_subdirectory(ForceField)
add_subdirectory(DistGeom)
add_subdirectory(Catalogs)

add_subdirectory(GraphMol)
add_subdirectory(Query)

add_subdirectory(DataManip)
add_subdirectory(SimDivPickers)

add_subdirectory(ML)
add_subdirectory(ChemicalFeatures)
