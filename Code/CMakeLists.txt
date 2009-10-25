set(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

#------- 
# Config variables:
set(RDKit_ExternalDir "${CMAKE_SOURCE_DIR}/External")
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
if(MSVC)
# I cannot, for the life of me, get this working without hard coding BOOST_ROOT
# here
#set(BOOST_ROOT "c:\\boost\\boost_1_40_0")
SET(Boost_ADDITIONAL_VERSIONS "1.40" "1.40.0")
endif(MSVC)
find_package(Boost COMPONENTS python thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# -----------
# pull BLAS and LAPACK

# Things are easy if we have fortran:
workaround_9220(Fortran W9220_Fortran_COMPILER_WORKS)
IF (W9220_Fortran_COMPILER_WORKS)
enable_language(Fortran)
find_package(BLAS)
find_package(LAPACK)
ELSE()
# otherwise just assume there are lapack and blas libraries out there somewhere:
if(MSVC)
# on windows look for our own copy:
ADD_DEFINITIONS("/D BOOST_NUMERIC_BINDINGS_USE_CLAPACK")
LINK_DIRECTORIES(${RDKit_ExternalDir}/Lapack/win32)
endif(MSVC)
set(LAPACK_LIBRARIES lapack blas)
ENDIF()

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing -fPIC")
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
