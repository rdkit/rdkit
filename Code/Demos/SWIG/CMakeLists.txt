find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(RDKFuncs.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(RDKFuncs.i PROPERTIES SWIG_FLAGS "-module RDKFuncs")
set(RDKit_Wrapper_Libs FileParsers_static SmilesParse_static
Depictor_static 
ForceFieldHelpers_static ForceField_static MolAlign_static  Optimizer_static Alignment_static EigenSolvers_static
DistGeomHelpers_static DistGeometry_static
SubstructMatch_static ChemReactions_static Fingerprints_static Descriptors_static
Subgraphs_static GraphMol_static DataStructs_static
RDGeometryLib_static RDGeneral_static )

#add_subdirectory(python_example)
add_subdirectory(java_example)
