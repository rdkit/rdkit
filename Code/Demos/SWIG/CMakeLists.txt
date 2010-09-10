find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(RDKFuncs.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(RDKFuncs.i PROPERTIES SWIG_FLAGS "-module RDKFuncs")
set(RDKit_Wrapper_Libs FileParsers SmilesParse 
Depictor SubstructMatch ChemReactions Fingerprints Subgraphs
GraphMol DataStructs
RDGeometryLib RDGeneral )

#add_subdirectory(python_example)
add_subdirectory(java_example)
