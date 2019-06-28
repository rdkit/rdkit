
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLTRANSFORMS_BUILD)
rdkit_library(MolTransforms MolTransforms.cpp
              LINK_LIBRARIES GraphMol EigenSolvers)

rdkit_headers(MolTransforms.h DEST GraphMol/MolTransforms)

rdkit_test(testMolTransforms test1.cpp
           LINK_LIBRARIES MolTransforms
	   SubstructMatch
           FileParsers SmilesParse
	   EigenSolvers GraphMol RDGeometryLib RDGeneral )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
