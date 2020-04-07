
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_DISTGEOMHELPERS_BUILD)
rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp
              LINK_LIBRARIES ForceFieldHelpers GraphMol DistGeometry Alignment
                )

rdkit_headers(BoundsMatrixBuilder.h
              Embedder.h DEST GraphMol/DistGeomHelpers)

rdkit_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers MolAlign MolTransforms FileParsers SmilesParse SubstructMatch )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()