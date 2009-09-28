rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp
              LINK_LIBRARIES ForceFieldHelpers DistGeometry Alignment)

rdkit_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers FileParsers MolAlign)

add_subdirectory(Wrap)


