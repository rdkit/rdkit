rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp
              LINK_LIBRARIES ForceFieldHelpers DistGeometry Alignment)

rdkit_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse SubstructMatch ForceFieldHelpers ForceField MolAlign Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral )

add_subdirectory(Wrap)


