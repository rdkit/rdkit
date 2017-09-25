rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp
              LINK_LIBRARIES ForceFieldHelpers DistGeometry Alignment
                ${RDKit_THREAD_LIBS})

rdkit_headers(BoundsMatrixBuilder.h
              Embedder.h DEST GraphMol/DistGeomHelpers)

rdkit_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers DistGeometry ForceFieldHelpers ForceField MolAlign MolTransforms  FileParsers SmilesParse SubstructMatch Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS} )

add_subdirectory(Wrap)


