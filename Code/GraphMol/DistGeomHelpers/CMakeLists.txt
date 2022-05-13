
rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp EmbedderUtils.cpp
              LINK_LIBRARIES ForceFieldHelpers SubstructMatch GraphMol DistGeometry Alignment
                )
target_compile_definitions(DistGeomHelpers PRIVATE RDKIT_DISTGEOMHELPERS_BUILD)

rdkit_headers(BoundsMatrixBuilder.h
              Embedder.h DEST GraphMol/DistGeomHelpers)

rdkit_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers MolAlign MolTransforms FileParsers SmilesParse  )

rdkit_catch_test(distGeomHelpersCatch catch_tests.cpp 
LINK_LIBRARIES DistGeomHelpers MolAlign MolTransforms FileParsers SmilesParse )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()