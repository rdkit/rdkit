
rdkit_library(DistGeomHelpers BoundsMatrixBuilder.cpp Embedder.cpp EmbedderUtils.cpp
              LINK_LIBRARIES MolAlign ForceFieldHelpers SubstructMatch GraphMol DistGeometry Alignment
                )
target_compile_definitions(DistGeomHelpers PRIVATE RDKIT_DISTGEOMHELPERS_BUILD)

rdkit_headers(BoundsMatrixBuilder.h
              Embedder.h DEST GraphMol/DistGeomHelpers)

rdkit_catch_test(testDistGeomHelpers testDgeomHelpers.cpp
           LINK_LIBRARIES
           DistGeomHelpers MolTransforms FileParsers SmilesParse  )

rdkit_catch_test(distGeomHelpersCatch catch_tests.cpp 
LINK_LIBRARIES DistGeomHelpers MolTransforms FileParsers SmilesParse )

rdkit_catch_test(testDgeomHelpersGH3019_0 testGithub3019.cpp
           LINK_LIBRARIES DistGeomHelpers SmilesParse
           CLIARGS --shard-count=2 --shard-index=0 )

rdkit_catch_test(testDgeomHelpersGH3019_1 testGithub3019.cpp
           LINK_LIBRARIES DistGeomHelpers SmilesParse 
           CLIARGS --shard-count=2 --shard-index=1 )


if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()