
rdkit_library(MolTransforms MolTransforms.cpp
              LINK_LIBRARIES GraphMol EigenSolvers RDGeneral )
target_compile_definitions(MolTransforms PRIVATE RDKIT_MOLTRANSFORMS_BUILD)

rdkit_headers(MolTransforms.h DEST GraphMol/MolTransforms)

rdkit_test(testMolTransforms test1.cpp
           LINK_LIBRARIES MolTransforms
           FileParsers SmilesParse )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
