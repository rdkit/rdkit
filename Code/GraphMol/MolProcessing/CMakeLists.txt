rdkit_library(MolProcessing
              MolProcessing.cpp
             LINK_LIBRARIES FileParsers Fingerprints RDStreams RDGeneral)
target_compile_definitions(MolProcessing PRIVATE RDKIT_MOLPROCESSING_BUILD)

rdkit_headers(MolProcessing.h DEST GraphMol/MolProcessing)

rdkit_catch_test(testMolProcessing processing_tests.cpp 
  LINK_LIBRARIES MolProcessing )

if(RDK_BUILD_PYTHON_WRAPPERS)
 add_subdirectory(Wrap)
endif()
