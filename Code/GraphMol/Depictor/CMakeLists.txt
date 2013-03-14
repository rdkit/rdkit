rdkit_library(Depictor
              RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp
              LINK_LIBRARIES GraphMol)

rdkit_headers(DepictUtils.h
              EmbeddedFrag.h
              RDDepictor.h DEST GraphMol/Depictor)

rdkit_test(testDepictor testDepictor.cpp 
           LINK_LIBRARIES Depictor ChemTransforms FileParsers SmilesParse SubstructMatch GraphMol RDGeometryLib RDGeneral )

FILE(GLOB TEST_OUTPUT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/test_data/*.out.sdf")
get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
        "${extra_clean_files};${TEST_OUTPUT_FILES}")


add_subdirectory(Wrap)



