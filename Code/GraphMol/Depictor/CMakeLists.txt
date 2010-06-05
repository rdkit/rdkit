rdkit_library(Depictor
              RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp
              LINK_LIBRARIES GraphMol)

rdkit_test(testDepictor testDepictor.cpp 
           LINK_LIBRARIES Depictor FileParsers SmilesParse SubstructMatch GraphMol RDGeometryLib RDGeneral )

add_subdirectory(Wrap)



