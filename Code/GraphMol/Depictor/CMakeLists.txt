rdkit_library(Depictor
              RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp
              LINK_LIBRARIES GraphMol)

rdkit_test(testDepictor testDepictor.cpp 
           LINK_LIBRARIES Depictor SmilesParse FileParsers)

add_subdirectory(Wrap)



