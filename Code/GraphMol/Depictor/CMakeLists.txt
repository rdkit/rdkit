rdkit_library(Depictor RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp)

rdkit_test(testDepictor testDepictor.cpp LINK_LIBRARIES
  Depictor SmilesParse FileParsers GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)



