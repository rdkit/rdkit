rdkit_library(Alignment AlignPoints.cpp LINK_LIBRARIES RDGeometryLib)

rdkit_test(testAlignment testAlignment.cpp LINK_LIBRARIES Alignment RDGeometryLib RDGeneral )

add_subdirectory(Wrap)
