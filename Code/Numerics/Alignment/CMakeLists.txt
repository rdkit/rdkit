rdkit_library(Alignment AlignPoints.cpp LINK_LIBRARIES RDGeometry)

rdkit_test(alignment testAlignment.cpp LINK_LIBRARIES Alignment)

add_subdirectory(Wrap)
