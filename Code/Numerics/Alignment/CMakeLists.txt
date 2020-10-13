
rdkit_library(Alignment AlignPoints.cpp LINK_LIBRARIES RDGeometryLib)
target_compile_definitions(Alignment PRIVATE RDKIT_ALIGNMENT_BUILD)

rdkit_headers(AlignPoints.h DEST Numerics/Alignment)

rdkit_test(testAlignment testAlignment.cpp LINK_LIBRARIES Alignment )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
