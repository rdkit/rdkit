
rdkit_library(InfoTheory InfoBitRanker.cpp LINK_LIBRARIES DataStructs RDGeneral)
target_compile_definitions(InfoTheory PRIVATE RDKIT_INFOTHEORY_BUILD)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
