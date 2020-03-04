
add_definitions(-DRDKIT_INFOTHEORY_BUILD)
rdkit_library(InfoTheory InfoBitRanker.cpp LINK_LIBRARIES DataStructs RDGeneral)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
