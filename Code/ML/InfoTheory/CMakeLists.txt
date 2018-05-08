
add_definitions(-DRDKIT_INFOTHEORY_BUILD)
rdkit_library(InfoTheory InfoBitRanker.cpp LINK_LIBRARIES DataStructs RDGeneral)

add_subdirectory(Wrap)
