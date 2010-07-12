add_library(InfoTheory STATIC InfoBitRanker.cpp)

rdkit_python_extension(cEntropy cEntropy.cpp 
                       DEST ML/InfoTheory
                       LINK_LIBRARIES
                       InfoTheory RDGeneral RDBoost)

add_subdirectory(Wrap)
