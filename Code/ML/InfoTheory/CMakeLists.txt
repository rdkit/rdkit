add_library(InfoTheory InfoBitRanker.cpp)

rdkit_python_extension(cEntropy cEntropy.cpp 
                       DEST ML/InfoTheory
                       LINK_LIBRARIES
                       InfoTheory RDBoost)

add_subdirectory(Wrap)
