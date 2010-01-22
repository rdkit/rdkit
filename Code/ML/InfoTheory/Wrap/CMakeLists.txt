rdkit_python_extension(rdInfoTheory 
                       InfoBitRanker.cpp BitCorrMatGenerator.cpp 
                       rdInfoTheory.cpp 
                       DEST ML/InfoTheory
                       LINK_LIBRARIES
                       InfoTheory RDGeneral DataStructs RDBoost)

add_pytest(pyRanker ${CMAKE_CURRENT_SOURCE_DIR}/testRanker.py)
