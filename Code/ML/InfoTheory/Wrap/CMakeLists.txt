remove_definitions(-DRDKIT_INFOTHEORY_BUILD)
rdkit_python_extension(rdInfoTheory 
                       InfoBitRanker.cpp BitCorrMatGenerator.cpp 
                       rdInfoTheory.cpp 
                       DEST ML/InfoTheory
                       LINK_LIBRARIES
                       InfoTheory)

add_pytest(pyRanker ${CMAKE_CURRENT_SOURCE_DIR}/testRanker.py)
