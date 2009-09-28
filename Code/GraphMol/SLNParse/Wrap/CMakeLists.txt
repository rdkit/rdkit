rdkit_python_extension(rdSLNParse rdSLNParse.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SLNParse SmilesParse SubstructMatch GraphMol
                       RDGeometry RDGeneral DataStructs RDBoost)

add_test(pySLNParse python ${CMAKE_CURRENT_SOURCE_DIR}/testSLN.py)

