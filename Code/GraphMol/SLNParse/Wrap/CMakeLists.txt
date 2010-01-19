rdkit_python_extension(rdSLNParse rdSLNParse.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SLNParse SmilesParse SubstructMatch GraphMol
                       RDGeometryLib RDGeneral DataStructs RDBoost)

add_test(pySLNParse ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testSLN.py)

