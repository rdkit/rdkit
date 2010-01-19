rdkit_python_extension(rdMolAlign rdMolAlign.cpp
                       DEST Chem
                       LINK_LIBRARIES
MolAlign MolTransforms
EigenSolvers Alignment Optimizer SubstructMatch
GraphMol RDGeneral RDGeometryLib DataStructs
                       RDBoost)

add_test(pyMolAlign ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testMolAlign.py)





