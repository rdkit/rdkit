rdkit_python_extension(rdMolAlign rdMolAlign.cpp
                       DEST Chem
                       LINK_LIBRARIES
MolAlign MolTransforms
EigenSolvers Alignment Optimizer SubstructMatch
GraphMol RDGeneral RDGeometryLib DataStructs
                       RDBoost)

add_pytest(pyMolAlign ${CMAKE_CURRENT_SOURCE_DIR}/testMolAlign.py)





