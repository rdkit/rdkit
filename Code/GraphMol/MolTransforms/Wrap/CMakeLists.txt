rdkit_python_extension(rdMolTransforms 
                       rdMolTransforms.cpp
                       DEST Chem 
                       LINK_LIBRARIES MolTransforms 
EigenSolvers 
GraphMol DataStructs RDGeometryLib RDGeneral 
RDBoost)

add_pytest(pyMolTransforms ${CMAKE_CURRENT_SOURCE_DIR}/testMolTransforms.py)
