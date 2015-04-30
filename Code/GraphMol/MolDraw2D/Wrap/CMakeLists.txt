rdkit_python_extension(rdMolDraw2D
                       rdMolDraw2D.cpp
                       DEST Chem/Draw 
                       LINK_LIBRARIES MolDraw2D
GraphMol DataStructs RDGeometryLib RDGeneral 
RDBoost)

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2D.py)
