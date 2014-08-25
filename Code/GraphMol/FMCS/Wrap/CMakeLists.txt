rdkit_python_extension(rdFMCS
                       rdFMCS.cpp
                       DEST Chem 
                       LINK_LIBRARIES FMCS
GraphMol DataStructs RDGeometryLib RDGeneral 
RDBoost)

add_pytest(pyFMCS ${CMAKE_CURRENT_SOURCE_DIR}/testFMCS.py)
