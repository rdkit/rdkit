rdkit_python_extension(rdPartialCharges
                       rdPartialCharges.cpp
                       DEST Chem 
                       LINK_LIBRARIES PartialCharges 
GraphMol RDGeometryLib RDGeneral DataStructs
RDBoost) 

add_pytest(pyPartialCharges ${CMAKE_CURRENT_SOURCE_DIR}/testPartialCharges.py)
