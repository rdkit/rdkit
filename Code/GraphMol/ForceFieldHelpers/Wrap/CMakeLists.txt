rdkit_python_extension(rdForceFieldHelpers rdForceFields.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ForceFieldHelpers 
SmilesParse SubstructMatch GraphMol ForceField Optimizer RDGeometryLib RDGeneral DataStructs
RDBoost) 

add_pytest(pyForceFieldHelpers ${CMAKE_CURRENT_SOURCE_DIR}/testHelpers.py)




