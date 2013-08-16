rdkit_python_extension(rdUFFForceFieldHelpers rdUFFForceFields.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ForceFieldHelpers 
SmilesParse SubstructMatch GraphMol ForceField Optimizer RDGeometryLib RDGeneral DataStructs
RDBoost) 

add_pytest(pyUFFForceFieldHelpers ${CMAKE_CURRENT_SOURCE_DIR}/testUFFHelpers.py)




