rdkit_python_extension(rdForceField ForceField.cpp
                       DEST ForceField
                       LINK_LIBRARIES ForceField Optimizer RDGeneral RDGeometryLib RDBoost)
add_pytest(pyForceFieldConstraints ${CMAKE_CURRENT_SOURCE_DIR}/testConstraints.py)
