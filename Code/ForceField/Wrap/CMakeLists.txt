remove_definitions(-DRDKIT_FORCEFIELD_BUILD)
if(RDK_BUILD_ANI)
rdkit_python_extension(rdForceField ForceField.cpp
                       DEST ForceField
                       LINK_LIBRARIES ForceFieldHelpers ANIForceField ANIForceFieldHelpers)
add_pytest(pyForceFieldANI ${CMAKE_CURRENT_SOURCE_DIR}/testANI.py)
else()
rdkit_python_extension(rdForceField ForceField.cpp
                       DEST ForceField
                       LINK_LIBRARIES ForceFieldHelpers)
endif(RDK_BUILD_ANI)
add_pytest(pyForceFieldConstraints ${CMAKE_CURRENT_SOURCE_DIR}/testConstraints.py)
