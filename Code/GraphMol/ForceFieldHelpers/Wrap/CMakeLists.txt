remove_definitions(-DRDKIT_FORCEFIELDHELPERS_BUILD)
if(RDK_BUILD_ANI)

rdkit_python_extension(rdForceFieldHelpers rdForceFields.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ForceFieldHelpers ANIForceFieldHelpers)
add_pytest(pyForceFieldHelpersANI ${CMAKE_CURRENT_SOURCE_DIR}/testANIHelpers.py)
else()

rdkit_python_extension(rdForceFieldHelpers rdForceFields.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ForceFieldHelpers ) 
endif(RDK_BUILD_ANI)

add_pytest(pyForceFieldHelpers ${CMAKE_CURRENT_SOURCE_DIR}/testHelpers.py)




