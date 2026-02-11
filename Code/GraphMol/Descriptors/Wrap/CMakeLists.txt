remove_definitions(-DRDKIT_DESCRIPTORS_BUILD)

rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES Descriptors Fingerprints)

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)
add_pytest(pyMolDescriptors3D ${CMAKE_CURRENT_SOURCE_DIR}/test3D.py)

if(RDK_BUILD_OSMORDRED)
  add_pytest(pyMolDescriptorsOsmordred ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_basic.py)
endif(RDK_BUILD_OSMORDRED)
