remove_definitions(-DRDKIT_DESCRIPTORS_BUILD)

rdkit_python_extension(rdMolDescriptors rdMolDescriptors.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES Descriptors Fingerprints)

add_pytest(pyMolDescriptors ${CMAKE_CURRENT_SOURCE_DIR}/testMolDescriptors.py)
add_pytest(pyMolDescriptors3D ${CMAKE_CURRENT_SOURCE_DIR}/test3D.py)

if(RDK_BUILD_OSMORDRED_SUPPORT)
  add_pytest(pyMolDescriptorsOsmordred ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_basic.py)
  add_pytest(pyMolDescriptorsOsmordredPorted ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_ported_mordred.py)
  add_pytest(pyMolDescriptorsOsmordredICInvariants ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_ic_invariants.py)
  add_pytest(pyMolDescriptorsOsmordredICBasak ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_ic_basak_polly.py)
endif(RDK_BUILD_OSMORDRED_SUPPORT)
