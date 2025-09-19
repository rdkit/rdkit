remove_definitions(-DRDKIT_OSMORDRED_BUILD)
rdkit_python_extension(rdOsmordred
                       rdOsmordred.cpp
                       DEST Chem
                       LINK_LIBRARIES Osmordred)


add_pytest(pyOsmordredBasic ${CMAKE_CURRENT_SOURCE_DIR}/test_osmordred_basic.py)

