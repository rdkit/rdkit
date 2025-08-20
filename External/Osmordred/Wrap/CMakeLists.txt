remove_definitions(-DRDKIT_OSMORDRED_BUILD)
rdkit_python_extension(rdOsmordred
                       rdOsmordred.cpp
                       DEST Chem
                       LINK_LIBRARIES Osmordred)

# install the osomordred python files


file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/Osmordred" DESTINATION "${RDKit_PythonDir}/Chem")


endif()
# add_pytest(pyOsmordred ${CMAKE_CURRENT_SOURCE_DIR}/testOsmordred.py)

