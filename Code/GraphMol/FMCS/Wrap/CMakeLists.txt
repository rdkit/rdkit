remove_definitions(-DRDKIT_FMCS_BUILD)
rdkit_python_extension(rdFMCS
                       rdFMCS.cpp
                       DEST Chem 
                       LINK_LIBRARIES FMCS)

add_pytest(pyFMCS ${CMAKE_CURRENT_SOURCE_DIR}/testFMCS.py)
