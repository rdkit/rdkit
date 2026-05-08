rdkit_python_extension(rdDetermineBonds
                       rdDetermineBonds.cpp
                       DEST Chem
                       LINK_LIBRARIES DetermineBonds 
GraphMol ) 

add_pytest(pyDetermineBonds ${CMAKE_CURRENT_SOURCE_DIR}/testDetermineBonds.py)

