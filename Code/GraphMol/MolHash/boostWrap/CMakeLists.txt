remove_definitions(-DRDKIT_MOLHASH_BUILD)
rdkit_python_extension(rdMolHash
                       rdMolHash.cpp
                       DEST Chem 
                       LINK_LIBRARIES MolHash)

add_pytest(pyMolHash ${CMAKE_CURRENT_SOURCE_DIR}/testMolHash.py)
