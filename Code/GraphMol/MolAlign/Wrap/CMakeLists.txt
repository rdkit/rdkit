remove_definitions(-DRDKIT_MOLALIGN_BUILD)
rdkit_python_extension(rdMolAlign rdMolAlign.cpp
                       DEST Chem
                       LINK_LIBRARIES
MolAlign O3AAlign DataStructs Descriptors SmilesParse)

add_pytest(pyMolAlign ${CMAKE_CURRENT_SOURCE_DIR}/testMolAlign.py)





