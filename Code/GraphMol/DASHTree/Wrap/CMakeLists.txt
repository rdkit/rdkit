remove_definitions(-DRDKIT_DASHTREE_BUILD)

rdkit_python_extension(rdDASHTree
                       rdDASHTree.cpp
                       DEST Chem
                       LINK_LIBRARIES DASHTree)
