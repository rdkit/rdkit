remove_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_python_extension(rdSubstructLibrary
                       rdSubstructLibrary.cpp
                       SubstructLibraryWrap.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SubstructLibrary )

add_pytest(pySubstructLibrary ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
