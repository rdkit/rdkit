remove_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_nanobind_extension(rdSubstructLibrary
                         rdSubstructLibrary.cpp
                         SubstructLibraryWrap.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         SubstructLibrary)

add_pytest(pySubstructLibrary ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
