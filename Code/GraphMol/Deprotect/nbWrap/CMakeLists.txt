remove_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_python_extension(rdDeprotect
                       rdDeprotect.cpp
                       DeprotectWrap.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       Deprotect )

add_pytest(pyDeprotect ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
