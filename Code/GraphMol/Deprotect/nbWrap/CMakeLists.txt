remove_definitions(-DRDKIT_SUBSTRUCTLIBRARY_BUILD)
rdkit_nanobind_extension(rdDeprotect
                         rdDeprotect.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         Deprotect)

add_pytest(pyDeprotect ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/rough_test.py)
