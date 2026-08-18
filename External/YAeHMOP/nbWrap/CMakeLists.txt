
rdkit_nanobind_extension(rdEHTTools
                         rdEHTTools.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         EHTLib )

add_pytest(pyEHTTools
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testEHTTools.py)
