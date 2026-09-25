rdkit_nanobind_extension(rdinchi
                        pyInchi.cpp
                        DEST Chem
                        LINK_LIBRARIES
                        RDInchiLib)

add_pytest(pyInchiWrap ${CMAKE_SOURCE_DIR}/rdkit/Chem/UnitTestInchi.py)
