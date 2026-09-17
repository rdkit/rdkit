remove_definitions(-DRDKIT_RASCALMCES_BUILD)

rdkit_nanobind_extension(rdRascalMCES
        rdRascalMCES.cpp
        DEST Chem
        LINK_LIBRARIES RascalMCES)

add_pytest(pyRascalMCES ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testRascalMCES.py)
