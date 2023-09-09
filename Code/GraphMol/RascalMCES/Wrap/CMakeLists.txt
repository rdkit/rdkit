remove_definitions(-DRDKIT_RASCALMCES_BUILD)

rdkit_python_extension(rdRascalMCES
        rdRascalMCES.cpp
        DEST Chem
        LINK_LIBRARIES RascalMCES)

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/testRascalMCES.py)
