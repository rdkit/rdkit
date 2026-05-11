remove_definitions(-DRDKIT_SLNPARSE_BUILD)
rdkit_python_extension(rdSLNParse rdSLNParse.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SLNParse )

add_pytest(pySLNParse ${CMAKE_CURRENT_SOURCE_DIR}/testSLN.py)
