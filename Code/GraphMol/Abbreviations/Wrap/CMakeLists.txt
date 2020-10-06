rdkit_python_extension(rdAbbreviations
                       rdAbbreviations.cpp
                       DEST Chem
                       LINK_LIBRARIES Abbreviations 
GraphMol ) 

add_pytest(pyAbbreviations ${CMAKE_CURRENT_SOURCE_DIR}/testAbbreviations.py)

