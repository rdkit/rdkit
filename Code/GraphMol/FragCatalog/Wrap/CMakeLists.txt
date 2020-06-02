remove_definitions(-DRDKIT_FRAGCATALOG_BUILD)
rdkit_python_extension(rdfragcatalog
                       rdfragcatalog.cpp FragCatalog.cpp FragCatParams.cpp 
                       FragCatGenerator.cpp FragFPGenerator.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
                       FragCatalog FileParsers) 

add_pytest(pyFragCatalog ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
