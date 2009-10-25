rdkit_python_extension(rdfragcatalog
                       rdfragcatalog.cpp FragCatalog.cpp FragCatParams.cpp 
                       FragCatGenerator.cpp FragFPGenerator.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
                       FragCatalog 
Subgraphs SubstructMatch SmilesParse Catalogs FileParsers GraphMol DataStructs RDGeometryLib RDGeneral 
${LAPACK_LIBRARIES} RDBoost) 

add_test(pyFragCatalog python ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
