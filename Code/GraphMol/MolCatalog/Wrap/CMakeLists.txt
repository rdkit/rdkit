rdkit_python_extension(rdMolCatalog rdMolCatalog.cpp
                       DEST Chem
                       LINK_LIBRARIES 
MolCatalog FileParsers Catalogs SmilesParse GraphMol DataStructs RDGeneral RDGeometryLib 
RDBoost)

add_test(pyMolCatalog python ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)



