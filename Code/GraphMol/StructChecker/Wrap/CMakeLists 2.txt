remove_definitions(-DRDKIT_STRUCTCHECKER_BUILD)
rdkit_python_extension(rdStructChecker
        structchecker.cpp
        DEST Chem
        LINK_LIBRARIES StructChecker ChemTransforms Depictor Descriptors
                          SubstructMatch FileParsers SmilesParse 
                          GraphMol RDGeneral RDGeometryLib)
add_pytest(pyStructChecker ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)
