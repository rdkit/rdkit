rdkit_python_extension(rdDepictor
                       rdDepictor.cpp
                       DEST Chem 
                       LINK_LIBRARIES 
           Depictor SmilesParse FileParsers SubstructMatch GraphMol RDGeometryLib RDGeneral RDBoost )

add_test(pyDepictor ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/testDepictor.py)
