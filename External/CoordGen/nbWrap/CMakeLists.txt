rdkit_nanobind_extension(rdCoordGen
                         rdCoordGen.cpp
                         DEST Chem
                         LINK_LIBRARIES
                         coordgen MolAlign SubstructMatch GraphMol DataStructs RDGeometryLib RDGeneral)
add_pytest(pyCoordGen
         ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testCoordGen.py)
