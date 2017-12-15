rdkit_python_extension(rdCoordGen
                       rdCoordGen.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ${RDK_COORDGEN_LIBS} MolAlign SubstructMatch GraphMol DataStructs RDGeometryLib RDGeneral RDBoost)
add_pytest(pyCoordGen
         ${CMAKE_CURRENT_SOURCE_DIR}/testCoordGen.py)
