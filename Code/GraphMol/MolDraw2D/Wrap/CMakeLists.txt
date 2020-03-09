remove_definitions(-DRDKIT_MOLDRAW2D_BUILD)
rdkit_python_extension(rdMolDraw2D
                       rdMolDraw2D.cpp
                       DEST Chem/Draw
                       LINK_LIBRARIES MolDraw2D
Depictor ChemReactions SmilesParse FileParsers SubstructMatch MolTransforms GraphMol
DataStructs EigenSolvers RDGeometryLib RDGeneral RDBoost ${EXTRA_LOCAL_LIBS})

add_pytest(pyMolDraw2D ${CMAKE_CURRENT_SOURCE_DIR}/testMolDraw2D.py)
