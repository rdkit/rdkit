rdkit_python_extension(rdMolStandardize rdMolStandardize.cpp
                       DEST Chem
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize SmilesParse
SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost )

add_pytest(pyMolStandardize ${CMAKE_CURRENT_SOURCE_DIR}/testMolStandardize.py)
