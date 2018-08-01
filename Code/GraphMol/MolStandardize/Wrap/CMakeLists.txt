rdkit_python_extension(rdMolStandardize rdMolStandardize.cpp
                       DEST Chem/MolStandardize
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize SmilesParse
SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost )

rdkit_python_extension(Metal Metal.cpp
                       DEST Chem/MolStandardize
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize SmilesParse
SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost )

rdkit_python_extension(Charge Charge.cpp
                       DEST Chem/MolStandardize
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize SmilesParse
SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost ) 

rdkit_python_extension(Fragment Fragment.cpp
                       DEST Chem/MolStandardize
                       LINK_LIBRARIES
LINK_LIBRARIES MolStandardize SmilesParse
SubstructMatch GraphMol RDGeneral RDGeometryLib RDBoost ) 

add_pytest(pyMolStandardize ${CMAKE_CURRENT_SOURCE_DIR}/testMolStandardize.py)
