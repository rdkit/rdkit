rdkit_python_extension(rdMIF 
                       rdMIF.cpp
                       DEST Chem
                       LINK_LIBRARIES MolInteractionFields ForceFieldHelpers
                       ForceField PartialCharges SubstructMatch SmilesParse GraphMol RDGeometryLib
                       RDGeneral DataStructs RDBoost)

add_pytest(pyMIF  ${CMAKE_CURRENT_SOURCE_DIR}/testMIF.py)
