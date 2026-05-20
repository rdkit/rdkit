remove_definitions(-DRDKIT_MOLINTERACTIONFIELDS_BUILD)
rdkit_nanobind_extension(rdMIF
                         rdMIF.cpp
                         DEST Chem
                         LINK_LIBRARIES MolInteractionFields ForceFieldHelpers
                         ForceField PartialCharges SubstructMatch SmilesParse GraphMol RDGeometryLib
                         RDGeneral DataStructs)

add_pytest(pyMIF ${CMAKE_CURRENT_SOURCE_DIR}/../Wrap/testMIF.py)
