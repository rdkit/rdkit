
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_MOLALIGN_BUILD)
rdkit_library(MolAlign AlignMolecules.cpp O3AAlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment RDGeometryLib RDGeneral ForceFieldHelpers)

rdkit_headers(AlignMolecules.h O3AAlignMolecules.h DEST GraphMol/MolAlign)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign Descriptors
FileParsers DistGeomHelpers ForceFieldHelpers
SmilesParse DistGeometry ForceField Optimizer SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})

add_subdirectory(Wrap)
