rdkit_library(MolAlign AlignMolecules.cpp O3AAlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment)

rdkit_headers(AlignMolecules.h O3AAlignMolecules.h DEST GraphMol/MolAlign)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign Descriptors
FileParsers DistGeomHelpers MolTransforms ForceFieldHelpers 
SmilesParse DistGeometry ForceField Optimizer EigenSolvers Alignment SubstructMatch
GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})

add_subdirectory(Wrap)




