rdkit_library(MolAlign AlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment)

rdkit_headers(AlignMolecules.h DEST GraphMol/MolAlign)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign 
FileParsers DistGeomHelpers MolTransforms ForceFieldHelpers 
SmilesParse DistGeometry ForceField Optimizer EigenSolvers Alignment SubstructMatch
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)




