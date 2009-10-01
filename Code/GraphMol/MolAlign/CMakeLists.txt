rdkit_library(MolAlign AlignMolecules.cpp
              LINK_LIBRARIES MolTransforms SubstructMatch Alignment)

rdkit_test(testMolAlign testMolAlign.cpp
           LINK_LIBRARIES MolAlign 
FileParsers DistGeomHelpers MolTransforms
SmilesParse DistGeometry ForceFieldHelpers ForceField Optimizer EigenSolvers Alignment SubstructMatch
GraphMol RDGeneral RDGeometryLib )

add_subdirectory(Wrap)




