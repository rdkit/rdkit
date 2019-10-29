rdkit_test(testCrystalForceField testCrystalFF.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral  ${RDKit_THREAD_LIBS})
