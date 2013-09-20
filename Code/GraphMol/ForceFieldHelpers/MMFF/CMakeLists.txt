rdkit_test(testMMFFForceFieldHelpers testHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral )
