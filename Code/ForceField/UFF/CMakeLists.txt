rdkit_test(testUFFForceField testUFFForceField.cpp 
LINK_LIBRARIES ForceField Optimizer RDGeometryLib RDGeneral
  EigenSolvers SmilesParse SubstructMatch GraphMol
  FileParsers MolTransforms ForceFieldHelpers)
