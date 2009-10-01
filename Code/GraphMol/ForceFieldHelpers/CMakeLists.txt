# moved library build directives to UFF (actual location of ForceFieldHelpers 
# modules.. might have been a bad idea).
add_subdirectory(UFF)

rdkit_test(testForceFieldHelpers testHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment GraphMol RDGeometryLib RDGeneral )

add_subdirectory(Wrap)

