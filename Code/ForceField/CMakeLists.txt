rdkit_library(ForceField
              ForceField.cpp 
              UFF/AngleBend.cpp UFF/BondStretch.cpp UFF/Nonbonded.cpp 
              UFF/TorsionAngle.cpp UFF/DistanceConstraint.cpp UFF/Params.cpp 
              LINK_LIBRARIES Optimizer)

rdkit_test(testForceField testForceField.cpp 
LINK_LIBRARIES ForceField Optimizer RDGeometryLib RDGeneral )

add_subdirectory(Wrap)




