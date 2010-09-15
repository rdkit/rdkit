rdkit_library(ForceField
              ForceField.cpp 
              UFF/AngleBend.cpp UFF/BondStretch.cpp UFF/Nonbonded.cpp 
              UFF/TorsionAngle.cpp UFF/DistanceConstraint.cpp UFF/Params.cpp 
              LINK_LIBRARIES Optimizer)

rdkit_headers(Contrib.h
              ForceField.h DEST ForceField)

rdkit_headers(UFF/AngleBend.h
              UFF/BondStretch.h
              UFF/Contribs.h
              UFF/DistanceConstraint.h
              UFF/Nonbonded.h
              UFF/Params.h
              UFF/TorsionAngle.h DEST ForceField/UFF)

rdkit_test(testForceField testForceField.cpp 
LINK_LIBRARIES ForceField Optimizer RDGeometryLib RDGeneral )

add_subdirectory(Wrap)




