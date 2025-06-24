
rdkit_library(ForceField
              ForceField.cpp AngleConstraint.cpp AngleConstraints.cpp 
              DistanceConstraint.cpp DistanceConstraints.cpp
              PositionConstraint.cpp TorsionConstraint.cpp
              UFF/AngleBend.cpp UFF/BondStretch.cpp UFF/Nonbonded.cpp
              UFF/Inversion.cpp UFF/Inversions.cpp UFF/TorsionAngle.cpp
              UFF/Params.cpp UFF/Utils.cpp
              MMFF/AngleBend.cpp MMFF/BondStretch.cpp MMFF/StretchBend.cpp
              MMFF/OopBend.cpp MMFF/TorsionAngle.cpp MMFF/Nonbonded.cpp
              MMFF/Params.cpp
              LINK_LIBRARIES Optimizer Trajectory)
target_compile_definitions(ForceField PRIVATE RDKIT_FORCEFIELD_BUILD)

rdkit_headers(Contrib.h
              ForceField.h
              AngleConstraint.h
              AngleConstraints.h
              DistanceConstraint.h
              DistanceConstraints.h
              PositionConstraint.h
              TorsionConstraint.h
              DEST ForceField)

rdkit_headers(UFF/AngleBend.h
              UFF/BondStretch.h
              UFF/Contribs.h
              UFF/DistanceConstraint.h
              UFF/AngleConstraint.h
              UFF/TorsionConstraint.h
              UFF/PositionConstraint.h
              UFF/Nonbonded.h
              UFF/Params.h
              UFF/Utils.h
              UFF/Inversion.h
              UFF/Inversions.h
              UFF/TorsionAngle.h DEST ForceField/UFF)

rdkit_headers(MMFF/AngleBend.h
              MMFF/BondStretch.h
              MMFF/StretchBend.h
              MMFF/OopBend.h
              MMFF/TorsionAngle.h
              MMFF/Nonbonded.h
              MMFF/Contribs.h
              MMFF/DistanceConstraint.h
              MMFF/AngleConstraint.h
              MMFF/TorsionConstraint.h
              MMFF/PositionConstraint.h
              MMFF/Params.h DEST ForceField/MMFF)

rdkit_catch_test(ForceFieldsCatch 
                 catch_tests.cpp
                 LINK_LIBRARIES 
                 ForceField ForceFieldHelpers)

add_subdirectory(UFF)
add_subdirectory(MMFF)
if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()