
rdkit_library(ForceField
              ForceField.cpp 
              UFF/AngleBend.cpp UFF/BondStretch.cpp UFF/Nonbonded.cpp
              UFF/Inversion.cpp UFF/TorsionAngle.cpp
              UFF/DistanceConstraint.cpp UFF/AngleConstraint.cpp
              UFF/TorsionConstraint.cpp UFF/PositionConstraint.cpp
              UFF/Params.cpp 
              MMFF/AngleBend.cpp MMFF/BondStretch.cpp MMFF/StretchBend.cpp
              MMFF/OopBend.cpp MMFF/TorsionAngle.cpp MMFF/Nonbonded.cpp
              MMFF/DistanceConstraint.cpp MMFF/AngleConstraint.cpp
              MMFF/TorsionConstraint.cpp MMFF/PositionConstraint.cpp
              MMFF/Params.cpp
              LINK_LIBRARIES Optimizer Trajectory)
target_compile_definitions(ForceField PRIVATE RDKIT_FORCEFIELD_BUILD)

if(RDK_BUILD_ANI)
rdkit_library(ANIForceField
             ANI/AtomicContrib.cpp
             LINK_LIBRARIES Optimizer Trajectory ANIDescriptors EigenSerializer ForceField)
target_compile_definitions(ANIForceField PRIVATE RDKIT_ANIFORCEFIELD_BUILD)
endif(RDK_BUILD_ANI)

rdkit_headers(Contrib.h
              ForceField.h DEST ForceField)

rdkit_headers(UFF/AngleBend.h
              UFF/BondStretch.h
              UFF/Contribs.h
              UFF/DistanceConstraint.h
              UFF/AngleConstraint.h
              UFF/TorsionConstraint.h
              UFF/PositionConstraint.h
              UFF/Nonbonded.h
              UFF/Params.h
              UFF/Inversion.h
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
if(RDK_BUILD_ANI)
rdkit_headers(ANI/AtomicContrib.h DEST ForceField/ANI)
endif(RDK_BUILD_ANI)

add_subdirectory(UFF)
add_subdirectory(MMFF)
if(RDK_BUILD_ANI)
add_subdirectory(ANI)
endif(RDK_BUILD_ANI)
if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()