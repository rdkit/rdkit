
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_TRAJECTORY_BUILD)
rdkit_library(Trajectory Trajectory.cpp
              LINK_LIBRARIES RDGeneral RDGeometryLib GraphMol)

rdkit_headers(Snapshot.h Trajectory.h DEST GraphMol/Trajectory)

rdkit_test(trajectoryTest trajectoryTest.cpp
           LINK_LIBRARIES Trajectory RDGeneral RDGeometryLib 
           GraphMol FileParsers
           SmilesParse SubstructMatch ForceFieldHelpers ForceField)
