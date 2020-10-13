
remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
rdkit_library(Trajectory Trajectory.cpp
              LINK_LIBRARIES GraphMol RDGeneral )
target_compile_definitions(Trajectory PRIVATE RDKIT_TRAJECTORY_BUILD)

rdkit_headers(Snapshot.h Trajectory.h DEST GraphMol/Trajectory)

rdkit_test(trajectoryTest trajectoryTest.cpp
           LINK_LIBRARIES Trajectory FileParsers
           SmilesParse SubstructMatch ForceFieldHelpers )
