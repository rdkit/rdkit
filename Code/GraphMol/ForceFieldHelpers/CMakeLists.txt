# moved back library build directives here

remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_FORCEFIELDHELPERS_BUILD)
rdkit_library(ForceFieldHelpers UFF/AtomTyper.cpp UFF/Builder.cpp
              MMFF/AtomTyper.cpp MMFF/Builder.cpp CrystalFF/TorsionAngleM6.cpp
              CrystalFF/TorsionPreferences.cpp
              LINK_LIBRARIES SmilesParse SubstructMatch ForceField)
rdkit_headers(UFF/AtomTyper.h
              UFF/Builder.h UFF/UFF.h DEST GraphMol/ForceFieldHelpers/UFF)
rdkit_headers(MMFF/AtomTyper.h
              MMFF/Builder.h MMFF/MMFF.h DEST GraphMol/ForceFieldHelpers/MMFF)
rdkit_headers(CrystalFF/TorsionAngleM6.h CrystalFF/TorsionPreferences.h
              DEST GraphMol/ForceFieldHelpers/CrystalFF)

add_subdirectory(MMFF)
add_subdirectory(UFF)
add_subdirectory(CrystalFF)
add_subdirectory(Wrap)
