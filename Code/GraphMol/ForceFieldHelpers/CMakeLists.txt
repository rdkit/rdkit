# moved back library build directives here
rdkit_library(ForceFieldHelpers UFF/AtomTyper.cpp UFF/Builder.cpp
              MMFF/AtomTyper.cpp MMFF/Builder.cpp
              LINK_LIBRARIES SmilesParse SubstructMatch ForceField)
rdkit_headers(UFF/AtomTyper.h
              UFF/Builder.h DEST GraphMol/ForceFieldHelpers/UFF)
rdkit_headers(MMFF/AtomTyper.h
              MMFF/Builder.h DEST GraphMol/ForceFieldHelpers/MMFF)

add_subdirectory(MMFF)
add_subdirectory(UFF)
add_subdirectory(Wrap)
