rdkit_library(ForceFieldHelpers AtomTyper.cpp Builder.cpp
              LINK_LIBRARIES SmilesParse SubstructMatch ForceField)

rdkit_headers(AtomTyper.h
              Builder.h DEST GraphMol/ForceFieldHelpers/UFF)

