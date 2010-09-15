rdkit_library(PartialCharges
              GasteigerCharges.cpp GasteigerParams.cpp
              LINK_LIBRARIES GraphMol)

rdkit_headers(GasteigerCharges.h
              GasteigerParams.h DEST GraphMol/PartialCharges)

add_subdirectory(Wrap)
