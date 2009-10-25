rdkit_library(PartialCharges
              GasteigerCharges.cpp GasteigerParams.cpp
              LINK_LIBRARIES GraphMol)

add_subdirectory(Wrap)
