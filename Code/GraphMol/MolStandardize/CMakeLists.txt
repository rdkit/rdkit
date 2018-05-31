rdkit_library(MolStandardize
              MolStandardize.cpp
	      Metal.cpp
	      Normalize.cpp
              LINK_LIBRARIES RDGeneral GraphMol) 

rdkit_headers(MolStandardize.h
	      Metal.h
	      Normalize.h
	      DEST GraphMol/MolStandardize)

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol SubstructMatch)
