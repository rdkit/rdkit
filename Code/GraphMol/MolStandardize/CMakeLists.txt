rdkit_library(MolStandardize
              MolStandardize.cpp
	      SHARED
              LINK_LIBRARIES RDGeneral
                 ${RDKit_THREAD_LIBS})

rdkit_headers(MolStandardize.h
	      DEST GraphMol/MolStandardize)

rdkit_test(molStandardizeTest test1.cpp LINK_LIBRARIES MolStandardize
	SmilesParse RDGeneral GraphMol)

