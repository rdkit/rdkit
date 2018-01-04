rdkit_library(MolInterchange
              MolInterchange.cpp
              LINK_LIBRARIES GraphMol)

rdkit_headers(MolInterchange.h
              DEST GraphMol/MolInterchange)

rdkit_test(molInterchangeTest1 test1.cpp
           LINK_LIBRARIES MolInterchange FileParsers SmilesParse GraphMol RDGeneral )
