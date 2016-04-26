rdkit_library(StructChecker
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp 
              LINK_LIBRARIES Depictor FileParsers ChemTransforms)

rdkit_headers(StructChecker.h
              DEST GraphMol/AVALON)

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker ChemTransforms Depictor FileParsers SmilesParse 
                           GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
