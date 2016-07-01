rdkit_library(StructChecker
              StructChecker.cpp StructCheckerOptions.cpp StructureFlags.cpp Pattern.cpp AtomSymbolMatch.cpp
              LINK_LIBRARIES Depictor FileParsers ChemTransforms)

rdkit_headers(StructChecker.h Pattern.h
              DEST GraphMol/StructChecker)

rdkit_test(testStructChecker  testStructChecker.cpp 
           LINK_LIBRARIES StructChecker ChemTransforms Depictor FileParsers SmilesParse 
                           GraphMol RDGeneral RDGeometryLib ${RDKit_THREAD_LIBS})
