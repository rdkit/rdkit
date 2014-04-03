rdkit_library(FMCS
              FMCS.cpp Seed.cpp MaximumCommonSubgraph.cpp SubstructMatchCustom.cpp
              LINK_LIBRARIES Depictor FileParsers ChemTransforms)

rdkit_headers(FMCS.h
              Seed.h
              SeedSet.h
              SubstructMatchCustom.h
              MaximumCommonSubgraph.h
              DEST GraphMol/FMCS)

rdkit_test(testFMCS testFMCS.cpp LINK_LIBRARIES
FMCS ChemTransforms Depictor FileParsers SmilesParse 
GraphMol RDGeneral RDGeometryLib )
