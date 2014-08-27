rdkit_library(FMCS
              FMCS.cpp Seed.cpp MaximumCommonSubgraph.cpp SubstructMatchCustom.cpp
              LINK_LIBRARIES Depictor FileParsers ChemTransforms SubstructMatch)

rdkit_headers(FMCS.h
              Seed.h
              SeedSet.h
			  MatchTable.h
			  RingMatchTableSet.h
			  SubstructureCache.h
              SubstructMatchCustom.h
              MaximumCommonSubgraph.h
              DEST GraphMol/FMCS)

rdkit_test(testFMCS testFMCS_Unit.cpp LINK_LIBRARIES
FMCS ChemTransforms Depictor FileParsers SmilesParse 
GraphMol RDGeneral RDGeometryLib SubstructMatch  ${RDKit_THREAD_LIBS})

add_subdirectory(Wrap)
add_subdirectory(Test)
