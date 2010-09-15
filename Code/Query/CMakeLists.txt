rdkit_headers(AndQuery.h
              EqualityQuery.h
              GreaterEqualQuery.h
              GreaterQuery.h
              LessEqualQuery.h
              LessQuery.h
              OrQuery.h
              Query.h
              QueryObjects.h
              RangeQuery.h
              SetQuery.h
              XOrQuery.h DEST Query)

rdkit_test(testQuery test.cpp LINK_LIBRARIES RDGeneral)

