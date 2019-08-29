remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
add_definitions(-DRDKIT_NMMOLHASH_BUILD)

set(molhashlib_srcs
  hashfunctions.cpp
  normalize.cpp
  )

rdkit_headers(
   molhash.h mf.h
   DEST GraphMol/NMMolHash)



# Static library
rdkit_library(NMMolHashLib ${molhashlib_srcs}
              LINK_LIBRARIES SmilesParse GraphMol
  ${RDKit_THREAD_LIBS})

rdkit_catch_test(molHashCatchTest catch_tests.cpp LINK_LIBRARIES NMMolHashLib SmilesParse ${TOOLKIT_LIBS} )

