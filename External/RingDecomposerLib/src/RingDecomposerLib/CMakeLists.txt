include_directories(${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

rdkit_library(RingDecomposerLib
              RDLutility.c
              RDLhandler.c
              RDLrelation.c
              RDLapsp.c
              RingDecomposerLib.c
              RDLgraph.c
              RDLinfo.c
              RDLcycleFams.c
              RDLstack.c
              RDLtesting.c
              RDLtarjan.c
              RDLdimacs.c
              RDLbitset.c
              SHARED LINK_LIBRARIES RDGeneral )
rdkit_headers(RingDecomposerLib.h DEST GraphMol)
