CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.h.cmake
               ${CMAKE_CURRENT_SOURCE_DIR}/versions.h )

rdkit_library(RDGeneral
              Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp SHARED)

rdkit_headers(BadFileException.h
              Dict.h
              FileParseException.h
              Invariant.h
              RDLog.h
              StreamOps.h
              types.h
              utils.h
              versions.h
              LocaleSwitcher.h
              DEST RDGeneral)
if (NOT RDK_INSTALL_INTREE)
  install(DIRECTORY hash DESTINATION ${RDKit_HdrDir}/RDGeneral/hash
        PATTERN ".svn" EXCLUDE)
endif (NOT RDK_INSTALL_INTREE)
rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)

