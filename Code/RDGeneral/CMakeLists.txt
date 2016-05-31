CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp.cmake
               ${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp )

rdkit_library(RDGeneral
              Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp
              LocaleSwitcher.cpp versions.cpp SHARED
              LINK_LIBRARIES ${RDKit_THREAD_LIBS})

target_link_libraries(RDGeneral ${RDKit_THREAD_LIBS})

rdkit_headers(Exceptions.h
              BadFileException.h
              BoostStartInclude.h
              BoostEndInclude.h
              Dict.h
              FileParseException.h
              Invariant.h
              RDAny.h
              RDValue.h
              RDValue-doublemagic.h
              RDValue-taggedunion.h
              RDLog.h
              RDProps.h
              StreamOps.h
              tags.h
              types.h
              utils.h
              versions.h
              LocaleSwitcher.h
              Ranking.h
              hanoiSort.h
              DEST RDGeneral)
if (NOT RDK_INSTALL_INTREE)
  install(DIRECTORY hash DESTINATION ${RDKit_HdrDir}/RDGeneral/hash
        PATTERN ".svn" EXCLUDE)
endif (NOT RDK_INSTALL_INTREE)
rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)
rdkit_test(testRDValue testRDValue.cpp LINK_LIBRARIES RDGeneral)
rdkit_test(testRDTags testRDTags.cpp LINK_LIBRARIES RDGeneral)

