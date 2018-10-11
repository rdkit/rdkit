CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp.cmake
               ${CMAKE_CURRENT_SOURCE_DIR}/versions.cpp )
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/versions.h.cmake
               ${CMAKE_CURRENT_SOURCE_DIR}/versions.h )
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/RDConfig.h.cmake
              ${CMAKE_CURRENT_SOURCE_DIR}/RDConfig.h )


add_definitions(-DRDKIT_RDGENERAL_BUILD)
rdkit_library(RDGeneral
              Invariant.cpp types.cpp utils.cpp RDLog.cpp Dict.cpp
              LocaleSwitcher.cpp versions.cpp SHARED
              LINK_LIBRARIES ${RDKit_THREAD_LIBS})

if(MINGW)
  set(hasMSVCRuntime "FALSE")
  if(RDK_BUILD_THREADSAFE_SSS)
    set(needMSVCRuntime "TRUE")
  else()
    set(needMSVCRuntime "FALSE")
  endif()
  if(needMSVCRuntime AND MSVC_RUNTIME_DLL)
    if(EXISTS "${MSVC_RUNTIME_DLL}")
      set(hasMSVCRuntime "TRUE")
    endif()
  endif()
  if(needMSVCRuntime AND (NOT hasMSVCRuntime))
    set(systemRootPath "$ENV{SystemRoot}")
    string(REGEX REPLACE "\\\\" "/" systemRootPath "${systemRootPath}")
    message(FATAL_ERROR "Please set MSVC_RUNTIME_DLL to the full path to a msvcrXXX.dll "
      "Microsoft Visual C Runtime DLL appropriate for your system; the most likely path is "
      "${systemRootPath}/System32, or ${systemRootPath}/SysWOW64 in case you are building "
      "32-bit RDKit under a 64-bit Windows OS.")
  endif()
endif(MINGW)

target_link_libraries(RDGeneral PUBLIC ${RDKit_THREAD_LIBS} ${MSVC_RUNTIME_DLL})

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
              RDThreads.h
              StreamOps.h
              types.h
              utils.h
              versions.h
              RDConfig.h
              LocaleSwitcher.h
              Ranking.h
              hanoiSort.h
              export.h
              test.h
              DEST RDGeneral)
if (NOT RDK_INSTALL_INTREE)
  install(DIRECTORY hash DESTINATION ${RDKit_HdrDir}/RDGeneral/hash
        PATTERN ".svn" EXCLUDE)
endif (NOT RDK_INSTALL_INTREE)
rdkit_test(testDict testDict.cpp LINK_LIBRARIES RDGeneral)
rdkit_test(testRDValue testRDValue.cpp LINK_LIBRARIES RDGeneral)
