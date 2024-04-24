project (GraphMolCSharp)

include_directories( ${RDKit_ExternalDir} )
SET_SOURCE_FILES_PROPERTIES(GraphMolCSharp.i PROPERTIES CPLUSPLUS ON )

# Setup a few variables for environment-specific things
if(MSVC)
  ADD_COMPILE_OPTIONS(/W3 /wd4716 /bigobj)
  SET(PATH_SEP ";")
  SET(COPY_CMD xcopy ${COPY_SOURCE} ${COPY_DEST} /Y /I)
else()
  SET(PATH_SEP ":")
  SET(COPY_CMD cp -p ${COPY_SOURCE} ${COPY_DEST})
endif()

# Coax SWIG into playing nicely with Apple/ARM environments
set(HAS_ARM_PROC FALSE)
string(TOLOWER ${CMAKE_SYSTEM_PROCESSOR} LOWERCASE_SYS_PROC)
if(LOWERCASE_SYS_PROC MATCHES "arm")
  set(HAS_ARM_PROC TRUE)
endif()
# Based on https://github.com/mono/mono/blob/5d2e3bc3b3c8184d35b2f7801e88d96470d367c4/mcs/mcs/settings.cs#L54
# Mono seems to only support 32-bit ARM
if(APPLE OR HAS_ARM_PROC)
  SET(CMAKE_SIZEOF_VOID_P 4)
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  SET(CMAKE_SWIG_FLAGS -namespace "GraphMolWrap")
  if(HAS_ARM_PROC)
    SET(PLATFORM "arm")
  else()
    SET(PLATFORM "x86")
  endif()
else()
  SET(PLATFORM "x64")
  if (WIN32)
    SET(CMAKE_SWIG_FLAGS -namespace "GraphMolWrap")
  else()
    SET(CMAKE_SWIG_FLAGS -namespace "GraphMolWrap" "-DSWIGWORDSIZE64")
  endif()
endif()
SET(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/swig_csharp )

if(RDK_BUILD_INCHI_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DRDK_BUILD_INCHI_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()
if(RDK_BUILD_AVALON_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DRDK_BUILD_AVALON_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()
if(RDK_USE_BOOST_IOSTREAMS)
  SET(CMAKE_SWIG_FLAGS "-DRDK_USE_BOOST_IOSTREAMS" ${CMAKE_SWIG_FLAGS} )
endif()
if (RDK_BUILD_CAIRO_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DRDK_BUILD_CAIRO_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()
  

FILE(GLOB SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../*.i")

# we added all source files, now remove the ones that we're not supporting in this build:
if(NOT RDK_BUILD_AVALON_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../AvalonLib.i")
endif()

if(NOT RDK_BUILD_INCHI_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../Inchi.i")
endif()

SET(SWIG_MODULE_RDKFuncs_EXTRA_DEPS ${SWIG_SRC_FILES} )

SWIG_ADD_LIBRARY(RDKFuncs TYPE MODULE LANGUAGE CSharp SOURCES GraphMolCSharp.i )


# it doesnt seem like the threading libs should need to be here, but
# as of Oct 2012 using boost 1.51 under at least ubuntu 12.04 we get a
# link error if they aren't there.
SWIG_LINK_LIBRARIES(RDKFuncs ${RDKit_Wrapper_Libs}
      rdkit_base ${RDKit_THREAD_LIBS} )

INSTALL(TARGETS RDKFuncs
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

if (RDK_BUILD_DOTNET_CSHARP_TESTS OR MSVC)

  find_program(DOTNET_EXE dotnet)
  if (NOT DOTNET_EXE)
    MESSAGE("dotnet (executable) is not found.  Pleae add it to PATH and rerun cmake.")
    MESSAGE(FATAL_ERROR "Cannot find required executable dotnet")
  endif (NOT DOTNET_EXE)

  ADD_CUSTOM_COMMAND(
    OUTPUT RDKitDotNetLib SYMBOLIC
  
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:RDKFuncs> ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${DOTNET_EXE} build RDKitDotNetCore.sln
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS RDKFuncs
  )
  
  ADD_CUSTOM_TARGET(RDKitDotNet ALL
     DEPENDS RDKFuncs RDKitDotNetLib
     COMMENT "Building RDKit DotNet project"
  )
  
  ADD_TEST(
    NAME CSharpTests
    COMMAND ${DOTNET_EXE} test RDKitDotNetCore.sln "-l:console$<SEMICOLON>verbosity=normal" 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
else(RDK_BUILD_DOTNET_CSHARP_TESTS OR MSVC)
  # find the mcs executable on non-windows systems, if dotnet is not used:
  find_program(GMCS_EXE mcs)
  if (NOT GMCS_EXE)
	MESSAGE ("mcs (executable) is not found. Please add it to PATH and rerun cmake.")
	MESSAGE(FATAL_ERROR "Cannot find required executable mcs")
  endif (NOT GMCS_EXE)

  # code adapted from the wrapper code for
  # GDCM: http://gdcm.svn.sf.net/viewvc/gdcm/trunk/Wrapping/Java/CMakeLists.txt?view=markup
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/RDKit2DotNet.dll
    ## 1. run this custom command only after swig has been run.
    COMMAND ${GMCS_EXE} -platform:${PLATFORM} -out:RDKit2DotNet.dll -t:library "swig_csharp/*.cs"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS "${swig_generated_file_fullname}"
  )
  ADD_CUSTOM_TARGET(RDKFuncsDLL ALL
    DEPENDS RDKFuncs ${CMAKE_CURRENT_SOURCE_DIR}/RDKit2DotNet.dll
    COMMENT "building mono dll"
  )
endif(RDK_BUILD_DOTNET_CSHARP_TESTS OR MSVC)

