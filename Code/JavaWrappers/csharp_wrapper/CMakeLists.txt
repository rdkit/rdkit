project (GraphMolCSharp)

include_directories( ${RDKit_ExternalDir} )

SET_SOURCE_FILES_PROPERTIES(GraphMolCSharp.i PROPERTIES CPLUSPLUS ON )

# Setup a few variables for environment-specific things
if(WIN32)
  ADD_DEFINITIONS("/W3 /wd4716 /bigobj")
  SET(PATH_SEP ";")
  SET(COPY_CMD xcopy ${COPY_SOURCE} ${COPY_DEST} /Y /I)
else()
  SET(PATH_SEP ":")
  SET(COPY_CMD cp -p ${COPY_SOURCE} ${COPY_DEST})
endif()

# Coax SWIG into playing nicely with Apple environments
if(APPLE)
  SET(CMAKE_SIZEOF_VOID_P 4)
endif(APPLE)

if(CMAKE_SIZEOF_VOID_P MATCHES 4)
  SET(CMAKE_SWIG_FLAGS -namespace "GraphMolWrap")
else()
  SET(CMAKE_SWIG_FLAGS -namespace "GraphMolWrap" "-DSWIGWORDSIZE64")
endif()
SET(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/swig_csharp )

if(RDK_BUILD_INCHI_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DBUILD_INCHI_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()
if(RDK_BUILD_AVALON_SUPPORT)
  SET(CMAKE_SWIG_FLAGS "-DBUILD_AVALON_SUPPORT" ${CMAKE_SWIG_FLAGS} )
endif()

FILE(GLOB SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../*.i")

# we added all source files, now remove the ones that we're not supporting in this build:
if(NOT RDK_BUILD_AVALON_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../AvalonLib.i")
endif()

if(NOT RDK_BUILD_INCHI_SUPPORT)
LIST(REMOVE_ITEM SWIG_SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../Inchi.i")
endif()

SET(SWIG_MODULE_GraphMolCSWrap_EXTRA_DEPS ${SWIG_SRC_FILES} )

SWIG_ADD_MODULE(GraphMolCSWrap "CSharp" GraphMolCSharp.i )

# it doesnt seem like the threading libs should need to be here, but
# as of Oct 2012 using boost 1.51 under at least ubuntu 12.04 we get a
# link error if they aren't there.
SWIG_LINK_LIBRARIES(GraphMolCSWrap ${RDKit_Wrapper_Libs}
      ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} )

INSTALL(TARGETS GraphMolCSWrap 
              DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )
