

if(RDK_USE_BOOST_IOSTREAMS)
if(WIN32)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set (link_iostreams ${Boost_LIBRARIES})
  if (NOT Boost_USE_STATIC_LIBS)
     add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
  endif()
  find_package(ZLIB)
  include_directories(${ZLIB_INCLUDE_DIRS})
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS zlib)
  if(Boost_zlib_FOUND)
    set(zlib_lib Boost::zlib)
  else()
    set(zlib_lib ${ZLIB_LIBRARIES})
  endif()


else()
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS system iostreams REQUIRED)
  set (link_iostreams ${Boost_LIBRARIES})
  if (NOT Boost_USE_STATIC_LIBS)
     add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
  endif()
  find_package(ZLIB)
  include_directories(${ZLIB_INCLUDE_DIRS})
  set(zlib_lib ${ZLIB_LIBRARIES})
endif()
ADD_DEFINITIONS("-DRDK_USE_BOOST_IOSTREAMS")
endif(RDK_USE_BOOST_IOSTREAMS)

if(RDK_USE_BOOST_REGEX)
  add_definitions(-DRDKIT_USE_BOOST_REGEX)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS regex REQUIRED)
  set(regex_lib Boost::regex)
endif()

rdkit_library(MarvinParser
              MarvinParser.cpp
              MarvinWriter.cpp
              MarvinDefs.cpp
              ${maesupplier}
              LINK_LIBRARIES GenericGroups Depictor GraphMol ChemReactions ${RDK_MAEPARSER_LIBS} ${regex_lib} ${link_iostreams} ${zlib_lib})
target_compile_definitions(MarvinParser PRIVATE RDKIT_MARVINPARSER_BUILD)

rdkit_headers(MarvinParser.h
              MarvinDefs.h
              DEST GraphMol/MarvinParser)

rdkit_test(testMrvToMol testMrvToMol.cpp
           LINK_LIBRARIES MarvinParser FileParsers SubstructMatch )


