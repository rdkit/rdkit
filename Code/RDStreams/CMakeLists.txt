if(RDK_USE_BOOST_IOSTREAMS)
ADD_DEFINITIONS("-DRDK_USE_BOOST_IOSTREAMS")

find_package(Boost 1.56.0 COMPONENTS system iostreams REQUIRED)
set (link_iostreams ${Boost_LIBRARIES})
if (NOT Boost_USE_STATIC_LIBS)
   add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
endif()

if (WIN32)
  find_package(Boost 1.58.0 COMPONENTS zlib)
  if(Boost_zlib_FOUND)
    set(zlib_lib Boost::zlib)
  endif()
else()
  if(Boost_USE_STATIC_LIBS)
    # when we're doing static linkage of boost, at least on windows, we need
    # to also link against zlib (due to iostreams)
    find_package(ZLIB)
    set(zlib_lib ${ZLIB_LIBRARIES})
  endif()
endif()

add_definitions(-DRDKIT_RDSTREAMS_BUILD)
endif(RDK_USE_BOOST_IOSTREAMS)
rdkit_library(RDStreams streams.cpp
              LINK_LIBRARIES ${Boost_LIBRARIES} ${link_iostreams} ${zlib_lib}
              DEST RDStreams)
