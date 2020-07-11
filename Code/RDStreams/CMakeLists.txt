if(RDK_USE_BOOST_IOSTREAMS)
ADD_DEFINITIONS("-DRDK_USE_BOOST_IOSTREAMS")

find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS iostreams REQUIRED)
if (NOT Boost_USE_STATIC_LIBS)
   add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
endif()

if (WIN32)
  find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS zlib)
  if(Boost_zlib_FOUND)
    set(zlib_lib Boost::zlib)
  endif()
else()
  if(Boost_USE_STATIC_LIBS)
    # when we're doing static linkage of boost
    # to also link against zlib (due to iostreams)
    find_package(ZLIB)
    set(zlib_lib ${ZLIB_LIBRARIES})
  endif()
endif()

add_definitions(-DRDKIT_RDSTREAMS_BUILD)
set(boost_iostreams Boost::iostreams)
endif(RDK_USE_BOOST_IOSTREAMS)
rdkit_library(RDStreams streams.cpp
              LINK_LIBRARIES ${boost_iostrems} ${zlib_lib})
rdkit_headers(streams.h DEST RDStreams)
