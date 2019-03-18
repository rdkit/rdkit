find_package(Boost 1.56.0 COMPONENTS system iostreams REQUIRED)
set (link_iostreams ${Boost_LIBRARIES})
if (NOT Boost_USE_STATIC_LIBS)
   add_definitions("-DBOOST_IOSTREAMS_DYN_LINK")
endif()

if (WIN32)
  find_package(Boost 1.58.0 COMPONENTS zlib REQUIRED)
  set(zlib_lib Boost::zlib)
endif()

add_definitions(-DRDKIT_RDSTREAMS_BUILD)
rdkit_library(RDStreams streams.cpp
              LINK_LIBRARIES ${Boost_LIBRARIES} ${link_iostreams} ${zlib_lib}
              DEST RDStreams)


