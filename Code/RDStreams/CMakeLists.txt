
# gzip is backed by zlib (required); bzip2 by libbz2 when RDK_USE_BZIP2 is set
# (libbz2 is located at the top level). Linking these directly means RDStreams
# needs no compiled compression-stream library.
find_package(ZLIB REQUIRED)
set(RDStreams_dep_libs ZLIB::ZLIB)
if(RDK_USE_BZIP2)
  list(APPEND RDStreams_dep_libs BZip2::BZip2)
endif()

rdkit_library(RDStreams streams.cpp LINK_LIBRARIES ${RDStreams_dep_libs})
target_compile_definitions(RDStreams PRIVATE RDKIT_RDSTREAMS_BUILD)

rdkit_headers(streams.h DEST RDStreams)
