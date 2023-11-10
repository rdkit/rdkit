
rdkit_library(RDStreams streams.cpp)
target_compile_definitions(RDStreams PRIVATE RDKIT_RDSTREAMS_BUILD)

rdkit_headers(streams.h DEST RDStreams)
