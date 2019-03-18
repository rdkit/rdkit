#include "streams.h"

namespace RDKit
{
gzstream::gzstream(const std::string &fname) :
    boost::iostreams::filtering_istream(),
    is(fname.c_str(), std::ios_base::binary) {
  push(boost::iostreams::gzip_decompressor());
  push(is);
}
}
