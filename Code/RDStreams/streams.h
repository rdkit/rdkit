//
#include <RDGeneral/export.h>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace RDKit {
// gzstream from a file
class RDKIT_RDSTREAMS_EXPORT gzstream
    : public boost::iostreams::filtering_istream {
  std::ifstream is;

 public:
  gzstream(const std::string &fname);
};
}  // namespace RDKit
