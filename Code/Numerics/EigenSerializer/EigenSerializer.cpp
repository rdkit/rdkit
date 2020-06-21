#include <Numerics/EigenSerializer/EigenSerializer.h>

namespace RDNumeric {
namespace EigenSerializer {
template <typename T>
bool serialize(const T& data, const std::string& filename) {
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open()) return false;
  {
    boost::archive::binary_oarchive oa(ofs);
    oa << data;
  }
  ofs.close();
  return true;
}

template <typename T>
bool deSerialize(T& data, const std::string& filename) {
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open()) return false;
  {
    boost::archive::binary_iarchive ia(ifs);
    ia >> data;
  }
  ifs.close();
  return true;
}

template bool serialize<Eigen::ArrayXXd>(const Eigen::ArrayXXd&,
                                         const std::string&);
template bool serialize<Eigen::MatrixXd>(const Eigen::MatrixXd&,
                                         const std::string&);

template bool deSerialize<Eigen::ArrayXXd>(Eigen::ArrayXXd&,
                                           const std::string&);
template bool deSerialize<Eigen::MatrixXd>(Eigen::MatrixXd&,
                                           const std::string&);
}  // namespace EigenSerializer
}  // namespace RDNumeric