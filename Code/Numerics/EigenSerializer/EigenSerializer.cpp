// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
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

template <typename T>
bool deSerializeAll(std::vector<T>* weights, std::vector<T>* biases,
                    std::string& filename, std::string atomType) {
  std::ifstream ifs(filename.c_str(), std::ios::in);
  if (!ifs.is_open()) return false;
  {
    boost::archive::binary_iarchive ia(ifs);

    std::streampos archiveOffset = ifs.tellg();
    std::streampos streamEnd = ifs.seekg(0, std::ios_base::end).tellg();
    ifs.seekg(archiveOffset);
    while (ifs.tellg() < streamEnd) {
      std::string weightType;
      T weight;
      ia >> weightType;
      ia >> weight;
      if (weightType.find(atomType) != std::string::npos) {
        if (weightType.find("bias") != std::string::npos) {
          biases->push_back(weight);
        }
        if (weightType.find("weight") != std::string::npos) {
          weights->push_back(weight);
        }
      }
    }
  }
  ifs.close();
  return true;
}

template bool serialize<Eigen::ArrayXXd>(const Eigen::ArrayXXd&,
                                         const std::string&);
template bool serialize<Eigen::MatrixXd>(const Eigen::MatrixXd&,
                                         const std::string&);
template bool serialize<Eigen::ArrayXXf>(const Eigen::ArrayXXf&,
                                         const std::string&);
template bool serialize<Eigen::ArrayXd>(const Eigen::ArrayXd&,
                                        const std::string&);

template bool deSerialize<Eigen::ArrayXXd>(Eigen::ArrayXXd&,
                                           const std::string&);
template bool deSerialize<Eigen::MatrixXd>(Eigen::MatrixXd&,
                                           const std::string&);
template bool deSerialize<Eigen::ArrayXXf>(Eigen::ArrayXXf&,
                                           const std::string&);
template bool deSerialize<Eigen::ArrayXd>(Eigen::ArrayXd&, const std::string&);

template bool deSerializeAll<Eigen::ArrayXXd>(std::vector<Eigen::ArrayXXd>*,
                                              std::vector<Eigen::ArrayXXd>*,
                                              std::string&, std::string);
template bool deSerializeAll<Eigen::ArrayXXf>(std::vector<Eigen::ArrayXXf>*,
                                              std::vector<Eigen::ArrayXXf>*,
                                              std::string&, std::string);

}  // namespace EigenSerializer
}  // namespace RDNumeric