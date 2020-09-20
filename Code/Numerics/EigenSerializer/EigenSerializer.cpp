// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
#include <Numerics/EigenSerializer/EigenSerializer.h>
#include <RDGeneral/Invariant.h>
namespace RDNumeric {
namespace EigenSerializer {
template <typename T>
bool serialize(const T& data, std::ofstream& ofs) {
  if (!ofs.is_open()) {
    return false;
  }
  { // create enclosing scope
    portable_binary_oarchive oa(ofs);
    oa << data;
  }
  ofs.close();
  return true;
}

template <typename T>
bool serialize(const T& data, const std::string& filename) {
  std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
  return serialize(data, ofs);
}

template <typename T>
bool deserialize(T& data, std::ifstream& ifs) {
  if (!ifs.is_open()) {
    return false;
  }
  { // create enclosing scope
    portable_binary_iarchive ia(ifs);
    ia >> data;
  }
  ifs.close();
  return true;
}

template <typename T>
bool deserialize(T& data, const std::string& filename) {
  std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
  return deserialize(data, ifs);
}

template <typename T>
bool deserializeAll(std::vector<T>& data, std::vector<std::string>& labels,
                    const std::string& filename) {
  std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
  return deserializeAll(data, labels, ifs);
}

template <typename T>
bool deserializeAll(std::vector<T>& data, std::vector<std::string>& labels,
                    std::ifstream& ifs) {
  if (!ifs.is_open()) {
    return false;
  }
  { // create enclosing scope
    portable_binary_iarchive ia(ifs);
    std::streampos archiveOffset = ifs.tellg();
    std::streampos streamEnd = ifs.seekg(0, std::ios_base::end).tellg();
    ifs.seekg(archiveOffset);
    while (ifs.tellg() < streamEnd) {
      std::string next_label;
      T next_data;
      ia >> next_label;
      ia >> next_data;
      labels.push_back(next_label);
      data.push_back(next_data);
    }
  }
  ifs.close();
  return true;
}

template <typename T>
bool serializeAll(const std::vector<T>& data,
                  const std::vector<std::string>& labels, std::ofstream& ofs) {
  if (!ofs.is_open()) {
    return false;
  }
  { // create enclosing scope
    portable_binary_oarchive oa(ofs);
    for (size_t i = 0; i < labels.size(); i++) {
      oa << labels[i];
      oa << data[i];
    }
  }
  return true;
}

template <typename T>
bool serializeAll(const std::vector<T>& data,
                  const std::vector<std::string>& labels,
                  const std::string& filename) {
  std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
  return serializeAll(data, labels, ofs);
}

// template specifications for Eigen types
template bool serialize<Eigen::ArrayXXd>(const Eigen::ArrayXXd&,
                                         const std::string&);
template bool serialize<Eigen::MatrixXd>(const Eigen::MatrixXd&,
                                         const std::string&);
template bool serialize<Eigen::ArrayXXf>(const Eigen::ArrayXXf&,
                                         const std::string&);
template bool serialize<Eigen::ArrayXd>(const Eigen::ArrayXd&,
                                        const std::string&);
template bool serialize<Eigen::MatrixXf>(const Eigen::MatrixXf&,
                                         const std::string&);

template bool deserialize<Eigen::ArrayXXd>(Eigen::ArrayXXd&,
                                           const std::string&);
template bool deserialize<Eigen::MatrixXd>(Eigen::MatrixXd&,
                                           const std::string&);
template bool deserialize<Eigen::ArrayXXf>(Eigen::ArrayXXf&,
                                           const std::string&);
template bool deserialize<Eigen::ArrayXd>(Eigen::ArrayXd&,
                                          const std::string&);
template bool deserialize<Eigen::MatrixXf>(Eigen::MatrixXf&,
                                           const std::string&);

template bool serialize<Eigen::ArrayXXd>(const Eigen::ArrayXXd&,
                                         std::ofstream& ofs);
template bool serialize<Eigen::MatrixXd>(const Eigen::MatrixXd&,
                                         std::ofstream& ofs);
template bool serialize<Eigen::ArrayXXf>(const Eigen::ArrayXXf&,
                                         std::ofstream& ofs);
template bool serialize<Eigen::ArrayXd>(const Eigen::ArrayXd&,
                                        std::ofstream& ofs);
template bool serialize<Eigen::MatrixXf>(const Eigen::MatrixXf&,
                                         std::ofstream& ofs);

template bool deserialize<Eigen::ArrayXXd>(Eigen::ArrayXXd&,
                                           std::ifstream& ifs);
template bool deserialize<Eigen::MatrixXd>(Eigen::MatrixXd&,
                                           std::ifstream& ifs);
template bool deserialize<Eigen::ArrayXXf>(Eigen::ArrayXXf&,
                                           std::ifstream& ifs);
template bool deserialize<Eigen::ArrayXd>(Eigen::ArrayXd&,
                                          std::ifstream& ifs);
template bool deserialize<Eigen::MatrixXf>(Eigen::MatrixXf&,
                                           std::ifstream& ifs);

template bool deserializeAll(std::vector<Eigen::ArrayXXd>& data,
                             std::vector<std::string>& labels,
                             const std::string& filename);
template bool deserializeAll(std::vector<Eigen::ArrayXXd>& data,
                             std::vector<std::string>& labels,
                             std::ifstream& ifs);
template bool serializeAll(const std::vector<Eigen::ArrayXXd>& data,
                           const std::vector<std::string>& labels,
                           const std::string& filename);
template bool serializeAll(const std::vector<Eigen::ArrayXXd>& data,
                           const std::vector<std::string>& labels,
                           std::ofstream& ofs);

template bool deserializeAll(std::vector<Eigen::ArrayXXf>& data,
                             std::vector<std::string>& labels,
                             const std::string& filename);
template bool deserializeAll(std::vector<Eigen::ArrayXXf>& data,
                             std::vector<std::string>& labels,
                             std::ifstream& ifs);
template bool serializeAll(const std::vector<Eigen::ArrayXXf>& data,
                           const std::vector<std::string>& labels,
                           const std::string& filename);
template bool serializeAll(const std::vector<Eigen::ArrayXXf>& data,
                           const std::vector<std::string>& labels,
                           std::ofstream& ofs);

template bool deserializeAll(std::vector<Eigen::MatrixXd>& data,
                             std::vector<std::string>& labels,
                             const std::string& filename);
template bool deserializeAll(std::vector<Eigen::MatrixXd>& data,
                             std::vector<std::string>& labels,
                             std::ifstream& ifs);
template bool serializeAll(const std::vector<Eigen::MatrixXd>& data,
                           const std::vector<std::string>& labels,
                           const std::string& filename);
template bool serializeAll(const std::vector<Eigen::MatrixXd>& data,
                           const std::vector<std::string>& labels,
                           std::ofstream& ofs);

template bool deserializeAll(std::vector<Eigen::MatrixXf>& data,
                             std::vector<std::string>& labels,
                             const std::string& filename);
template bool deserializeAll(std::vector<Eigen::MatrixXf>& data,
                             std::vector<std::string>& labels,
                             std::ifstream& ifs);
template bool serializeAll(const std::vector<Eigen::MatrixXf>& data,
                           const std::vector<std::string>& labels,
                           const std::string& filename);
template bool serializeAll(const std::vector<Eigen::MatrixXf>& data,
                           const std::vector<std::string>& labels,
                           std::ofstream& ofs);
}  // namespace EigenSerializer
}  // namespace RDNumeric