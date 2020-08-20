// https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
#include <Numerics/EigenSerializer/EigenSerializerNonPortable.h>
#include <RDGeneral/Invariant.h>
namespace RDNumeric {
namespace EigenSerializerNonPortable {
template <typename T>
bool serialize(const T& data, std::ofstream& ofs) {
  if (!ofs.is_open()) {
    return false;
  }
  {
    boost::archive::binary_oarchive oa(ofs);
    oa << data;
  }
  ofs.close();
  return true;
}

template <typename T>
bool serialize(const T& data, const std::string& filename) {
  std::ofstream ofs(filename.c_str(), std::ios::out);
  return serialize(data, ofs);
}

template <typename T>
bool deserialize(T& data, std::ifstream& ifs) {
  if (!ifs.is_open()) {
    return false;
  }
  {
    boost::archive::binary_iarchive ia(ifs);
    ia >> data;
  }
  ifs.close();
  return true;
}

template <typename T>
bool deserialize(T& data, const std::string& filename) {
  std::ifstream ifs(filename.c_str(), std::ios::in);
  return deserialize(data, ifs);
}

template <typename T>
bool deserializeAll(std::vector<T>* weights, std::vector<T>* biases,
                    std::string& filename, std::string atomType) {
  PRECONDITION(weights != nullptr, "Weights Array is NULL");
  PRECONDITION(biases != nullptr, "Biases Array is NULL");
  std::ifstream ifs(filename.c_str(), std::ios::in);
  return deserializeAll(weights, biases, ifs, atomType);
}

template <typename T>
bool deserializeAll(std::vector<T>* weights, std::vector<T>* biases,
                    std::ifstream& ifs, std::string atomType) {
  PRECONDITION(weights != nullptr, "Weights Array is NULL");
  PRECONDITION(biases != nullptr, "Biases Array is NULL");
  if (!ifs.is_open()) {
    return false;
  }
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

template <typename T>
bool serializeAll(
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, T>>>>*
        weightsAndBiasesForEachAtomType,
    std::ofstream& ofs) {
  PRECONDITION(weightsAndBiasesForEachAtomType != nullptr,
               "Array of Weights and Biases is NULL");
  if (!ofs.is_open()) {
    return false;
  }
  {
    boost::archive::binary_oarchive oa(ofs);
    for (unsigned int i = 0; i < weightsAndBiasesForEachAtomType->size(); i++) {
      auto atomType = (*weightsAndBiasesForEachAtomType)[i].first;
      auto weights = (*weightsAndBiasesForEachAtomType)[i].second;
      for (unsigned int j = 0; j < weights.size(); j++) {
        std::string identifier =
            atomType + "_" + std::to_string(j / 2) + "_" + weights[j].first;
        oa << identifier;
        oa << weights[j].second;
      }
    }
  }
  return true;
}

template <typename T>
bool serializeAll(
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, T>>>>*
        weightsAndBiasesForEachAtomType,
    const std::string& fileName) {
  std::ofstream ofs(fileName.c_str(), std::ios::out);
  return serializeAll(weightsAndBiasesForEachAtomType, ofs);
}

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
template bool deserialize<Eigen::ArrayXd>(Eigen::ArrayXd&, const std::string&);
template bool deserialize<Eigen::MatrixXf>(Eigen::MatrixXf&,
                                           const std::string&);

template bool deserializeAll<Eigen::ArrayXXd>(std::vector<Eigen::ArrayXXd>*,
                                              std::vector<Eigen::ArrayXXd>*,
                                              std::string&, std::string);
template bool deserializeAll<Eigen::ArrayXXf>(std::vector<Eigen::ArrayXXf>*,
                                              std::vector<Eigen::ArrayXXf>*,
                                              std::string&, std::string);
template bool deserializeAll<Eigen::MatrixXf>(std::vector<Eigen::MatrixXf>*,
                                              std::vector<Eigen::MatrixXf>*,
                                              std::string&, std::string);
template bool deserializeAll<Eigen::MatrixXd>(std::vector<Eigen::MatrixXd>*,
                                              std::vector<Eigen::MatrixXd>*,
                                              std::string&, std::string);

template bool deserializeAll<Eigen::ArrayXXd>(std::vector<Eigen::ArrayXXd>*,
                                              std::vector<Eigen::ArrayXXd>*,
                                              std::ifstream&, std::string);
template bool deserializeAll<Eigen::ArrayXXf>(std::vector<Eigen::ArrayXXf>*,
                                              std::vector<Eigen::ArrayXXf>*,
                                              std::ifstream&, std::string);
template bool deserializeAll<Eigen::MatrixXd>(std::vector<Eigen::MatrixXd>*,
                                              std::vector<Eigen::MatrixXd>*,
                                              std::ifstream&, std::string);
template bool deserializeAll<Eigen::MatrixXf>(std::vector<Eigen::MatrixXf>*,
                                              std::vector<Eigen::MatrixXf>*,
                                              std::ifstream&, std::string);

template bool serializeAll<Eigen::ArrayXXf>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::ArrayXXf>>>>*,
    const std::string&);
template bool serializeAll<Eigen::ArrayXXd>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::ArrayXXd>>>>*,
    const std::string&);
template bool serializeAll<Eigen::MatrixXf>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::MatrixXf>>>>*,
    const std::string&);
template bool serializeAll<Eigen::MatrixXd>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::MatrixXd>>>>*,
    const std::string&);

template bool serializeAll<Eigen::ArrayXXf>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::ArrayXXf>>>>*,
    std::ofstream&);
template bool serializeAll<Eigen::ArrayXXd>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::ArrayXXd>>>>*,
    std::ofstream&);
template bool serializeAll<Eigen::MatrixXf>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::MatrixXf>>>>*,
    std::ofstream&);
template bool serializeAll<Eigen::MatrixXd>(
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, Eigen::MatrixXd>>>>*,
    std::ofstream&);

}  // namespace EigenSerializer
}  // namespace RDNumeric