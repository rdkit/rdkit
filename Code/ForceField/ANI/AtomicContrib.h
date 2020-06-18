#include <RDGeneral/export.h>
#ifndef __RD_ANI_H__
#define __RD_ANI_H__

#include <ForceField/ForceField.h>
#include <ForceField/Contrib.h>

#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>
using namespace Eigen;

namespace ForceFields {
namespace ANI {

class RDKIT_FORCEFIELD_EXPORT ANIAtomContrib : public ForceFieldContrib {
 public:
  ANIAtomContrib(){};

  ANIAtomContrib(ForceField *owner, int atomType, size_t atomIdx,
                 VectorXi speciesVec, int numAtoms, size_t numLayers,
                 size_t ensembleSize, std::string modelType);
  double getEnergy(double *pos) const;
  double getEnergy(Eigen::ArrayXXd aev) const;
  void getGrad(double *pos, double *grad) const;

  double forwardProp(ArrayXXd aev) const;

  virtual ANIAtomContrib *copy() const { return new ANIAtomContrib(*this); };

 private:
  int d_atomType;
  int d_atomIdx;
  int d_numAtoms;
  VectorXi d_speciesVec;
  std::vector<std::vector<ArrayXXd>> d_weights;
  std::vector<std::vector<ArrayXXd>> d_biases;
  double d_selfEnergy;
  size_t d_ensembleSize;
  std::string d_modelType;
  std::map<int, std::string> d_atomEncoding = {
      {0, "H"}, {1, "C"}, {2, "N"}, {3, "O"}};
};

namespace Utils {

void CELU(ArrayXXd &input, double alpha);

void loadFromCSV(std::vector<ArrayXXd> *weights, size_t model,
                 std::string weightType, size_t layer, std::string atomType,
                 std::string modelType);

void loadSelfEnergy(double *energy, std::string atomType,
                    std::string modelType);

std::vector<std::string> tokenize(const std::string &s);

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
#endif
#endif