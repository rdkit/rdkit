#include <RDGeneral/export.h>
// #ifndef __RD_ANI_H__
// #define __RD_ANI_H__

#include <ForceField/ForceField.h>
#include <ForceField/Contrib.h>

// #ifdef RDK_HAS_EIGEN3
#include <eigen3/Eigen/Dense>
#include <boost/tokenizer.hpp>
using namespace Eigen;

namespace ForceFields {
namespace ANI {

class RDKIT_FORCEFIELD_EXPORT ANIAtomContrib : public ForceFieldContrib {
 public:
  ANIAtomContrib(){};

  ANIAtomContrib(ForceField *owner, int atomType, unsigned int atomIdx,
                 VectorXi speciesVec, int numAtoms);
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
  std::vector<ArrayXXd> d_weights;
  std::vector<ArrayXXd> d_biases;
};

namespace Utils {

void CELU(ArrayXXd &input, double alpha);

void loadFromCSV(ArrayXXd *param, unsigned int model, std::string type, unsigned int layer, char atomType);

std::vector<std::string> tokenize(const std::string &s);

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
   // #endif
   // #endif