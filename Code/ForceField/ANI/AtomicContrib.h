//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef __RD_ANI_H__
#define __RD_ANI_H__

#include <ForceField/ForceField.h>
#include <ForceField/Contrib.h>

// #ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>
using namespace Eigen;

namespace ForceFields {
namespace ANI {

class RDKIT_FORCEFIELD_EXPORT ANIAtomContrib : public ForceFieldContrib {
 public:
  ANIAtomContrib(){};

  ANIAtomContrib(ForceField *owner, int atomType, unsigned int atomIdx,
                 VectorXi speciesVec, unsigned int numAtoms,
                 unsigned int numLayers, unsigned int ensembleSize,
                 std::string modelType);
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
  unsigned int d_ensembleSize;
  std::string d_modelType;
  std::map<int, std::string> d_atomEncoding = {
      {0, "H"}, {1, "C"}, {2, "N"}, {3, "O"}};
};

namespace Utils {

void CELU(ArrayXXd &input, double alpha);

void loadFromCSV(std::vector<ArrayXXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType);


void loadFromBin(std::vector<ArrayXXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType);

void loadSelfEnergy(double *energy, std::string atomType,
                    std::string modelType);

std::vector<std::string> tokenize(const std::string &s);

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
#endif
// #endif