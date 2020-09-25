//
//  Copyright (C) 2020 Manan Goel, James Stevenson
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

#include <ForceField/Contrib.h>
#include <ForceField/ForceField.h>

#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>
using namespace Eigen;

namespace ForceFields {
namespace ANI {

//! Atomic contribution term from Atomic Environment Vectors using ANI models
class RDKIT_FORCEFIELD_EXPORT ANIAtomContrib : public ForceFieldContrib {
 public:
  ANIAtomContrib(){};
  //! Constructor
  /*!
    The contribution is found according to the AEV, atom type along with
    hyperparameters and weights, biases of the nueral network

    \param owner        pointer to the owning ForceField
    \param speciesVec   atomic number for each atom
    \param model        a parameter directory, such as one of those in
                        /Code/ForceField/ANI/Params/
  */
  ANIAtomContrib(ForceField *owner, VectorXi &speciesVec, std::string model);
  double getEnergy(double *pos) const;
  void getGrad(double *pos, double *grad) const;

  /*!
    Predict molecular energy by forward propagation through the neural network

    \param aev      Atomic Environment Vectors of the molecule

    \return         Molecule energy, not including atomic self-energy
  */
  double forwardProp(MatrixXd &aev) const;
  void forwardPropOneAtom(MatrixXd &layer_values, const int modelNo,
                            const int element_i) const;

  virtual ANIAtomContrib *copy() const { return new ANIAtomContrib(*this); };

 private:
  int d_aev_size; // AEV size, scales as n_elements^2
  // nested lists of parameters, in order: <ensemble<element<network>>>
  std::vector<std::vector<std::vector<MatrixXd>>> d_weights;
  std::vector<std::vector<std::vector<MatrixXd>>> d_biases;
  std::vector<double> d_selfEnergies;  // in element order
  VectorXi d_atomTypes;  // atom type for each atom (this, along with
                         // double *pos, represents the input molecule)
};

namespace Utils {}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
#endif
#endif