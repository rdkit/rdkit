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
    \param atomType     ANI currently supports H.C, N and O
    \param atomIdx      index of atom in ForceField's positions
    \param speciesVec   vector with atom wise encoding
    \param numAtoms     number of atoms in the molecule
    \param numLayers    number of layers in the neural network
    \param ensembleSize number of models in the ensemble
    \param modelType    model types like ANI-1x and ANI-1ccx
  */
  ANIAtomContrib(ForceField *owner, int atomType, unsigned int atomIdx,
                 VectorXi &speciesVec, unsigned int numAtoms,
                 unsigned int numLayers, unsigned int ensembleSize,
                 std::string modelType);
  double getEnergy(double *pos) const;
  double getEnergy(Eigen::ArrayXXd &aev) const;
  void getGrad(double *pos, double *grad) const;

  /*!
    Find atomic contribution according to atom's interactions with other atoms
    by forward prpogation through the neural network

    \param aev      Atomic Environment Vector of the atom

    \return         Contribtution to the total energy of molecule according to
                    interactions
  */
  double forwardProp(ArrayXXd &aev) const;

  virtual ANIAtomContrib *copy() const { return new ANIAtomContrib(*this); };

 private:
  int d_atomType;
  int d_atomIdx;
  int d_numAtoms;
  VectorXi d_speciesVec;
  std::vector<std::vector<MatrixXd>> d_weights;
  std::vector<std::vector<MatrixXd>> d_biases;
  double d_selfEnergy;
  unsigned int d_ensembleSize;
  std::string d_modelType;
  std::map<int, std::string> d_atomEncoding = {
      {0, "H"}, {1, "C"}, {2, "N"}, {3, "O"}};
  std::map<std::string, ArrayXXd> d_aevParams;
};

namespace Utils {

//! CELU activation function
/*!
  Continuously Differentiable Exponential Linear Unit Activation function
  \param input      Vector from hidden layer
  \param alpha      hyperparameter for CELU
*/
void CELU(MatrixXd &input, double alpha);

/*!
  Load model weights from CSV file
  \param weights    Pointer to array of weights of neural network
  \param model      Index of model in the ensemble
  \param weightType Type of weights "weight" or "bias"
  \param layer      Index of layer of NN
  \param atomType   Atomic Symbol of atom
  \param modelType  Architecture being used
*/
void loadFromCSV(std::vector<MatrixXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType);

/*!
  Load model weights boost serialized files
  \param weights    Pointer to array of weights of neural network
  \param model      Index of model in the ensemble
  \param weightType Type of weights "weight" or "bias"
  \param layer      Index of layer of NN
  \param atomType   Atomic Symbol of atom
  \param modelType  Architecture being used
*/
void loadFromBin(std::vector<MatrixXd> *weights, unsigned int model,
                 std::string weightType, unsigned int layer,
                 std::string atomType, std::string modelType);

/*!
  Load all model weights from a single boost serialized file
  \param weights    Pointer to array of weights of neural network
  \param biases     Pointer to array of biases of neural network
  \param model      Index of model in the ensemble
  \param atomType   Atomic Symbol of atom
  \param modelType  Architecture being used
*/
void loadFromBin(std::vector<MatrixXd> *weights, std::vector<MatrixXd> *biases,
                 unsigned int model, std::string atomType,
                 std::string modelType);

//! Load self energy of atom from modelType/selfEnergies file
void loadSelfEnergy(double *energy, std::string atomType,
                    std::string modelType);

std::vector<std::string> tokenize(const std::string &s);

}  // namespace Utils
}  // namespace ANI
}  // namespace ForceFields
#endif
#endif