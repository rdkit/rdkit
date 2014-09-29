/*
 * PMIDescriptors.cpp
 *
 *  Created on: Jun 25, 2014
 *      Author: hahnda6
 */

#include "PMIDescriptors.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/MolAlign/AlignMolecules.h>

namespace RDKit {
  void Descriptors::calcPMIDescriptors(const ROMol &mol, Descriptors::Moments &MoI,
                                                 int confId, RDNumeric::DoubleVector *weights,
                                                 unsigned int maxIterations){
    std::vector<double> eigenVals(3, 0.0);
    std::vector< std::vector<double> >eigenVecs(3, std::vector<double>(3, 0.0));

    MolAlign::getMomentsOfInertia(mol, eigenVals, eigenVecs, confId, weights, maxIterations);
    MoI.PMI1 = eigenVals[0];
    MoI.PMI2 = eigenVals[1];
    MoI.PMI3 = eigenVals[2];
    MoI.NPR1 = MoI.PMI1/MoI.PMI3;
    MoI.NPR2 = MoI.PMI2/MoI.PMI3;
  }
}
