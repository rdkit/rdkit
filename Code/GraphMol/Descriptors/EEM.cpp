//
//  Copyright (c) 2016, Guillaume GODIN
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Institue of Cancer Research.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// for build & set RDBASE! => export RDBASE=/Users/GVALMTGG/Github/rdkit_mine/

#include <GraphMol/RDKitBase.h>

#include "WHIM.h"
#include "MolData3Ddescriptors.h"

#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace Eigen;

namespace RDKit {
namespace Descriptors {
namespace {

const float kappa = 0.1960;

const float A1[] = {0.0,2.3594,0.0,0.0,0.0,0.0,2.4541,2.5908,2.7130,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.3833};
const float A2[] = {0.0,0.0,0.0,0.0,0.0,0.0,2.4726,2.5409,2.6766,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.4956};
const float A3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

const float B1[] = {0.0,0.5962,0.0,0.0,0.0,0.0,0.2591,0.3316,0.5028,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.4564};
const float B2[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.2268,0.2319,0.4992,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1493};
const float B3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

// function to retreive the atomtype value based on the highest bond type of an atom
// in the publication they don't have access to "Aromatic type"
double getAtomtype(const RDKit::ROMol mol, const RDKit::Atom *atom) {
      double t=1.0;
      double a;
      RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
      while (nbrIdx != endNbrs) {
        const RDKit::Bond *bond = mol.getBondBetweenAtoms(*nbrIdx, atom->getIdx());
        a = bond->getBondTypeAsDouble();
        if (a == 1.5) {
             t = 2.0;
            }
        if (a == 2.0) {
             t = 2.0;
            }
        if (a == 3.0) {
             t = 3.0;
            }
        ++nbrIdx;
  }
  return t;
}


double* getEEMMatrix(RDKit::ROMol mol, double * dist3D, int n) {
  int sizeArray = (n+1) * (n+1);
  double* EEM = new double[sizeArray];
  #define IDX(x, y) (x * (n+1) + y)
    /* Fill the full n * n block */
    for(long int i = 0; i < n; i++) {
        int idx=mol.getAtomWithIdx(i)->getAtomicNum();
        double t = getAtomtype(mol,mol.getAtomWithIdx(i));
        double v;
        if (t== 1.0) {
                v=  B1[idx];
        }
        if (t== 2.0) {
                v = B2[idx];
        }
        if (t== 3.0) {
                v = B3[idx];
        }
        //std::cout << "v: " << v <<  "\n";

        EEM[IDX(i, i)] = v;
        for(long int j = i + 1; j < n; j++) {
                EEM[IDX(i, j)] = kappa / dist3D[i*n+j];
                EEM[IDX(j, i)] = EEM[IDX(i, j)];
            }
        }
    /* Fill last column & row */
    for(long int i = 0; i < n; i++) {
        EEM[IDX(n, i)] = 1.0; // column
        EEM[IDX(i, n)] = -1.0; // row

    }
    /* Set the bottom right element to zero */
    EEM[(n+1)*(n+1)+n] = 0.0;
    #undef IDX
    return EEM;
}


double* getBVector(RDKit::ROMol mol, int n) {
       /* Fill vector b i.e. -A */
        double * b = new double[n+1];
        for(int j = 0; j < n; j++) {
            double t = getAtomtype(mol,mol.getAtomWithIdx(j));
            int idx=mol.getAtomWithIdx(j)->getAtomicNum();
            if (t== 1.0) {
                b[j] = -A1[idx];
            }
            if (t== 2.0) {
                b[j] = -A2[idx];
            }
            if (t== 3.0) {
                b[j] = -A3[idx];
            }
        }

        b[n] = RDKit::MolOps::getFormalCharge(mol); // sum of charges

    return b;
}






/* Calculate charges for a particular kappa_data structure */
std::vector<double>  calculate_charges(RDKit::ROMol mol, double * dist3D, int numAtoms) {
        double * A = getEEMMatrix(mol,dist3D,numAtoms);
        double * b = getBVector(mol,  numAtoms);

        Map<MatrixXd> AM(A, numAtoms+1, numAtoms+1);
        Map<VectorXd> bv(b, numAtoms+1);
        VectorXd Res(numAtoms+1);
        Res =  AM.jacobiSvd(ComputeThinU | ComputeThinV).solve(bv); // SDV linear equation solver

        std::vector<double> charges;
        charges.resize(Res.size());
        VectorXd::Map(&charges[0], Res.size()) = Res;

        for(int j = 0; j < numAtoms; j++) {
            std::cout << charges[j] <<  ", ";
        }
        std::cout << "\n";

   free(A);
   free(b);


    return charges;

}



void getEEMs(const ROMol &mol, std::vector<double> &result, int numAtoms) {
  std::vector<double> wu(numAtoms);

  int confId=-1;
  // 3D distance matrix
  double * dist3D = RDKit::MolOps::get3DDistanceMat(mol, confId, false, true);

  wu = calculate_charges(mol, dist3D, numAtoms);

  result.clear();
  result.resize(numAtoms);

  for (int i = 0; i < numAtoms; i++) {
    result[i] = wu[i];
  }
  wu.clear();
}


}  // end of anonymous namespace

void EEM(const ROMol &mol, std::vector<double> &res) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(numAtoms);
  getEEMs(mol, res, numAtoms);//, confId);
}
}  // end of Descriptors namespace
}  // end of RDKit namespace