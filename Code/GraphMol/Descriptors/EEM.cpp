//
//  Copyright (c) 2017, Guillaume GODIN
//  "Copyright 2013-2016 Tomas Racek (tom@krab1k.net)"
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

#include "EEM.h"
#include "MolData3Ddescriptors.h"
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>


using namespace Eigen;

namespace RDKit {
namespace Descriptors {
namespace {


struct EEM_arrays {
  int* Atomindex;
  double* EEMatomtype;
};


// Those Parameters change to adapted to the molecule dataset using the optimization method:
// The best parameters published in the NEEMP article (https://jcheminf.springeropen.com/articles/10.1186/s13321-016-0171-1),
// are CCD_gen_DE_RMSD_B3LYP_6311G_NPA.par from Additional file 6 was built using the larger dataset with B3LYP_6311G ab initio => 17769 molecules
// H 1, C 6, N 7, O 8, F 9, P 15, S 16, Cl 17, Br 35
// Caution: Iodine is not in this training set of molecules!!!
const float kappa = 0.5125;
                   //0   1      2   3   4   5   6       7      8      9     10  11  12  13  14  15      16    17    18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
const float A1[] = {0.0,2.5473,0.0,0.0,0.0,0.0,2.7221,2.9750,3.1503,2.9976,0.0,0.0,0.0,0.0,0.0,0.0,   2.6511,2.7026,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.6263};
const float A2[] = {0.0,0.0,   0.0,0.0,0.0,0.0,2.7667,2.8895,3.0486,0.0,   0.0,0.0,0.0,0.0,0.0,2.2933,2.6471,   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
const float A3[] = {0.0,0.0,   0.0,0.0,0.0,0.0,2.6944,3.0240,0.0,   0.0,   0.0,0.0,0.0,0.0,0.0,0.0,   0.0,      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                   //0   1      2   3   4   5   6       7      8      9     10  11  12  13  14  15      16    17     18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35
const float B1[] = {0.0,1.1641,0.0,0.0,0.0,0.0,0.6403,0.9083,1.0577,0.9983,0.0,0.0,0.0,0.0,0.0,0.0,   0.4897,1.1537,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.1105};
const float B2[] = {0.0,0.0,   0.0,0.0,0.0,0.0,0.6513,0.6647,0.8410,0.0,   0.0,0.0,0.0,0.0,0.0,0.5759,0.4512,   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
const float B3[] = {0.0,0.0,   0.0,0.0,0.0,0.0,0.6776,1.4240,0.0,   0.0,   0.0,0.0,0.0,0.0,0.0,0.0,   0.0,      0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

/* simple Parameters trained using the small 500 molecules dataset.
const float kappa = 0.1960;

const float A1[] = {0.0,2.3594,0.0,0.0,0.0,0.0,2.4541,2.5908,2.7130,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.3833};
const float A2[] = {0.0,0.0,0.0,0.0,0.0,0.0,2.4726,2.5409,2.6766,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.4956};
const float A3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

const float B1[] = {0.0,0.5962,0.0,0.0,0.0,0.0,0.2591,0.3316,0.5028,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.4564};
const float B2[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.2268,0.2319,0.4992,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1493};
const float B3[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
*/

// function to retreive the atomtype value based on the "highest" (e.g. max) bond type of an atom
// potential improvement : in the original publication they don't have access to "Aromatic type" like in RDKit
double getAtomtype(const ROMol &mol, const RDKit::Atom *atom) {
      double t=1.0;
      double a;
      PRECONDITION(atom != nullptr,"bad atom argument")
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

double* getEEMMatrix(double * dist3D, int n, EEM_arrays EEMatoms) {
  int sizeArray = (n+1) * (n+1);
  PRECONDITION(dist3D != nullptr,"bad dist3D argument")
  double* EEM = new double[sizeArray](); // declaration to set all elements to zeros!
    /* Fill the full n * n block */
    for(long int i = 0; i < n; i++) {
        double t = EEMatoms.EEMatomtype[i];
        int idx = EEMatoms.Atomindex[i];
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

        EEM[i * (n+1) + i] = v;
        for(long int j = i + 1; j < n; j++) {
                EEM[i * (n+1) + j] = kappa / dist3D[i*n+j];
                EEM[j * (n+1) + i] = EEM[i * (n+1) + j];
            }
        }
    /* Fill last column & row */
    for(long int i = 0; i < n; i++) {
        EEM[n* (n+1) + i] = 1.0; // column
        EEM[i* (n+1) + n] = -1.0; // row

    }
    return EEM;
}

double* getBVector(const ROMol &mol, int n, EEM_arrays EEMatoms) {
       /* Fill vector b i.e. -A */
        double * b = new double[n+1];
        for(int j = 0; j < n; j++) {
          double t = EEMatoms.EEMatomtype[j];
          int idx = EEMatoms.Atomindex[j];
            if ( t== 1.0) {
                b[j] = -A1[idx];
            }
            if (t== 2.0) {
                b[j] = -A2[idx];
            }
            if (t== 3.0) {
                b[j] = -A3[idx];
            }
        }

        b[n] = MolOps::getFormalCharge(mol); // sum of charges

    return b;
}


EEM_arrays calculate_EEM_parameters(ROMol mol, int n) {
       /* Fill vector b i.e. -A */
        int * a = new int[n];
        double * b = new double[n];

        for(int j = 0; j < n; j++) {
            b[j]=getAtomtype(mol,mol.getAtomWithIdx(j));
            a[j]=mol.getAtomWithIdx(j)->getAtomicNum();
        }

    return EEM_arrays {a,b};
}


/* Calculate charges for a particular kappa_data structure */
void calculate_charges(ROMol mol, double * dist3D, int numAtoms, EEM_arrays EEMatoms, std::vector<double> &res) {
        double * A = getEEMMatrix(dist3D,numAtoms, EEMatoms);
        double * b = getBVector(mol, numAtoms, EEMatoms);

        MatrixXd AM = Map<MatrixXd>(A, numAtoms+1, numAtoms+1);
        VectorXd bv = Map<VectorXd>(b, numAtoms+1);
        VectorXd Res(numAtoms+1);

        FullPivLU<MatrixXd> lu(AM);
        Res = lu.solve(bv);

        std::vector<double>charges(Res.data(), Res.data() + Res.size());

        for(long int aix = 0; aix < numAtoms; aix++) {
            mol.getAtomWithIdx(aix)->setProp("_EEMCharge", charges[aix], true);
        }

        res =  charges;
       delete[] A;
       delete[] b;

}


void getEEMs(const ROMol &mol, std::vector<double> &result, int numAtoms, int confId) {
  // 3D distance matrix
  double * dist3D = MolOps::get3DDistanceMat(mol, confId, false, true);
  EEM_arrays EEMatoms = calculate_EEM_parameters(mol,  numAtoms);

  result.clear();
  result.resize(numAtoms);
  calculate_charges(mol, dist3D, numAtoms, EEMatoms, result);

}

}  // end of anonymous namespace

void EEM(ROMol &mol, std::vector<double> &res, int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  int numAtoms = mol.getNumAtoms();

  res.clear();
  res.resize(numAtoms);
  // cast ROMol to RWMol
  RWMol &wmol = static_cast<RWMol &>(mol);
  // kekulize is currenlty required but it could be remove if and only if:
  // we use "Aromatic type"  in RDKit retrain the model without Kekulize
  // that would be part of a future release if it's really important
  MolOps::Kekulize(wmol, true);
  // recast to ROMol
  ROMol *remol = static_cast<ROMol *>(&wmol);

  getEEMs(*remol, res, numAtoms, confId);
}
}  // end of Descriptors namespace
}  // end of RDKit namespace
