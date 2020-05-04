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
// Created by Guillaume GODIN 2016

#include <GraphMol/RDKitBase.h>

#include "MORSE.h"
#include "MolData3Ddescriptors.h"

#include <cmath>

// data checked using book Todeschini R., Consonni V. - Molecular Descriptors
// for Chemoinformatics 2009 atomic properties page 21/22

namespace RDKit {
namespace Descriptors {
namespace {

MolData3Ddescriptors moldata3D;

std::vector<double> getG(int n) {
  std::vector<double> res(n);
  for (int i = 0; i < n; i++) {
    res[i] = i;
  }
  return res;
}

std::vector<double> prepareIState(const ROMol &mol) {
  std::vector<double> IState = moldata3D.GetIState(mol);
  return IState;
}

void getMORSEDesc(const double *DM, const ROMol &mol, const Conformer &conf,
                  std::vector<double> &res) {
  int numAtoms = conf.getNumAtoms();

  std::vector<double> R = getG(32);
  std::vector<double> R1(32);
  std::vector<double> R2(32);
  std::vector<double> R3(32);
  std::vector<double> R4(32);
  std::vector<double> R5(32);
  std::vector<double> R6(32);
  std::vector<double> R7(32);

  std::vector<double> Mass = moldata3D.GetRelativeMW(mol);
  std::vector<double> RelativePol = moldata3D.GetRelativePol(mol);
  std::vector<double> IonPol = moldata3D.GetRelativeIonPol(mol);
  std::vector<double> RelativeElectroNeg = moldata3D.GetRelativeENeg(mol);
  std::vector<double> RelativeVdW = moldata3D.GetRelativeVdW(mol);
  std::vector<double> IState = prepareIState(mol);

  double p;
  for (size_t i = 0; i < R.size(); i++) {
    double res1 = 0.0;
    double res2 = 0.0;
    double res3 = 0.0;
    double res4 = 0.0;
    double res5 = 0.0;
    double res6 = 0.0;
    double res7 = 0.0;

    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        if (i == 0) {
          p = 1;
        } else {
          p = sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
        }
        res1 += p;
        res2 += Mass[j] * Mass[k] * p;
        res3 += RelativeVdW[j] * RelativeVdW[k] * p;
        res4 += RelativeElectroNeg[j] * RelativeElectroNeg[k] * p;
        res5 += RelativePol[j] * RelativePol[k] * p;
        res6 += IonPol[j] * IonPol[k] * p;
        res7 += IState[j] * IState[k] * p;
      }
    }
    R1[i] = std::round(1000 * res1) / 1000;
    R2[i] = std::round(1000 * res2) / 1000;
    R3[i] = std::round(1000 * res3) / 1000;
    R4[i] = std::round(1000 * res4) / 1000;
    R5[i] = std::round(1000 * res5) / 1000;
    R6[i] = std::round(1000 * res6) / 1000;
    R7[i] = std::round(1000 * res7) / 1000;
  }

  R1.insert(R1.end(), R2.begin(), R2.end());
  R1.insert(R1.end(), R3.begin(), R3.end());
  R1.insert(R1.end(), R4.begin(), R4.end());
  R1.insert(R1.end(), R5.begin(), R5.end());
  R1.insert(R1.end(), R6.begin(), R6.end());
  R1.insert(R1.end(), R7.begin(), R7.end());

  res = R1;
}

void getMORSEDescCustom(const double *DM, const ROMol &mol,
                        const Conformer &conf, std::vector<double> &res,
                        const std::string &customAtomPropName) {
  int numAtoms = conf.getNumAtoms();

  std::vector<double> R = getG(32);
  std::vector<double> R1(32);
  std::vector<double> customAtomArray =
      moldata3D.GetCustomAtomProp(mol, customAtomPropName);

  double p;
  for (size_t i = 0; i < R.size(); i++) {
    double res1 = 0.0;

    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        if (i == 0) {
          p = 1;
        } else {
          p = sin(R[i] * DM[j * numAtoms + k]) / (R[i] * DM[j * numAtoms + k]);
        }
        res1 += customAtomArray[j] * customAtomArray[k] * p;  // "custom"
      }
    }
    R1[i] = std::round(1000 * res1) / 1000;
  }
  res = R1;
}

void GetMORSE(double *dist3D, const ROMol &mol, const Conformer &conf,
              std::vector<double> &res) {
  getMORSEDesc(dist3D, mol, conf, res);
}

void GetMORSEone(double *dist3D, const ROMol &mol, const Conformer &conf,
                 std::vector<double> &res,
                 const std::string &customAtomPropName) {
  getMORSEDescCustom(dist3D, mol, conf, res, customAtomPropName);
}

}  // end of anonymous namespace

void MORSE(const ROMol &mol, std::vector<double> &res, int confId,
           const std::string &customAtomPropName) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")

  // Mor01u Mor02u  Mor03u  Mor04u  Mor05u  Mor06u  Mor07u  Mor08u  Mor09u
  // Mor10u  Mor11u  Mor12u  Mor13u  Mor14u  Mor15u  Mor16u  Mor17u  Mor18u
  // Mor19u  Mor20u  Mor21u  Mor22u  Mor23u  Mor24u  Mor25u  Mor26u  Mor27u
  // Mor28u  Mor29u  Mor30u  Mor31u  Mor32u
  // Mor01m  Mor02m  Mor03m  Mor04m  Mor05m  Mor06m  Mor07m  Mor08m  Mor09m
  // Mor10m  Mor11m  Mor12m  Mor13m  Mor14m  Mor15m  Mor16m  Mor17m  Mor18m
  // Mor19m  Mor20m  Mor21m  Mor22m  Mor23m  Mor24m  Mor25m  Mor26m  Mor27m
  // Mor28m  Mor29m  Mor30m  Mor31m  Mor32m
  // Mor01v  Mor02v  Mor03v  Mor04v  Mor05v  Mor06v  Mor07v  Mor08v  Mor09v
  // Mor10v  Mor11v  Mor12v  Mor13v  Mor14v  Mor15v  Mor16v  Mor17v  Mor18v
  // Mor19v  Mor20v  Mor21v  Mor22v  Mor23v  Mor24v  Mor25v  Mor26v  Mor27v
  // Mor28v  Mor29v  Mor30v  Mor31v  Mor32v
  // Mor01e  Mor02e  Mor03e  Mor04e  Mor05e  Mor06e  Mor07e  Mor08e  Mor09e
  // Mor10e  Mor11e  Mor12e  Mor13e  Mor14e  Mor15e  Mor16e  Mor17e  Mor18e
  // Mor19e  Mor20e  Mor21e  Mor22e  Mor23e  Mor24e  Mor25e  Mor26e  Mor27e
  // Mor28e  Mor29e  Mor30e  Mor31e  Mor32e
  // Mor01p  Mor02p  Mor03p  Mor04p  Mor05p  Mor06p  Mor07p  Mor08p  Mor09p
  // Mor10p  Mor11p  Mor12p  Mor13p  Mor14p  Mor15p  Mor16p  Mor17p  Mor18p
  // Mor19p  Mor20p  Mor21p  Mor22p  Mor23p  Mor24p  Mor25p  Mor26p  Mor27p
  // Mor28p  Mor29p  Mor30p  Mor31p  Mor32p
  // Mor01i  Mor02i  Mor03i  Mor04i  Mor05i  Mor06i  Mor07i  Mor08i  Mor09i
  // Mor10i  Mor11i  Mor12i  Mor13i  Mor14i  Mor15i  Mor16i  Mor17i  Mor18i
  // Mor19i  Mor20i  Mor21i  Mor22i  Mor23i  Mor24i  Mor25i  Mor26i  Mor27i
  // Mor28i  Mor29i  Mor30i  Mor31i  Mor32i
  // Mor01s  Mor02s  Mor03s  Mor04s  Mor05s  Mor06s  Mor07s  Mor08s  Mor09s
  // Mor10s  Mor11s  Mor12s  Mor13s  Mor14s  Mor15s  Mor16s  Mor17s  Mor18s
  // Mor19s  Mor20s  Mor21s  Mor22s  Mor23s  Mor24s  Mor25s  Mor26s  Mor27s
  // Mor28s  Mor29s  Mor30s  Mor31s  Mor32s

  const Conformer &conf = mol.getConformer(confId);

  double *dist3D =
      MolOps::get3DDistanceMat(mol, confId, false, true);  // 3D distance matrix
  if (!customAtomPropName.empty()) {
    res.clear();
    res.resize(32);  // 7 * 32
    GetMORSEone(dist3D, mol, conf, res, customAtomPropName);

  } else {
    res.clear();
    res.resize(224);  // 7 * 32
    GetMORSE(dist3D, mol, conf, res);
  }
}

}  // namespace Descriptors
}  // namespace RDKit
