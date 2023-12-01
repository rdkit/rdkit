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
// Adding RBF descriptors to 3D descriptors by Guillaume Godin

#include <GraphMol/RDKitBase.h>

#include "RDF.h"
#include "MolData3Ddescriptors.h"

#include <cmath>

namespace RDKit {
namespace Descriptors {
namespace {
MolData3Ddescriptors moldata3D;

std::vector<double> getG(int n) {
  std::vector<double> res(n);
  for (int i = 0; i < n; i++) {
    res[i] = 1 + i * 0.5;
  }
  return res;
}

std::vector<double> prepareIState(const ROMol& mol) {
  std::vector<double> IState = moldata3D.GetIStateDrag(
      mol);  // get the real IState value not the EState!
  return IState;
}

void getRDFDesc(double* DM, const ROMol& mol, const Conformer& conf,
                std::vector<double>& res) {
  // std::vector<double> reserror(210, 0); need to avoid the res NaNs
  // if (numAtoms < 4) return reserror;
  // if (!conf.is3D()) return reserror;

  int numAtoms = conf.getNumAtoms();

  std::vector<double> R = getG(30);
  std::vector<double> R1(210);
  std::vector<double> R2(30);
  std::vector<double> R3(30);
  std::vector<double> R4(30);
  std::vector<double> R5(30);
  std::vector<double> R6(30);
  std::vector<double> R7(30);

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
        p = exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
        res1 += p;                                                  // "u"
        res2 += Mass[j] * Mass[k] * p;                              // "m"
        res3 += RelativeVdW[j] * RelativeVdW[k] * p;                // "v"
        res4 += RelativeElectroNeg[j] * RelativeElectroNeg[k] * p;  //"e"
        res5 += RelativePol[j] * RelativePol[k] * p;                // "p"
        res6 += IonPol[j] * IonPol[k] * p;                          // "i"
        res7 += IState[j] * IState[k] * p;                          // "s"
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

  std::copy(R2.begin(), R2.end(), R1.begin() + 30);
  std::copy(R3.begin(), R3.end(), R1.begin() + 60);
  std::copy(R4.begin(), R4.end(), R1.begin() + 90);
  std::copy(R5.begin(), R5.end(), R1.begin() + 120);
  std::copy(R6.begin(), R6.end(), R1.begin() + 150);
  std::copy(R7.begin(), R7.end(), R1.begin() + 180);

  res = R1;
}

void getRDFDescCustom(double* DM, const ROMol& mol, const Conformer& conf,
                      std::vector<double>& res,
                      const std::string& customAtomPropName) {
  int numAtoms = conf.getNumAtoms();

  std::vector<double> R = getG(30);
  std::vector<double> R1(30);
  std::vector<double> customAtomArray =
      moldata3D.GetCustomAtomProp(mol, customAtomPropName);

  double p;
  for (size_t i = 0; i < R.size(); i++) {
    double res = 0.0;
    for (int j = 0; j < numAtoms - 1; j++) {
      for (int k = j + 1; k < numAtoms; k++) {
        p = exp(-100 * pow(R[i] - DM[j * numAtoms + k], 2));
        res += customAtomArray[j] * customAtomArray[k] * p;  // "custom"
      }
    }
    R1[i] = std::round(1000 * res) / 1000;
  }
  res = R1;
}

void GetRDF(double* dist3D, const ROMol& mol, const Conformer& conf,
            std::vector<double>& res) {
  getRDFDesc(dist3D, mol, conf, res);
}

void GetRDFone(double* dist3D, const ROMol& mol, const Conformer& conf,
               std::vector<double>& res, const std::string customAtomPropName) {
  getRDFDescCustom(dist3D, mol, conf, res, customAtomPropName);
}

}  // end of anonymous namespace

void RDF(const ROMol& mol, std::vector<double>& res, int confId,
         const std::string& customAtomPropName) {
  // RDF010u RDF015u RDF020u RDF025u RDF030u RDF035u RDF040u RDF045u RDF050u
  // RDF055u RDF060u RDF065u RDF070u RDF075u RDF080u RDF085u RDF090u RDF095u
  // RDF100u RDF105u RDF110u RDF115u RDF120u RDF125u RDF130u RDF135u RDF140u
  // RDF145u RDF150u RDF155u
  // RDF010m RDF015m RDF020m RDF025m RDF030m RDF035m RDF040m RDF045m RDF050m
  // RDF055m RDF060m RDF065m RDF070m RDF075m RDF080m RDF085m RDF090m RDF095m
  // RDF100m RDF105m RDF110m RDF115m RDF120m RDF125m RDF130m RDF135m RDF140m
  // RDF145m RDF150m RDF155m
  // RDF010v RDF015v RDF020v RDF025v RDF030v RDF035v RDF040v RDF045v RDF050v
  // RDF055v RDF060v RDF065v RDF070v RDF075v RDF080v RDF085v RDF090v RDF095v
  // RDF100v RDF105v RDF110v RDF115v RDF120v RDF125v RDF130v RDF135v RDF140v
  // RDF145v RDF150v RDF155v
  // RDF010e RDF015e RDF020e RDF025e RDF030e RDF035e RDF040e RDF045e RDF050e
  // RDF055e RDF060e RDF065e RDF070e RDF075e RDF080e RDF085e RDF090e RDF095e
  // RDF100e RDF105e RDF110e RDF115e RDF120e RDF125e RDF130e RDF135e RDF140e
  // RDF145e RDF150e RDF155e
  // RDF010p RDF015p RDF020p RDF025p RDF030p RDF035p RDF040p RDF045p RDF050p
  // RDF055p RDF060p RDF065p RDF070p RDF075p RDF080p RDF085p RDF090p RDF095p
  // RDF100p RDF105p RDF110p RDF115p RDF120p RDF125p RDF130p RDF135p RDF140p
  // RDF145p RDF150p RDF155p
  // RDF010i RDF015i RDF020i RDF025i RDF030i RDF035i RDF040i RDF045i RDF050i
  // RDF055i RDF060i RDF065i RDF070i RDF075i RDF080i RDF085i RDF090i RDF095i
  // RDF100i RDF105i RDF110i RDF115i RDF120i RDF125i RDF130i RDF135i RDF140i
  // RDF145i RDF150i RDF155i
  // RDF010s RDF015s RDF020s RDF025s RDF030s RDF035s RDF040s RDF045s RDF050s
  // RDF055s RDF060s RDF065s RDF070s RDF075s RDF080s RDF085s RDF090s RDF095s
  // RDF100s RDF105s RDF110s RDF115s RDF120s RDF125s RDF130s RDF135s RDF140s
  // RDF145s RDF150s RDF155s

  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
  // int numAtoms = mol.getNumAtoms();
  // if (numAtoms < 4) return reserror;

  const Conformer& conf = mol.getConformer(confId);
  // if (!conf.is3D()) return reserror;

  double* dist3D =
      MolOps::get3DDistanceMat(mol, confId, false, true);  // 3D distance matrix

  if (customAtomPropName != "") {
    // std::cout << " Using CustomAtomPropertie\n";
    // do something
    res.clear();
    res.resize(30);  // 7 * 30
    GetRDFone(dist3D, mol, conf, res, customAtomPropName);
  } else {
    res.clear();
    res.resize(210);  // 7 * 30
    GetRDF(dist3D, mol, conf, res);
  }
}
}  // namespace Descriptors
}  // namespace RDKit
