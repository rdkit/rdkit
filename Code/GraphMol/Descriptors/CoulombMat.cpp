//
//  Copyright (c) 2018, Guillaume GODIN
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
//  Modified by Greg Landrum, March 2020
//

#include <GraphMol/RDKitBase.h>
#include <cmath>
#include "CoulombMat.h"

namespace RDKit {
namespace Descriptors {

void CoulombMat(const ROMol &mol, std::vector<std::vector<double>> &res,
                int confId) {
  PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers");

  unsigned int numAtoms = mol.getNumAtoms();

  const auto conf = mol.getConformer(confId);

  res.resize(numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    res[i].resize(numAtoms);
    const auto at = mol.getAtomWithIdx(i);
    double Zi = at->getAtomicNum();
    res[i][i] = 0.5 * pow(Zi, 2.4);

    const auto Pi = conf.getAtomPos(i);
    for (unsigned int j = 0; j < i; ++j) {
      const auto Pj = conf.getAtomPos(j);
      double Zj = mol.getAtomWithIdx(j)->getAtomicNum();
      res[i][j] = res[j][i] = Zi * Zj / (Pi - Pj).length();
    }
  }
}
}  // namespace Descriptors
}  // namespace RDKit