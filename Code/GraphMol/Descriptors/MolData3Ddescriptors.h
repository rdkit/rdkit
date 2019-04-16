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
//
#include <RDGeneral/export.h>
#ifndef MOLDATA3DDESCRIPTORS_2017
#define MOLDATA3DDESCRIPTORS_2017

#include <GraphMol/RDKitBase.h>
#include "Data3Ddescriptors.h"

class RDKIT_DESCRIPTORS_EXPORT MolData3Ddescriptors {
 private:
  Data3Ddescriptors data3D;

 public:
  MolData3Ddescriptors();
  std::vector<double> GetCharges(const RDKit::ROMol& mol);
  std::vector<double> GetRelativeMW(const RDKit::ROMol& mol);
  std::vector<double> GetCustomAtomProp(const RDKit::ROMol& mol,
                                        const std::string& customAtomPropName);
  std::vector<double> GetRelativePol(const RDKit::ROMol& mol);
  std::vector<double> GetRelativeRcov(const RDKit::ROMol& mol);
  std::vector<double> GetRelativeENeg(const RDKit::ROMol& mol);
  std::vector<double> GetRelativeIonPol(const RDKit::ROMol& mol);
  std::vector<double> GetRelativeVdW(const RDKit::ROMol& mol);
  std::vector<double> GetUn(int numAtoms);
  int GetPrincipalQuantumNumber(int AtomicNum);
  std::vector<double> GetIState(const RDKit::ROMol& mol);
  std::vector<double> GetIStateDrag(const RDKit::ROMol& mol);
  std::vector<double> GetEState(const RDKit::ROMol& mol);
  std::vector<double> GetEState2(const RDKit::ROMol& mol);
};
#endif
