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
// Created by Guillaume Godin 2020

#include <RDGeneral/export.h>
#ifndef ATOMFEATRDKIT_H_MARC2020
#define ATOMFEATRDKIT_H_MARC2020

namespace RDKit {
class ROMol;
namespace Descriptors {
const std::string ATOMFEATVersion = "1.0.0";

std::vector <Atom::ChiralType> RS {  Atom::CHI_TETRAHEDRAL_CW, Atom::CHI_TETRAHEDRAL_CCW, Atom::CHI_OTHER };
std::vector <std::string> Symbols {"B", "C", "N", "O", "S", "F", "Si", "P", "Cl", "Br", "I", "H", "*"};
std::vector <Atom::HybridizationType> HS {  Atom::SP, Atom::SP2, Atom::SP3, Atom::SP3D,  Atom::SP3D2};


RDKIT_DESCRIPTORS_EXPORT void ATOMFEAT(
    const ROMol &, std::vector<double> &res, int atomid = 0,  bool addchiral = false );
}  // namespace Descriptors
}  // namespace RDKit
#endif
