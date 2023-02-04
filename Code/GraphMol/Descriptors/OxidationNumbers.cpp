//
//  Copyright (c) 2023, David Cosgrove, CozChemIx Limitied
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met,
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
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
// Calculate the oxidation numbers (states) of the atoms in a molecule.
// Based on the code at
// https,//github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/src/linchemin/cheminfo/functions.py#L544
// and therefore also subject to this licence from
// https,//github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/LICENSE
//
// MIT License
//
// Copyright (c) 2022 Syngenta Group Co. Ltd.
//
//    Permission is hereby granted, free of charge, to any person obtaining a
//    copy of this software and associated documentation files (the "Software"),
//    to deal in the Software without restriction, including without limitation
//    the rights to use, copy, modify, merge, publish, distribute, sublicense,
//    and/or sell
//                                                              copies of the
//                                                              Software, and to
//                                                              permit persons
//                                                              to whom the
//                                                              Software is
//    furnished to do so, subject to the following conditions,
//
//    The above copyright notice and this permission notice shall be included in
//    all copies or substantial portions of the Software.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//    DEALINGS IN THE SOFTWARE.

#include <GraphMol/RDKitBase.h>

namespace RDKit {
namespace Descriptors {

namespace {
int calcOxidationNumberByEN(const Atom *atom) {
  const static std::map<int, float> pauling_en_map = {
      {1, 2.2},   {3, 0.98},  {4, 1.57},  {5, 2.04},  {6, 2.55},  {7, 3.04},
      {8, 3.44},  {9, 3.98},  {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.9},
      {15, 2.19}, {16, 2.58}, {17, 3.16}, {19, 0.82}, {20, 1},    {21, 1.36},
      {22, 1.54}, {23, 1.63}, {24, 1.66}, {25, 1.55}, {26, 1.83}, {27, 1.88},
      {28, 1.91}, {29, 1.9},  {30, 1.65}, {31, 1.81}, {32, 2.01}, {33, 2.18},
      {34, 2.55}, {35, 2.96}, {36, 3.0},  {37, 0.82}, {38, 0.95}, {39, 1.22},
      {40, 1.33}, {41, 1.6},  {42, 2.16}, {43, 1.9},  {44, 2.2},  {45, 2.28},
      {46, 2.2},  {47, 1.93}, {48, 1.69}, {49, 1.78}, {50, 1.96}, {51, 2.05},
      {52, 2.1},  {53, 2.66}, {54, 2.6},  {55, 0.79}, {56, 0.89}, {57, 1.1},
      {58, 1.12}, {59, 1.13}, {60, 1.14}, {62, 1.17}, {64, 1.2},  {66, 1.22},
      {67, 1.23}, {68, 1.24}, {69, 1.25}, {71, 1.27}, {72, 1.3},  {73, 1.5},
      {74, 2.36}, {75, 1.9},  {76, 2.2},  {77, 2.2},  {78, 2.28}, {79, 2.54},
      {80, 2.0},  {81, 1.62}, {82, 2.33}, {83, 2.02}, {84, 2.0},  {85, 2.2},
      {88, 0.9},  {89, 1.1},  {90, 1.3},  {91, 1.5},  {92, 1.38}, {93, 1.36},
      {94, 1.28}, {95, 1.3},  {96, 1.3},  {97, 1.3},  {98, 1.3},  {99, 1.3},
      {100, 1.3}, {101, 1.3}, {102, 1.39}};
  PRECONDITION(atom->hasOwningMol(), "atom must have owning molecule")
  auto get_en = [&](int atomicNum) -> float {
    auto res = pauling_en_map.find(atomicNum);
    return res == pauling_en_map.end() ? 0.0 : res->second;
  };
  auto sf = [](float enDiff) -> int {
    if (enDiff > 0.0) {
      return -1;
    } else if (enDiff < 0.0) {
      return 1;
    } else {
      return 0;
    }
  };

  int oxNum = 0;
  float parEN = get_en(atom->getAtomicNum());

  for (const auto &bond : atom->getOwningMol().atomBonds(atom)) {
    auto otherAtom = bond->getOtherAtom(atom);
    if (otherAtom->getAtomicNum() > 1) {
      float en_diff = parEN - get_en(otherAtom->getAtomicNum());
      double bondType = bond->getBondTypeAsDouble();
      // Make sure this is a kekulized mol i.e. no bond type of 1.5.  This
      // shouldn't happen if called from calcOxidationNumbers, but who knows
      // what might happen in future.
      if (bondType > 1.0 && bondType < 2.0) {
        throw ValueErrorException(
            "Molecule appears not to be Kekulized,"
            " oxidation number calculation fails.");
      }
      oxNum += bondType * sf(en_diff);
    }
  }
  oxNum += sf(parEN - get_en(1)) * atom->getTotalNumHs();
  oxNum += atom->getFormalCharge();
  return oxNum;
}
}  // namespace

void calcOxidationNumbers(const ROMol &mol) {
  RWMol molCp(mol);
  RDKit::MolOps::Kekulize(molCp);
  for (const auto &atom : mol.atoms()) {
    auto cpAtom = molCp.getAtomWithIdx(atom->getIdx());
    int oxNum = calcOxidationNumberByEN(cpAtom);
    atom->setProp<int>(common_properties::OxidationNumber, oxNum);
  }
}

}  // end of namespace Descriptors
}  // end of namespace RDKit
