//
//  Copyright (c) 2023, David Cosgrove, CozChemIx Limitied
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
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/src/linchemin/cheminfo/functions.py#L544
// and therefore also subject to this licence from
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/LICENSE
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
//    furnished to do so, subject to the following conditions:
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

#include <RDGeneral/export.h>

/*!
  Calculates a molecule's exact molecular weight

\param mol        the molecule of interest
  \param onlyHeavy  (optional) if this is true (the default is false),
    only heavy atoms will be included in the MW calculation

*/

#ifndef __RD_OXIDATION_NUMBERS__
namespace RDKit {
class Atom;
class ROMol;
namespace Descriptors {

/*!
 * Calculates the oxidation numbers (states) of the atoms in a molecule
 * and stores them in the prop _OxidationNumber on the atoms.  Uses Pauling
 * electronegativies.
 *
 * @param mol the molecule of interest
 */
RDKIT_DESCRIPTORS_EXPORT void calculateOxidationNumbers(const ROMol &mol);

/*!
 * Calculate the oxidation number (state) of the atom of interest, which
 * must be in a parent molecule, preferably in kekulized form.
 * @param atom the atom of interest
 * @return the oxidation state as an integer
 */
RDKIT_DESCRIPTORS_EXPORT int calculateOxidationNumberByEN(const Atom *atom);

}  // end of namespace Descriptors
}  // end of namespace RDKit
#endif
