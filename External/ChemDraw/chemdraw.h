//
//  Copyright (c) 2024, Glysade Inc
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
#ifndef RDKIT_CHEMDRAW_H
#define RDKIT_CHEMDRAW_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <string>

namespace RDKit {
namespace v2 {
enum class CDXFormat {
  CDX = 1,
  CDXML = 2
};

struct RDKIT_RDCHEMDRAWLIB_EXPORT ChemDrawParserParams {
  bool sanitize;
  bool removeHs;
  CDXFormat format;
  ChemDrawParserParams() : sanitize(true), removeHs(true), format(CDXFormat::CDXML) {}
  ChemDrawParserParams(bool sanitize, bool removeHs, CDXFormat format) :
     sanitize(sanitize), removeHs(removeHs), format(format) {}
  
};

std::vector<std::unique_ptr<RWMol>> RDKIT_RDCHEMDRAWLIB_EXPORT
MolsFromChemDrawDataStream(std::istream &inStream,
               const ChemDrawParserParams &params = ChemDrawParserParams());

std::vector<std::unique_ptr<RWMol>> RDKIT_RDCHEMDRAWLIB_EXPORT
MolsFromChemDrawFile(const std::string &filename,
               const ChemDrawParserParams &params = ChemDrawParserParams());

std::vector<std::unique_ptr<RWMol>> RDKIT_RDCHEMDRAWLIB_EXPORT
MolsFromChemDrawBlock(const std::string &block,
               const ChemDrawParserParams &params = ChemDrawParserParams());

std::string RDKIT_RDCHEMDRAWLIB_EXPORT
MolToChemDrawBlock(const ROMol &mol, CDXFormat format = CDXFormat::CDXML);
}
}  // namespace RDKit
#endif
