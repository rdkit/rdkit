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

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

namespace RDKit {
enum CDXFormat {
  CDX = 1,
  CDXML = 2
};
  
struct RDKIT_RDCHEMDRAWLIB_EXPORT ChemDrawParserParams {
  bool sanitize = true;
  bool removeHs = true;
  CDXFormat format = CDXML;
};

std::unique_ptr<CDXDocument> RDKIT_RDCHEMDRAWLIB_EXPORT
ChemDrawToDocument(std::istream &inStream, CDXFormat format);

std::unique_ptr<CDXDocument> RDKIT_RDCHEMDRAWLIB_EXPORT
ChemDrawToDocument(const std::string &filename);

std::vector<std::unique_ptr<RWMol>> RDKIT_RDCHEMDRAWLIB_EXPORT
ChemDrawToMols(std::istream &inStream,
               const ChemDrawParserParams &params = ChemDrawParserParams());

std::vector<std::unique_ptr<RWMol>> RDKIT_RDCHEMDRAWLIB_EXPORT
ChemDrawToMols(const std::string &filename,
               const ChemDrawParserParams &params = ChemDrawParserParams());

std::string RDKIT_RDCHEMDRAWLIB_EXPORT
MolToChemDraw(const ROMol &mol, CDXFormat format = CDXFormat::CDXML);

}  // namespace RDKit
#endif
