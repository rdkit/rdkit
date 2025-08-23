/* 
* $Id$
*
*  Copyright (c) 2025, Glysade Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%{
#include <ChemDraw/chemdraw.h>
%}
%ignore RDKit::v2::MolsFromChemDrawDataStream;
%ignore RDKit::v2::MolsFromChemDrawFile;
%ignore RDKit::v2::MolsFromChemDrawBlock;
%rename("ChemDraw") RDKit::v2;
%rename(CDXFormat) RDKit::v2::CDXFormat;
%rename(ChemDrawParserParams) RDKit::v2::ChemDrawParserParams;

%include <ChemDraw/chemdraw.h>

%{
std::vector<RDKit::RWMOL_SPTR> MolsFromChemDrawBlockHelper(
 const std::string &text,
 const RDKit::v2::ChemDrawParserParams &params=RDKit::v2::ChemDrawParserParams()) {
  auto res = RDKit::v2::MolsFromChemDrawBlock(text, params);
  std::vector<RDKit::RWMOL_SPTR> mols;
  for(auto &mol: res) {
    mols.emplace_back(mol.release());
  }
  return mols;

}

std::vector<RDKit::RWMOL_SPTR> MolsFromChemDrawFileHelper(
 const std::string &filename,
 const RDKit::v2::ChemDrawParserParams &params=RDKit::v2::ChemDrawParserParams()) {
  auto res = RDKit::v2::MolsFromChemDrawFile(filename, params);
  std::vector<RDKit::RWMOL_SPTR> mols;
  for(auto &mol: res) {
    mols.emplace_back(mol.release());
  }
  return mols;
}

%}

%rename("MolsFromChemDrawBlock") MolsFromChemDrawBlockHelper;
%rename("MolsFromChemDrawFile") MolsFromChemDrawFileHelper;

std::vector<RDKit::RWMOL_SPTR> MolsFromChemDrawBlockHelper(
 const std::string &text,
 const RDKit::v2::ChemDrawParserParams &params=RDKit::v2::ChemDrawParserParams());

std::vector<RDKit::RWMOL_SPTR> MolsFromChemDrawFileHelper(
 const std::string &filename,
 const RDKit::v2::ChemDrawParserParams &params=RDKit::v2::ChemDrawParserParams());

