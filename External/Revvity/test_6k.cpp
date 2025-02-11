//
//  Copyright (c) 2025 Glysade Inc and other RDkit contributors
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

#include "chemdraw.h"
#include <catch2/catch_all.hpp>
#include "RDGeneral/test.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>
#include <CDXStdObjects.h>

using namespace RDKit;
namespace {
std::string replace(std::string& istr, const std::string& from, const std::string& to) {
  std::string str(istr);
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return str;
  str.replace(start_pos, from.length(), to);
  return str;
}

bool hasNonPolymerFragment(CDXDocument &document) {
  for (auto node : document.ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)node.second->GetTag();
    switch (id) {
      case kCDXObj_Page:
        for (auto frag : node.second->ContainedObjects()) {
          CDXDatumID id = (CDXDatumID)frag.second->GetTag();
          if (id == kCDXObj_Fragment) {
            CDXFragment &fragment = (CDXFragment &)(*frag.second);
            if(fragment.m_sequenceType == kCDXSeqType_Unknown)
              return true;
          }
          }
        }
    }
  return false;
}

bool hasNonPolymerFragment(const std::string &fname) {
  auto doc = ChemDrawToDocument(fname);
  return hasNonPolymerFragment(*doc);
}

TEST_CASE("Round TRIP") {
  std::string path = std::string(getenv("RDBASE")) + "External/Revvity/test_data/CDXML6K/";

  SECTION("round trip") {
    int failed = 0;
    int saniFailed = 0;
    int nomol = 0;
    int badparse = 0;
    int total = 0;
    int success = 0;
    int smimatches = 0;
    int skippedBiopolymer = 0;
    RDLog::LogStateSetter blocker;
    std::string cdxpath = path + "CDXML/";
    std::string molpath = path + "mol/";
    std::string smipath = path + "smiles/";
    
    std::string failpath = path + "FAILED/";
    std::string nomolpath = path + "NOMOL/";
    std::string badparsepath = path + "BADPARSE/";
    std::string sanitizationpath = path + "SANI/";
    for (auto p : {failpath, nomolpath, badparsepath, sanitizationpath}) {
      if(std::filesystem::exists(p)) {
        std::filesystem::remove_all(p);
      }
      std::filesystem::create_directory(p);
    }
    
    for (const auto & entry : std::filesystem::recursive_directory_iterator(cdxpath)) {
      if (entry.is_regular_file()) {
        std::string fname = entry.path().filename();
        //if (fname != "INDMUMLL1117_2025-01-24-17-25-54_22.cdxml") {
        //  continue;
        //}
        //if(fname != "INDMUMLL1117_2025-01-24-17-28-09_11375.cdxml") continue;
        //if(fname != "INDMUMLL1117_2024-04-12-15-02-47_26.cdxml") continue;
        if(fname != "INDMUMLL1117_2025-01-24-17-23-13_257.cdxml") continue;
        
        auto molfname = molpath + replace(fname, ".cdxml", ".mol");
        auto smifname = smipath + replace(fname, ".cdxml", ".smi");
        if (!std::filesystem::exists(molfname) || !std::filesystem::exists(smifname)) {
          continue;
        }
        total++;
        // Read the cdxml
        std::vector<std::unique_ptr<ROMol>> mols;
        bool santizationFailure = false;
        try {
          mols = ChemDrawToMols(entry.path());
          if(mols.size() == 0) {
            ChemDrawParserParams params;
            params.sanitize = false;
            mols = ChemDrawToMols(entry.path(), params);
            santizationFailure = true;
          }
          if(!mols.size()) {
            if(hasNonPolymerFragment(entry.path())) {
              std::cerr << "[NOMOL]: " << entry.path() << std::endl;
              std::filesystem::copy(entry.path(), nomolpath + (std::string)entry.path().filename());
              nomol++;
            } else {
              skippedBiopolymer++;
            }
            continue;
          }
        } catch (...) {
            std::cerr << "[BADPARSE]: " << entry.path() << std::endl;
            std::filesystem::copy(entry.path(), badparsepath + (std::string)entry.path().filename());
            badparse++;
            continue;
        }
        RWMol *mol = nullptr;
        try {
          mol = MolFileToMol(molfname);
        } catch (...) {
          continue;
        }
        //REQUIRE(mols.size());
        std::ifstream  ifs(smifname);
        std::string smiles_in;
        ifs >> smiles_in;
        std::string smiles;
        {
          try {
            auto smimol = SmilesToMol(smiles_in);
            if(!smimol) smiles = smiles_in;
            else {
              smiles = MolToSmiles(*smimol);
              delete smimol;
            }
          } catch (...) {
            smiles = smiles_in;
          }
        }
        std::unique_ptr<ROMol> m = std::make_unique<ROMol>(*mols[0]);
        for(size_t i=1; i<mols.size(); i++) {
          m.reset(combineMols(*m, *mols[i]));
        }
        
        auto rdkit_smi = MolToSmiles(*m);
        auto mol_smi = MolToSmiles(*mol);
        delete mol;
        if(mol_smi != rdkit_smi) {
          // Do we match chemscripts smiles output at least?
          if(rdkit_smi == smiles) {
            smimatches++;
            continue;
          }
          
          if(!hasNonPolymerFragment(entry.path())) {
            continue; // biopolymers not supported yet
          }
          if (santizationFailure)
          {
            std::cerr << "[SANI]: " << entry.path() << std::endl;
            std::filesystem::copy(entry.path(), sanitizationpath + (std::string)entry.path().filename());
            saniFailed++;
          }
          else {
            std::cerr << "[FAIL]: " << entry.path() << std::endl;
            std::filesystem::copy(entry.path(), failpath + (std::string)entry.path().filename());
            failed++;
          }
          std::cerr << "rdkit:               " << rdkit_smi << std::endl;
          std::cerr << "chemscript (mol):    " << mol_smi << std::endl;
          std::cerr << "chemscript (smiles): " << smiles << std::endl;
          std::cerr << molfname << std::endl;
          std::cerr << smifname << std::endl;
        }
        else {
          success++;
        }
      }
    }
    std::cerr << "Success:" << success + smimatches << std::endl;
    std::cerr << "Chemscript smiles matches not chemscript mol: " << smimatches << std::endl;
    std::cerr << "Failed:" << failed << std::endl;
    std::cerr << "Sanitization:" << saniFailed << std::endl;
    std::cerr << "Nomol:" << nomol << std::endl;
    std::cerr << "Badparse:" << badparse <<  std::endl;
    REQUIRE(failed == 0);
  }
}
}
