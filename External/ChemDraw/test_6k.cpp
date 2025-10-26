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
#include "chemdraw_doc.h"
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

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

#include <filesystem>

using namespace RDKit;
using namespace RDKit::v2;
namespace {
std::string replace(std::string &istr, const std::string &from,
                    const std::string &to) {
  std::string str(istr);
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos) {
    return str;
  }
  str.replace(start_pos, from.length(), to);
  return str;
}

bool hasNonSupportedFeatures(CDXDocument &document, const std::string &fname) {
  // check for monomers
  std::ifstream ifs(fname);
  std::stringstream xml;
  xml << ifs.rdbuf();
  // We should be able to figure this out from the node but...
  if (xml.str().find("monomerAttachmentStructure_") != std::string::npos ||
      xml.str().find("Name=\"monomerAttachments") != std::string::npos) {
    return true;
  }

  for (auto node : document.ContainedObjects()) {
    auto id = (CDXDatumID)node.second->GetTag();
    switch (id) {
      case kCDXObj_Page:
        for (auto frag : node.second->ContainedObjects()) {
          auto id = (CDXDatumID)frag.second->GetTag();
          if (id == kCDXObj_Fragment) {
            auto &fragment = (CDXFragment &)(*frag.second);
            if (fragment.m_sequenceType == kCDXSeqType_Unknown) {
              return true;
            }
          } else if (id == kCDXObj_BracketAttachment ||
                     id == kCDXObj_BracketedGroup) {
            return true;
          }
        }
        break;
      case kCDXObj_ObjectTag: {
        CDXObject &object = *((CDXObject *)node.second);
        id = (CDXDatumID)object.GetTag();
        // Check for monomers
        break;
      }
      default:
        break;
    }
  }
  return false;
}

bool hasNonSupportedFeatures(const std::string &fname) {
  auto doc = ChemDraw::ChemDrawToDocument(fname);
  return hasNonSupportedFeatures(*doc, fname);
}

TEST_CASE("Round TRIP") {
  std::string path =
      std::string(getenv("RDBASE")) + "/External/ChemDraw/test_data/CDXML6K/";

  SECTION("round trip") {
    // if we can't find the CDXML6K path, then don't run the test
    if (!std::filesystem::exists(path)) {
      return;
    }
    int failed = 0;
    int saniFailed = 0;
    int total = 0;
    int parseable = 0;
    int nomol = 0;
    int badparse = 0;
    int success = 0;
    int smimatches = 0;
    int nonSupported = 0;
    int no_mol_in_doc = 0;
    int bad_chemdraw_mol = 0;
    RDLog::LogStateSetter blocker;
    std::string cdxpath = path + "CDXML/";
    std::string molpath = path + "mol/";
    std::string smipath = path + "smiles/";

    std::string failpath = path + "FAILED/";
    std::string nomolpath = path + "NOMOL/";
    std::string badparsepath = path + "BADPARSE/";
    std::string sanitizationpath = path + "SANI/";

    std::set<std::string> known_failures{
        "INDMUMLL1117_2025-01-24-17-23-14_304.cdxml",  // Dative oxygen gets set
                                                       // to a radical
        "INDMUMLL1117_2025-01-24-17-26-06_1010.cdxml",  // The next batch has a
                                                        // type of stereochem I
                                                        // don't know how to
                                                        // parse yet
        "INDMUMLL1117_2025-01-24-17-26-06_1012.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1022.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1024.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1026.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1032.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1034.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1036.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1040.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1042.cdxml",
        "INDMUMLL1117_2025-01-24-17-26-06_1048.cdxml",  // Stereo chem batch
                                                        // ends here
        "INDMUMLL1117_2025-01-24-17-26-13_1690.cdxml",  // RDKit shows a radical
                                                        // for the dative ->[O]
        "INDMUMLL1117_2025-01-24-17-27-11_6877.cdxml",  // The next batch has a
                                                        // type of stereochem I
                                                        // don't know how to
                                                        // parse yet (same as
                                                        // before)
        "INDMUMLL1117_2025-01-24-17-27-11_6878.cdxml",
        "INDMUMLL1117_2025-01-24-17-27-11_6883.cdxml",
        "INDMUMLL1117_2025-01-24-17-27-11_6884.cdxml",
        "INDMUMLL1117_2025-01-24-17-27-11_6889.cdxml",
        "INDMUMLL1117_2025-01-24-17-27-11_6896.cdxml",
        "INDMUMLL1117_2025-01-24-17-27-30_8574.cdxml",   // Stereo chem batch
                                                         // ends here
        "INDMUMLL1117_2025-01-24-17-27-31_8633.cdxml",   // RDkit is missing a
                                                         // dummy atom molecule
        "INDMUMLL1117_2025-01-24-17-27-31_8651.cdxml",   // RDkit is missing a
                                                         // dummy atom molecule
        "INDMUMLL1117_2025-01-24-17-27-53_10330.cdxml",  // 2D projection of 3D
                                                         // stereo, we fail this
                                                         // one
        "INDMUMLL1117_2025-01-24-17-27-53_10332.cdxml",  // 2D projection of 3D
                                                         // stereo, we fail this
                                                         // one
        "INDMUMLL1117_2025-01-24-17-27-54_10336.cdxml",  // RDKit Smiles keeps
                                                         // any bonds ~,
                                                         // ChemDraw doesn't
        "INDMUMLL1117_2025-01-24-17-28-02_10942.cdxml",  // Chemdraw smiles
                                                         // doesn't support
                                                         // quadruple bond $
        "INDMUMLL1117_2025-01-24-17-28-15_11666.cdxml",  // RDKit Smiles keeps
                                                         // any bonds ~,
                                                         // ChemDraw doesn't
        "INDMUMLL1117_2025-01-24-17-28-20_12011.cdxml",  // RDKit gets stereo
                                                         // from the 3D data and
                                                         // the wedging
        "INDMUMLL1117_2025-01-24-17-28-20_12012.cdxml",  // RDKit gets stereo
                                                         // from the 3D data and
                                                         // the wedging
        "INDMUMLL1117_2025-01-24-17-28-21_12031.cdxml",  // 2D projection of 3D
                                                         // stereo, we fail this
                                                         // one
        "INDMUMLL1117_2025-01-24-17-28-30_12568.cdxml",  // 2D projection of 3D
                                                         // stereo, we fail this
                                                         // one
        "INDMUMLL1117_2025-01-24-17-29-06_14654.cdxml",  // Dative oxygen gets
                                                         // set to a radical
        "INDMUMLL1117_2025-01-24-17-29-08_14775.cdxml",  // RDKit Smiles keeps
                                                         // any bonds ~,
                                                         // ChemDraw doesn't
        "INDMUMLL1117_2025-01-24-17-29-09_14896.cdxml",  // We apparently do a
                                                         // bit of a better job
                                                         // than chemdraw here
                                                         // in parsing R/S
        "INDMUMLL1117_2025-01-24-17-29-09_14897.cdxml"   // RDKit just gets very
                                                         // different stereo
                                                         // chem, no idea why
    };

    for (auto p : {failpath, nomolpath, badparsepath, sanitizationpath}) {
      if (std::filesystem::exists(p)) {
        std::filesystem::remove_all(p);
      }
      std::filesystem::create_directory(p);
    }

    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(cdxpath)) {
      if (entry.is_regular_file()) {
        std::string fname = entry.path().filename().string();
        // issue here - graphite nanotube
        if (fname == "INDMUMLL1117_2025-01-24-17-28-02_10946.cdxml") {
          continue;  // nanotube takes forever
        }
        auto molfname = molpath + replace(fname, ".cdxml", ".mol");
        auto smifname = smipath + replace(fname, ".cdxml", ".smi");
        // if chemscript couldn't make an output, ignore it
        total++;

        if (!std::filesystem::exists(molfname) ||
            !std::filesystem::exists(smifname)) {
          no_mol_in_doc++;
          continue;
        }

        // Get the ChemScript mol and smiles
        std::unique_ptr<RWMol> mol;
        //= nullptr;
        try {
          mol.reset(MolFileToMol(molfname));
        } catch (...) {
          bad_chemdraw_mol++;
          continue;
        }
        // REQUIRE(mols.size());
        std::ifstream ifs(smifname);
        std::string smiles_in;
        ifs >> smiles_in;
        std::string smiles;
        {
          try {
            auto smimol = SmilesToMol(smiles_in);
            if (!smimol) {
              smiles = smiles_in;
            } else {
              smiles = MolToSmiles(*smimol);
              delete smimol;
            }
          } catch (...) {
            smiles = smiles_in;
          }
        }

        parseable++;
        // Read the cdxml
        std::vector<std::unique_ptr<RWMol>> mols;
        bool santizationFailure = false;
        try {
          mols = MolsFromChemDrawFile(entry.path().string());
          if (mols.size() == 0) {
            ChemDrawParserParams params;
            params.sanitize = false;
            mols = MolsFromChemDrawFile(entry.path().string(), params);
            santizationFailure = true;
          }
          if (!mols.size()) {
            if (smiles.size() == 0) {
              // At least we match the chemscript non-mol
              success++;
            } else if (hasNonSupportedFeatures(entry.path().string())) {
              // std::cerr << "[NOMOL (Unsupported)]: " << entry.path().string()
              //           << std::endl;
              nonSupported++;
            } else {
              std::cerr << "[NOMOL]: " << entry.path().string() << std::endl;
              std::filesystem::copy(
                  entry.path().string(),
                  nomolpath + entry.path().filename().string());
              nomol++;
            }
            continue;
          }
        } catch (...) {
          std::cerr << "[BADPARSE]: " << entry.path().string() << std::endl;
          std::filesystem::copy(
              entry.path(), badparsepath + entry.path().filename().string());
          badparse++;
          continue;
        }
        std::unique_ptr<ROMol> m = std::make_unique<ROMol>(*mols[0]);
        for (size_t i = 1; i < mols.size(); i++) {
          m.reset(combineMols(*m, *mols[i]));
        }

        auto rdkit_smi = MolToSmiles(*m);
        auto mol_smi = mol.get() ? MolToSmiles(*mol) : "";

        if (mol_smi != rdkit_smi) {
          // Do we match chemscripts smiles output at least?
          if (rdkit_smi == smiles) {
            smimatches++;
            continue;
          }

          if (hasNonSupportedFeatures(entry.path().string())) {
            nonSupported++;
            continue;  // has unsupported features
          }
          if (santizationFailure) {
            std::cerr << "[SANI]: " << entry.path() << std::endl;
            std::filesystem::copy(
                entry.path(),
                sanitizationpath + entry.path().filename().string());
            saniFailed++;
          } else {
            if (known_failures.find(entry.path().filename().string()) !=
                known_failures.end()) {
              continue;  // we know this failure and it's ok for now
            }

            std::cerr << "[FAIL]: " << entry.path() << std::endl;
            std::filesystem::copy(entry.path(),
                                  failpath + entry.path().filename().string());
            failed++;
          }
          std::cerr << "rdkit:               " << rdkit_smi << std::endl;
          std::cerr << "chemscript (mol):    " << mol_smi << std::endl;
          std::cerr << "chemscript (smiles): " << smiles << std::endl;
          std::cerr << molfname << std::endl;
          std::cerr << smifname << std::endl;
        } else {
          success++;
        }
      }
    }
    std::cerr << "Total:" << total << std::endl;
    std::cerr << "Parseable (has chemscript output):" << total << std::endl;
    std::cerr << "Success:" << success + smimatches << std::endl;
    std::cerr << "skipped (non supported features):" << nonSupported
              << std::endl;
    std::cerr << "skipped (no mol in doc):" << no_mol_in_doc << std::endl;
    std::cerr << "Chemscript smiles matches not chemscript mol: " << smimatches
              << std::endl;
    std::cerr << "Failed:" << failed << std::endl;
    std::cerr << "Sanitization:" << saniFailed << std::endl;
    std::cerr << "Nomol:" << nomol << std::endl;
    std::cerr << "Badparse:" << badparse << std::endl;
    std::cerr << "Bad ChemDraw Mol:" << bad_chemdraw_mol << std::endl;
    REQUIRE(failed == 0);
  }
}
}  // namespace
