//
//  Copyright (c) 2010-2013, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/MolHash/MolHash.h>
#include <GraphMol/FMCS/FMCS.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/integer_traits.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <RDGeneral/BoostEndInclude.h>

#ifdef BUILD_INCHI_SUPPORT
#include <INCHI-API/inchi.h>
#endif
#ifdef BUILD_AVALON_SUPPORT
#include <AvalonTools/AvalonTools.h>
#endif
#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

#include "rdkit.h"
#include "guc.h"
#include "bitstring.h"

using namespace std;
using namespace RDKit;

class ByteA : public std::string {
 public:
  ByteA() : string(){};
  ByteA(bytea *b) : string(VARDATA(b), VARSIZE(b) - VARHDRSZ){};
  ByteA(string &s) : string(s){};

  /*
   * Convert string to bytea. Convertaion is in pgsql's memory
   */
  bytea *toByteA() {
    bytea *res;
    int len;

    len = this->size();
    res = (bytea *)palloc(VARHDRSZ + len);
    memcpy(VARDATA(res), this->data(), len);
    SET_VARSIZE(res, VARHDRSZ + len);

    return res;
  };

  /* Just the copy of string's method */
  ByteA &operator=(const string &__str) {
    return (ByteA &)this->assign(__str);
  };
};

/*
 * Constant io
 */
static string StringData;

/*
 * Real sparse vector
 */

typedef SparseIntVect<boost::uint32_t> SparseFP;

/*******************************************
 *        ROMol transformation             *
 *******************************************/

extern "C" void freeCROMol(CROMol data) {
  ROMol *mol = (ROMol *)data;
  delete mol;
}

extern "C" CROMol constructROMol(Mol *data) {
  ROMol *mol = new ROMol();

  try {
    ByteA b(data);
    MolPickler::molFromPickle(b, mol);
  } catch (MolPicklerException &e) {
    elog(ERROR, "molFromPickle: %s", e.message());
  } catch (...) {
    elog(ERROR, "constructROMol: Unknown exception");
  }

  return (CROMol)mol;
}

extern "C" Mol *deconstructROMol(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ByteA b;

  try {
    MolPickler::pickleMol(mol, b);
  } catch (MolPicklerException &e) {
    elog(ERROR, "pickleMol: %s", e.message());
  } catch (...) {
    elog(ERROR, "deconstructROMol: Unknown exception");
  }

  return (Mol *)b.toByteA();
}

extern "C" CROMol parseMolText(char *data, bool asSmarts, bool warnOnFail,
                               bool asQuery) {
  RWMol *mol = NULL;

  try {
    if (!asSmarts) {
      if (!asQuery) {
        mol = SmilesToMol(data);
      } else {
        mol = SmilesToMol(data, 0, false);
        MolOps::sanitizeMol(*mol);
        MolOps::mergeQueryHs(*mol);
      }
    } else {
      mol = SmartsToMol(data, 0, false);
    }
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create molecule from SMILES '%s'", data)));
    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create molecule from SMILES '%s'", data)));
    }
  }

  return (CROMol)mol;
}

extern "C" CROMol parseMolBlob(char *data, int len) {
  ROMol *mol = NULL;

  try {
    string binStr(data, len);
    mol = new ROMol(binStr);
  } catch (...) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("problem generating molecule from blob data")));
  }
  if (mol == NULL) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("blob data could not be parsed")));
  }

  return (CROMol)mol;
}

extern "C" CROMol parseMolCTAB(char *data, bool keepConformer, bool warnOnFail,
                               bool asQuery) {
  RWMol *mol = NULL;

  try {
    if (!asQuery) {
      mol = MolBlockToMol(data);
    } else {
      mol = MolBlockToMol(data, true, false);
      MolOps::mergeQueryHs(*mol);
    }
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create molecule from CTAB '%s'", data)));

    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create molecule from CTAB '%s'", data)));
    }
  } else {
    if (!keepConformer) mol->clearConformers();
  }

  return (CROMol)mol;
}

extern "C" bool isValidSmiles(char *data) {
  RWMol *mol = NULL;
  bool res;
  try {
    string str(data);
    if (str.empty()) {
      // Pass the test - No-Structure input is allowed. No cleanup necessary.
      return true;
    }
    mol = SmilesToMol(str, 0, 0);
    if (mol) {
      MolOps::cleanUp(*mol);
      mol->updatePropertyCache();
      MolOps::Kekulize(*mol);
      MolOps::assignRadicals(*mol);
      MolOps::setAromaticity(*mol);
      MolOps::adjustHs(*mol);
    }
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    res = false;
  } else {
    res = true;
    delete mol;
  }
  return res;
}

extern "C" bool isValidSmarts(char *data) {
  ROMol *mol = NULL;
  bool res;
  try {
    mol = SmartsToMol(data);
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    res = false;
  } else {
    res = true;
    delete mol;
  }
  return res;
}

extern "C" bool isValidCTAB(char *data) {
  RWMol *mol = NULL;
  bool res;
  try {
    mol = MolBlockToMol(data, false, false);
    if (mol) {
      MolOps::cleanUp(*mol);
      mol->updatePropertyCache();
      MolOps::Kekulize(*mol);
      MolOps::assignRadicals(*mol);
      MolOps::setAromaticity(*mol);
      MolOps::adjustHs(*mol);
    }
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    res = false;
  } else {
    res = true;
    delete mol;
  }
  return res;
}

extern "C" bool isValidMolBlob(char *data, int len) {
  ROMol *mol = NULL;
  bool res = false;
  try {
    string binStr(data, len);
    mol = new ROMol(binStr);
  } catch (...) {
    mol = NULL;
  }
  if (mol == NULL) {
    res = false;
  } else {
    delete mol;
    res = true;
  }
  return res;
}

extern "C" char *makeMolText(CROMol data, int *len, bool asSmarts) {
  ROMol *mol = (ROMol *)data;

  try {
    if (!asSmarts) {
      StringData = MolToSmiles(*mol, true);
    } else {
      StringData = MolToSmarts(*mol, false);
    }
  } catch (...) {
    ereport(
        WARNING,
        (errcode(ERRCODE_WARNING),
         errmsg("makeMolText: problems converting molecule to SMILES/SMARTS")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" char *makeCtabText(CROMol data, int *len,
                              bool createDepictionIfMissing) {
  ROMol *mol = (ROMol *)data;

  try {
    if (createDepictionIfMissing && mol->getNumConformers() == 0) {
      RDDepict::compute2DCoords(*mol);
    }
    StringData = MolToMolBlock(*mol);
  } catch (...) {
    ereport(WARNING,
            (errcode(ERRCODE_WARNING),
             errmsg("makeCtabText: problems converting molecule to CTAB")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" char *makeMolBlob(CROMol data, int *len) {
  ROMol *mol = (ROMol *)data;
  StringData.clear();
  try {
    MolPickler::pickleMol(*mol, StringData);
  } catch (...) {
    elog(ERROR, "makeMolBlob: Unknown exception");
  }

  *len = StringData.size();
  return (char *)StringData.data();
}

extern "C" bytea *makeMolSignature(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  bytea *ret = NULL;

  try {
    res = RDKit::PatternFingerprintMol(*mol, getSubstructFpSize());
    // res =
    // RDKit::LayeredFingerprintMol(*mol,RDKit::substructLayers,1,5,SSS_FP_SIZE);

    if (res) {
      std::string sres = BitVectToBinaryText(*res);
      
      unsigned int varsize = VARHDRSZ + sres.size();
      ret = (bytea *) palloc0(varsize);
      memcpy(VARDATA(ret), sres.data(), sres.size());
      SET_VARSIZE(ret, varsize);
      
      delete res;
      res = 0;
    }
  } catch (...) {
    elog(ERROR, "makeMolSignature: Unknown exception");
    if (res) delete res;
  }

  return ret;
}

extern "C" int molcmp(CROMol i, CROMol a) {
  ROMol *im = (ROMol *)i;
  ROMol *am = (ROMol *)a;

  if (!im) {
    if (!am) return 0;
    return -1;
  }
  if (!am) return 1;

  int res = im->getNumAtoms() - am->getNumAtoms();
  if (res) return res;

  res = im->getNumBonds() - am->getNumBonds();
  if (res) return res;

  res = int(RDKit::Descriptors::calcAMW(*im, false)) -
        int(RDKit::Descriptors::calcAMW(*am, false));
  if (res) return res;

  res = im->getRingInfo()->numRings() - am->getRingInfo()->numRings();
  if (res) return res;

  RDKit::MatchVectType matchVect;
  bool ss1 = RDKit::SubstructMatch(*im, *am, matchVect, false, false);
  bool ss2 = RDKit::SubstructMatch(*am, *im, matchVect, false, false);
  if (ss1 != ss2) return ss1;

  // the above can still fail in some chirality cases
  std::string smi1 = MolToSmiles(*im, true);
  std::string smi2 = MolToSmiles(*am, true);
  return smi1 == smi2 ? 0 : (smi1 < smi2 ? -1 : 1);
}

extern "C" int MolSubstruct(CROMol i, CROMol a) {
  ROMol *im = (ROMol *)i;
  ROMol *am = (ROMol *)a;
  RDKit::MatchVectType matchVect;

  return RDKit::SubstructMatch(*im, *am, matchVect, true, getDoChiralSSS());
}

extern "C" int MolSubstructCount(CROMol i, CROMol a, bool uniquify) {
  ROMol *im = (ROMol *)i;
  ROMol *am = (ROMol *)a;
  std::vector<RDKit::MatchVectType> matchVect;

  return static_cast<int>(RDKit::SubstructMatch(*im, *am, matchVect, uniquify,
                                                true, getDoChiralSSS()));
}

/*******************************************
 *     Molecule operations                 *
 *******************************************/
#define MOLDESCR(name, func, ret)      \
  extern "C" ret Mol##name(CROMol i) { \
    const ROMol *im = (ROMol *)i;      \
    return func(*im);                  \
  }
MOLDESCR(FractionCSP3, RDKit::Descriptors::calcFractionCSP3, double)
MOLDESCR(TPSA, RDKit::Descriptors::calcTPSA, double)
MOLDESCR(AMW, RDKit::Descriptors::calcAMW, double)
MOLDESCR(HBA, RDKit::Descriptors::calcLipinskiHBA, int)
MOLDESCR(HBD, RDKit::Descriptors::calcLipinskiHBD, int)
MOLDESCR(NumHeteroatoms, RDKit::Descriptors::calcNumHeteroatoms, int)
MOLDESCR(NumRings, RDKit::Descriptors::calcNumRings, int)
MOLDESCR(NumAromaticRings, RDKit::Descriptors::calcNumAromaticRings, int)
MOLDESCR(NumAliphaticRings, RDKit::Descriptors::calcNumAliphaticRings, int)
MOLDESCR(NumSaturatedRings, RDKit::Descriptors::calcNumSaturatedRings, int)
MOLDESCR(NumAromaticHeterocycles,
         RDKit::Descriptors::calcNumAromaticHeterocycles, int)
MOLDESCR(NumAliphaticHeterocycles,
         RDKit::Descriptors::calcNumAliphaticHeterocycles, int)
MOLDESCR(NumSaturatedHeterocycles,
         RDKit::Descriptors::calcNumSaturatedHeterocycles, int)
MOLDESCR(NumAromaticCarbocycles, RDKit::Descriptors::calcNumAromaticCarbocycles,
         int)
MOLDESCR(NumAliphaticCarbocycles,
         RDKit::Descriptors::calcNumAliphaticCarbocycles, int)
MOLDESCR(NumSaturatedCarbocycles,
         RDKit::Descriptors::calcNumSaturatedCarbocycles, int)
MOLDESCR(NumHeterocycles, RDKit::Descriptors::calcNumHeterocycles, int)

MOLDESCR(NumRotatableBonds, RDKit::Descriptors::calcNumRotatableBonds, int)
MOLDESCR(Chi0v, RDKit::Descriptors::calcChi0v, double)
MOLDESCR(Chi1v, RDKit::Descriptors::calcChi1v, double)
MOLDESCR(Chi2v, RDKit::Descriptors::calcChi2v, double)
MOLDESCR(Chi3v, RDKit::Descriptors::calcChi3v, double)
MOLDESCR(Chi4v, RDKit::Descriptors::calcChi4v, double)
MOLDESCR(Chi0n, RDKit::Descriptors::calcChi0n, double)
MOLDESCR(Chi1n, RDKit::Descriptors::calcChi1n, double)
MOLDESCR(Chi2n, RDKit::Descriptors::calcChi2n, double)
MOLDESCR(Chi3n, RDKit::Descriptors::calcChi3n, double)
MOLDESCR(Chi4n, RDKit::Descriptors::calcChi4n, double)
MOLDESCR(Kappa1, RDKit::Descriptors::calcKappa1, double)
MOLDESCR(Kappa2, RDKit::Descriptors::calcKappa2, double)
MOLDESCR(Kappa3, RDKit::Descriptors::calcKappa3, double)
MOLDESCR(NumSpiroAtoms, RDKit::Descriptors::calcNumSpiroAtoms, int)
MOLDESCR(NumBridgeheadAtoms, RDKit::Descriptors::calcNumBridgeheadAtoms, int)

extern "C" double MolLogP(CROMol i) {
  double logp, mr;
  RDKit::Descriptors::calcCrippenDescriptors(*(ROMol *)i, logp, mr);
  return logp;
}
extern "C" int MolNumAtoms(CROMol i) {
  const ROMol *im = (ROMol *)i;
  return im->getNumAtoms(false);
}
extern "C" int MolNumHeavyAtoms(CROMol i) {
  const ROMol *im = (ROMol *)i;
  return im->getNumHeavyAtoms();
}

extern "C" char *makeMolFormulaText(CROMol data, int *len,
                                    bool separateIsotopes,
                                    bool abbreviateHIsotopes) {
  ROMol *mol = (ROMol *)data;

  try {
    StringData = RDKit::Descriptors::calcMolFormula(*mol, separateIsotopes,
                                                    abbreviateHIsotopes);
  } catch (...) {
    ereport(WARNING,
            (errcode(ERRCODE_WARNING),
             errmsg("makeMolFormulaText: problems converting molecule to "
                    "sum formula")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" const char *MolInchi(CROMol i, const char *opts) {
  std::string inchi = "InChI not available";
#ifdef BUILD_INCHI_SUPPORT
  const ROMol *im = (ROMol *)i;
  ExtraInchiReturnValues rv;
  try {
    std::string sopts = "/AuxNone /WarnOnEmptyStructure";
    if (strlen(opts)) {
      sopts += std::string(" ") + std::string(opts);
    }
    inchi = MolToInchi(*im, rv, sopts.c_str());
  } catch (MolSanitizeException &e) {
    inchi = "";
    elog(ERROR, "MolInchi: cannot kekulize molecule");
  } catch (...) {
    inchi = "";
    elog(ERROR, "MolInchi: Unknown exception");
  }
#endif
  return strdup(inchi.c_str());
}
extern "C" const char *MolInchiKey(CROMol i, const char *opts) {
  std::string key = "InChI not available";
#ifdef BUILD_INCHI_SUPPORT
  const ROMol *im = (ROMol *)i;
  ExtraInchiReturnValues rv;
  try {
    std::string sopts = "/AuxNone /WarnOnEmptyStructure";
    if (strlen(opts)) {
      sopts += std::string(" ") + std::string(opts);
    }
    std::string inchi = MolToInchi(*im, rv, sopts.c_str());
    key = InchiToInchiKey(inchi);
  } catch (MolSanitizeException &e) {
    key = "";
    elog(ERROR, "MolInchiKey: cannot kekulize molecule");
  } catch (...) {
    key = "";
    elog(ERROR, "MolInchiKey: Unknown exception");
  }
#endif
  return strdup(key.c_str());
}

extern "C" CROMol MolMurckoScaffold(CROMol i) {
  const ROMol *im = (ROMol *)i;
  ROMol *mol = MurckoDecompose(*im);
  if (mol && !mol->getNumAtoms()) {
    delete mol;
    mol = 0;
  } else {
    try {
      MolOps::sanitizeMol(*(RWMol *)mol);
    } catch (...) {
      delete mol;
      mol = 0;
    }
  }
  return (CROMol)mol;
}

namespace {
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

unsigned int parseWhichString(const std::string &txt) {
  unsigned int res = MolOps::ADJUST_IGNORENONE;
  boost::char_separator<char> sep("|");
  tokenizer tokens(txt, sep);
  tokenizer::iterator token = tokens.begin();
  while (token != tokens.end()) {
    std::string v = *token;
    ++token;
    if (v == "IGNORENONE") {
      res |= MolOps::ADJUST_IGNORENONE;
    } else if (v == "IGNORERINGATOMS") {
      res |= MolOps::ADJUST_IGNORERINGATOMS;
    } else if (v == "IGNORECHAINATOMS") {
      res |= MolOps::ADJUST_IGNORECHAINATOMS;
    } else if (v == "IGNOREDUMMIES") {
      res |= MolOps::ADJUST_IGNOREDUMMIES;
    } else if (v == "IGNORENONDUMMIES") {
      res |= MolOps::ADJUST_IGNORENONDUMMIES;
    } else if (v == "IGNOREALL") {
      res |= MolOps::ADJUST_IGNOREALL;
    } else {
      elog(ERROR, "bad which string provided '%s'", v.c_str());
    }
  }
  return res;
}

void parseAdjustQueryParameters(MolOps::AdjustQueryParameters &p,
                                const char *json) {
  PRECONDITION(json && strlen(json), "empty json");
  std::istringstream ss;
  ss.str(json);

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  p.adjustDegree = pt.get("adjustDegree", p.adjustDegree);
  p.adjustRingCount = pt.get("adjustRingCount", p.adjustRingCount);
  p.makeDummiesQueries = pt.get("makeDummiesQueries", p.makeDummiesQueries);
  std::string which;
  which = boost::to_upper_copy<std::string>(pt.get("adjustDegreeFlags", ""));
  if (which != "") p.adjustDegreeFlags = parseWhichString(which);
  which = boost::to_upper_copy<std::string>(pt.get("adjustRingCountFlags", ""));
  if (which != "") p.adjustRingCountFlags = parseWhichString(which);
}
}

extern "C" CROMol MolAdjustQueryProperties(CROMol i, const char *params) {
  const ROMol *im = (ROMol *)i;

  MolOps::AdjustQueryParameters p;

  if (params && strlen(params)) {
    // try {
    parseAdjustQueryParameters(p, params);
    // } catch (...) {
    //   ereport(
    //       WARNING,
    //       (errcode(ERRCODE_WARNING),
    //        errmsg(
    //            "adjustQueryProperties: Invalid argument \'params\'
    //            ignored")));
    // }
  }
  ROMol *mol = MolOps::adjustQueryProperties(*im, &p);
  return (CROMol)mol;
}

/*******************************************
 *     CBfp transformation                 *
 *******************************************/

extern "C" void freeCBfp(CBfp data) {
  std::string *fp = (std::string *)data;
  delete fp;
}

extern "C" CBfp constructCBfp(Bfp *data) {
  std::string *ebv = NULL;

  try {
    ebv = new std::string(VARDATA(data), VARSIZE(data) - VARHDRSZ);
  } catch (...) {
    elog(ERROR, "constructMolFingerPrint: Unknown exception");
  }

  return (CBfp)ebv;
}

extern "C" Bfp *deconstructCBfp(CBfp data) {
  std::string *ebv = (std::string *)data;
  ByteA b;

  try {
    b = *ebv;
  } catch (...) {
    elog(ERROR, "deconstructMolFingerPrint: Unknown exception");
  }

  return b.toByteA();
}

extern "C" BfpSignature *makeBfpSignature(CBfp data) {
  std::string *ebv = (std::string *)data;
  int siglen = ebv->size();
  
  unsigned int varsize = sizeof(BfpSignature) + siglen;
  BfpSignature * res = (BfpSignature *)palloc0(varsize);
  SET_VARSIZE(res, varsize);

  res->weight = bitstringWeight(siglen, (uint8 *)ebv->data());
  memcpy(res->fp, ebv->data(), siglen);
  
  return res;
}

extern "C" int CBfpSize(CBfp a) {
  std::string *ebv = (std::string *)a;
  int numBits = ebv->size() * 8;
  return numBits;
}

extern "C" double calcBitmapTanimotoSml(CBfp a, CBfp b) {
  std::string *abv = (std::string *)a;
  std::string *bbv = (std::string *)b;
  const unsigned char *afp = (const unsigned char *)abv->c_str();
  const unsigned char *bfp = (const unsigned char *)bbv->c_str();
  /* return CalcBitmapTanimoto(afp, bfp, abv->size()); */
  return bitstringTanimotoSimilarity(abv->size(), (uint8 *)afp, (uint8 *)bfp);
}

extern "C" double calcBitmapDiceSml(CBfp a, CBfp b) {
  std::string *abv = (std::string *)a;
  std::string *bbv = (std::string *)b;
  const unsigned char *afp = (const unsigned char *)abv->c_str();
  const unsigned char *bfp = (const unsigned char *)bbv->c_str();
  return CalcBitmapDice(afp, bfp, abv->size());
}

double calcBitmapTverskySml(CBfp a, CBfp b, float ca, float cb) {
  std::string *abv = (std::string *)a;
  std::string *bbv = (std::string *)b;
  const unsigned char *afp = (const unsigned char *)abv->c_str();
  const unsigned char *bfp = (const unsigned char *)bbv->c_str();
  return CalcBitmapTversky(afp, bfp, abv->size(), ca, cb);
}

/*******************************************
 *     CSfp transformation                 *
 *******************************************/

extern "C" void freeCSfp(CSfp data) {
  SparseFP *fp = (SparseFP *)data;
  delete fp;
}

extern "C" CSfp constructCSfp(Sfp *data) {
  SparseFP *ebv = NULL;

  try {
    ebv = new SparseFP(VARDATA(data), VARSIZE(data) - VARHDRSZ);
  } catch (...) {
    elog(ERROR, "constructMolFingerPrint: Unknown exception");
  }

  return (CSfp)ebv;
}

extern "C" Sfp *deconstructCSfp(CSfp data) {
  SparseFP *ebv = (SparseFP *)data;
  ByteA b;

  try {
    b = ebv->toString();
  } catch (...) {
    elog(ERROR, "deconstructMolFingerPrint: Unknown exception");
  }

  return b.toByteA();
}

extern "C" bytea *makeSfpSignature(CSfp data, int numBits) {
  SparseFP *v = (SparseFP *)data;
  int n, numBytes;
  bytea *res;
  unsigned char *s;
  SparseFP::StorageType::const_iterator iter;

  numBytes = VARHDRSZ + (numBits / 8);
  if ((numBits % 8) != 0) numBytes++;

  res = (bytea *)palloc0(numBytes);
  SET_VARSIZE(res, numBytes);
  s = (unsigned char *)VARDATA(res);

  for (iter = v->getNonzeroElements().begin();
       iter != v->getNonzeroElements().end(); iter++) {
    n = iter->first % numBits;
    s[n / 8] |= 1 << (n % 8);
  }

  return res;
}

extern "C" bytea *makeLowSparseFingerPrint(CSfp data, int numInts) {
  SparseFP *v = (SparseFP *)data;
  int numBytes;
  bytea *res;
  IntRange *s;
  int n;
  SparseFP::StorageType::const_iterator iter;

  numBytes = VARHDRSZ + (numInts * sizeof(IntRange));

  res = (bytea *)palloc0(numBytes);
  SET_VARSIZE(res, numBytes);
  s = (IntRange *)VARDATA(res);

  for (iter = v->getNonzeroElements().begin();
       iter != v->getNonzeroElements().end(); iter++) {
    uint32 iterV = (uint32)iter->second;
    n = iter->first % numInts;

    if (iterV > INTRANGEMAX) {
#if 0
        elog(ERROR, "sparse fingerprint is too big, increase INTRANGEMAX in rdkit.h");
#else
      iterV = INTRANGEMAX;
#endif
    }

    if (s[n].low == 0 || s[n].low > iterV) s[n].low = iterV;
    if (s[n].high < iterV) s[n].high = iterV;
  }

  return res;
}

extern "C" void countOverlapValues(bytea *sign, CSfp data,
                                   int numBits, int *sum, int *overlapSum,
                                   int *overlapN) {
  SparseFP *v = (SparseFP *)data;
  SparseFP::StorageType::const_iterator iter;

  *sum = *overlapSum = *overlapN = 0;

  if (sign) {
    unsigned char *s = (unsigned char *)VARDATA(sign);
    int n;

    for (iter = v->getNonzeroElements().begin();
         iter != v->getNonzeroElements().end(); iter++) {
      *sum += iter->second;
      n = iter->first % numBits;
      if (s[n / 8] & (1 << (n % 8))) {
        *overlapSum += iter->second;
        *overlapN += 1;
      }
    }
  } else {
    /* Assume, sign has only true bits */
    for (iter = v->getNonzeroElements().begin();
         iter != v->getNonzeroElements().end(); iter++)
      *sum += iter->second;

    *overlapSum = *sum;
    *overlapN = v->getNonzeroElements().size();
  }
}

extern "C" void countLowOverlapValues(bytea *sign, CSfp data,
                                      int numInts, int *querySum, int *keySum,
                                      int *overlapUp, int *overlapDown) {
  SparseFP *v = (SparseFP *)data;
  SparseFP::StorageType::const_iterator iter;
  IntRange *s = (IntRange *)VARDATA(sign);
  int n;

  *querySum = *keySum = *overlapUp = *overlapDown = 0;

  for (iter = v->getNonzeroElements().begin();
       iter != v->getNonzeroElements().end(); iter++) {
    *querySum += iter->second;
    n = iter->first % numInts;
    if (s[n].low == 0) {
      Assert(s[n].high == 0);
      continue;
    }

    *overlapDown += Min(s[n].low, (uint32)iter->second);
    *overlapUp += Min(s[n].high, (uint32)iter->second);
  }

  Assert(*overlapDown <= *overlapUp);

  for (n = 0; n < numInts; n++) {
    *keySum += s[n].low;
    if (s[n].low != s[n].high)
      *keySum +=
          s[n].high; /* there is at least two key mapped into current backet */
  }

  Assert(*overlapUp <= *keySum);
}

extern "C" double calcSparseTanimotoSml(CSfp a, CSfp b) {
  double res = -1.0;

  /*
   * Nsame / (Na + Nb - Nsame)
   */

  try {
    res = TanimotoSimilarity(*(SparseFP *)a, *(SparseFP *)b);
  } catch (ValueErrorException &e) {
    elog(ERROR, "TanimotoSimilarity: %s", e.message().c_str());
  } catch (...) {
    elog(ERROR, "calcSparseTanimotoSml: Unknown exception");
  }

  return res;
}

extern "C" double calcSparseDiceSml(CSfp a, CSfp b) {
  double res = -1.0;

  /*
   * 2 * Nsame / (Na + Nb)
   */

  try {
    res = DiceSimilarity(*(SparseFP *)a, *(SparseFP *)b);
  } catch (ValueErrorException &e) {
    elog(ERROR, "DiceSimilarity: %s", e.message().c_str());
  } catch (...) {
    elog(ERROR, "calcSparseDiceSml: Unknown exception");
  }

  return res;
}

extern "C" double calcSparseStringDiceSml(const char *a, unsigned int sza,
                                          const char *b, unsigned int szb) {
  const unsigned char *t1 = (const unsigned char *)a;
  const unsigned char *t2 = (const unsigned char *)b;

  boost::uint32_t tmp;
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != (boost::uint32_t)ci_SPARSEINTVECT_VERSION) {
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 1");
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2 += sizeof(boost::uint32_t);
  if (tmp != (boost::uint32_t)ci_SPARSEINTVECT_VERSION) {
    elog(ERROR, "calcSparseStringDiceSml: could not convert argument 2");
  }

  // check the element size:
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != sizeof(boost::uint32_t)) {
    elog(ERROR,
         "calcSparseStringDiceSml: could not convert argument 1 -> uint32_t");
  }
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2 += sizeof(boost::uint32_t);
  if (tmp != sizeof(boost::uint32_t)) {
    elog(ERROR,
         "calcSparseStringDiceSml: could not convert argument 2 -> uint32_t");
  }

  double res = 0.;
  // start reading:
  boost::uint32_t len1, len2;
  len1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  len2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2 += sizeof(boost::uint32_t);
  if (len1 != len2) {
    elog(ERROR, "attempt to compare fingerprints of different length");
  }

  boost::uint32_t nElem1, nElem2;
  nElem1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  nElem2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2 += sizeof(boost::uint32_t);

  if (!nElem1 || !nElem2) {
    return 0.0;
  }

  double v1Sum = 0, v2Sum = 0, numer = 0;
  boost::uint32_t idx1 = 0;
  boost::int32_t v1;
  boost::uint32_t idx2 = 0;
  boost::int32_t v2;
  idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
  t1 += sizeof(boost::int32_t);
  nElem1--;
  v1Sum += v1;

  idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
  t2 += sizeof(boost::uint32_t);
  v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
  t2 += sizeof(boost::int32_t);
  nElem2--;
  v2Sum += v2;

  while (1) {
    while (nElem2 && idx2 < idx1) {
      idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
      t2 += sizeof(boost::uint32_t);
      v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
      t2 += sizeof(boost::int32_t);
      nElem2--;
      v2Sum += v2;
    }
    if (idx2 == idx1) {
      // std::cerr<<"   --- "<<idx1<<" "<<v1<<" - "<<idx2<<" "<<v2<<std::endl;
      numer += std::min(v1, v2);
    }
    if (nElem1) {
      idx1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
      t1 += sizeof(boost::uint32_t);
      v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
      t1 += sizeof(boost::int32_t);
      nElem1--;
      v1Sum += v1;
    } else {
      break;
    }
  }
  while (nElem2) {
    idx2 = *(reinterpret_cast<const boost::uint32_t *>(t2));
    t2 += sizeof(boost::uint32_t);
    v2 = *(reinterpret_cast<const boost::int32_t *>(t2));
    t2 += sizeof(boost::int32_t);
    nElem2--;
    v2Sum += v2;
  }
  double denom = v1Sum + v2Sum;
  if (fabs(denom) < 1e-6) {
    res = 0.0;
  } else {
    res = 2. * numer / denom;
  }

  return res;
}

extern "C" bool calcSparseStringAllValsGT(const char *a, unsigned int sza,
                                          int tgt) {
  const unsigned char *t1 = (const unsigned char *)a;

  boost::uint32_t tmp;
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != (boost::uint32_t)ci_SPARSEINTVECT_VERSION) {
    elog(ERROR, "calcSparseStringAllValsGT: could not convert argument 1");
  }
  // check the element size:
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != sizeof(boost::uint32_t)) {
    elog(ERROR,
         "calcSparseStringAllValsGT: could not convert argument 1 -> uint32_t");
  }

  // boost::uint32_t len1;
  // len1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);

  boost::uint32_t nElem1;
  nElem1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);

  while (nElem1) {
    --nElem1;
    // skip the index:
    t1 += sizeof(boost::uint32_t);
    boost::int32_t v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
    t1 += sizeof(boost::int32_t);

    if (v1 <= tgt) return false;
  }
  return true;
}
extern "C" bool calcSparseStringAllValsLT(const char *a, unsigned int sza,
                                          int tgt) {
  const unsigned char *t1 = (const unsigned char *)a;

  boost::uint32_t tmp;
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != (boost::uint32_t)ci_SPARSEINTVECT_VERSION) {
    elog(ERROR, "calcSparseStringAllValsGT: could not convert argument 1");
  }
  // check the element size:
  tmp = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);
  if (tmp != sizeof(boost::uint32_t)) {
    elog(ERROR,
         "calcSparseStringAllValsGT: could not convert argument 1 -> uint32_t");
  }

  // boost::uint32_t len1;
  // len1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);

  boost::uint32_t nElem1;
  nElem1 = *(reinterpret_cast<const boost::uint32_t *>(t1));
  t1 += sizeof(boost::uint32_t);

  while (nElem1) {
    --nElem1;
    // skip the index:
    t1 += sizeof(boost::uint32_t);
    boost::int32_t v1 = *(reinterpret_cast<const boost::int32_t *>(t1));
    t1 += sizeof(boost::int32_t);

    if (v1 >= tgt) return false;
  }
  return true;
}

extern "C" CSfp addSFP(CSfp a, CSfp b) {
  SparseFP *res = NULL;
  try {
    SparseFP tmp = (*(SparseFP *)a + *(SparseFP *)b);
    res = (SparseFP *)new SparseFP(tmp);
  } catch (...) {
    elog(ERROR, "addSFP: Unknown exception");
  }
  return (CSfp)res;
}

extern "C" CSfp subtractSFP(CSfp a, CSfp b) {
  SparseFP *res = NULL;
  try {
    SparseFP tmp = (*(SparseFP *)a - *(SparseFP *)b);
    res = (SparseFP *)new SparseFP(tmp);
  } catch (...) {
    elog(ERROR, "addSFP: Unknown exception");
  }
  return (CSfp)res;
}

/*
 * Mol -> fp
 */
extern "C" CBfp makeLayeredBFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;

  try {
    res = RDKit::LayeredFingerprintMol(*mol, 0xFFFFFFFF, 1, 7,
                                       getLayeredFpSize());
  } catch (...) {
    elog(ERROR, "makeLayeredBFP: Unknown exception");
    if (res) delete res;
    res = NULL;
  }
  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CBfp makeRDKitBFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;

  try {
    res = RDKit::RDKFingerprintMol(*mol, 1, 6, getRDKitFpSize(), 2);
  } catch (...) {
    elog(ERROR, "makeRDKitBFP: Unknown exception");
    if (res) delete res;
    res = NULL;
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CSfp makeMorganSFP(CROMol data, int radius) {
  ROMol *mol = (ROMol *)data;
  SparseFP *res = NULL;
  std::vector<boost::uint32_t> invars(mol->getNumAtoms());
  try {
    RDKit::MorganFingerprints::getConnectivityInvariants(*mol, invars, true);
    res = (SparseFP *)RDKit::MorganFingerprints::getFingerprint(*mol, radius,
                                                                &invars);
  } catch (...) {
    elog(ERROR, "makeMorganSFP: Unknown exception");
  }

  return (CSfp)res;
}

extern "C" CBfp makeMorganBFP(CROMol data, int radius) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  std::vector<boost::uint32_t> invars(mol->getNumAtoms());
  try {
    RDKit::MorganFingerprints::getConnectivityInvariants(*mol, invars, true);
    res = RDKit::MorganFingerprints::getFingerprintAsBitVect(
        *mol, radius, getMorganFpSize(), &invars);
  } catch (...) {
    elog(ERROR, "makeMorganBFP: Unknown exception");
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CSfp makeFeatMorganSFP(CROMol data, int radius) {
  ROMol *mol = (ROMol *)data;
  SparseFP *res = NULL;
  std::vector<boost::uint32_t> invars(mol->getNumAtoms());
  try {
    RDKit::MorganFingerprints::getFeatureInvariants(*mol, invars);
    res = (SparseFP *)RDKit::MorganFingerprints::getFingerprint(*mol, radius,
                                                                &invars);
  } catch (...) {
    elog(ERROR, "makeMorganSFP: Unknown exception");
  }

  return (CSfp)res;
}

extern "C" CBfp makeFeatMorganBFP(CROMol data, int radius) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  std::vector<boost::uint32_t> invars(mol->getNumAtoms());
  try {
    RDKit::MorganFingerprints::getFeatureInvariants(*mol, invars);
    res = RDKit::MorganFingerprints::getFingerprintAsBitVect(
        *mol, radius, getFeatMorganFpSize(), &invars);
  } catch (...) {
    elog(ERROR, "makeMorganBFP: Unknown exception");
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CSfp makeAtomPairSFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  SparseFP *res = NULL;
#ifdef UNHASHED_PAIR_FPS
  try {
    SparseIntVect<boost::int32_t> *afp =
        RDKit::AtomPairs::getAtomPairFingerprint(*mol);
    res = new SparseFP(1 << RDKit::AtomPairs::numAtomPairFingerprintBits);
    for (SparseIntVect<boost::int32_t>::StorageType::const_iterator iter =
             afp->getNonzeroElements().begin();
         iter != afp->getNonzeroElements().end(); ++iter) {
      res->setVal(iter->first, iter->second);
    }
    delete afp;
  } catch (...) {
    elog(ERROR, "makeAtomPairSFP: Unknown exception");
  }
#else
  try {
    SparseIntVect<boost::int32_t> *afp =
        RDKit::AtomPairs::getHashedAtomPairFingerprint(
            *mol, getHashedAtomPairFpSize());
    res = new SparseFP(getHashedAtomPairFpSize());
    for (SparseIntVect<boost::int32_t>::StorageType::const_iterator iter =
             afp->getNonzeroElements().begin();
         iter != afp->getNonzeroElements().end(); ++iter) {
      res->setVal(iter->first, iter->second);
    }
    delete afp;
  } catch (...) {
    elog(ERROR, "makeAtomPairSFP: Unknown exception");
  }
#endif
  return (CSfp)res;
}

extern "C" CSfp makeTopologicalTorsionSFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  SparseFP *res = NULL;

#ifdef UNHASHED_PAIR_FPS
  try {
    SparseIntVect<boost::int64_t> *afp =
        RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(
            *mol, boost::integer_traits<boost::uint32_t>::const_max);
    res = new SparseFP(boost::integer_traits<boost::uint32_t>::const_max);
    for (SparseIntVect<boost::int64_t>::StorageType::const_iterator iter =
             afp->getNonzeroElements().begin();
         iter != afp->getNonzeroElements().end(); ++iter) {
      res->setVal(iter->first, iter->second);
    }
    delete afp;
  } catch (...) {
    elog(ERROR, "makeTopologicalTorsionSFP: Unknown exception");
  }
#else
  try {
    SparseIntVect<boost::int64_t> *afp =
        RDKit::AtomPairs::getHashedTopologicalTorsionFingerprint(
            *mol, getHashedTorsionFpSize());
    res = new SparseFP(getHashedTorsionFpSize());
    for (SparseIntVect<boost::int64_t>::StorageType::const_iterator iter =
             afp->getNonzeroElements().begin();
         iter != afp->getNonzeroElements().end(); ++iter) {
      res->setVal(iter->first, iter->second);
    }
    delete afp;
  } catch (...) {
    elog(ERROR, "makeTopologicalTorsionSFP: Unknown exception");
  }
#endif
  return (CSfp)res;
}

extern "C" CBfp makeAtomPairBFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  try {
    res = RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(
        *mol, getHashedAtomPairFpSize());
  } catch (...) {
    elog(ERROR, "makeAtomPairBFP: Unknown exception");
  }
  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CBfp makeTopologicalTorsionBFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  try {
    res = RDKit::AtomPairs::getHashedTopologicalTorsionFingerprintAsBitVect(
        *mol, getHashedTorsionFpSize());
  } catch (...) {
    elog(ERROR, "makeTopologicalTorsionBFP: Unknown exception");
  }
  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CBfp makeMACCSBFP(CROMol data) {
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  try {
    res = RDKit::MACCSFingerprints::getFingerprintAsBitVect(*mol);
  } catch (...) {
    elog(ERROR, "makeMACCSBFP: Unknown exception");
  }
  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" CBfp makeAvalonBFP(CROMol data, bool isQuery,
			      unsigned int bitFlags) {
#ifdef BUILD_AVALON_SUPPORT
  ROMol *mol = (ROMol *)data;
  ExplicitBitVect *res = NULL;
  try {
    res = new ExplicitBitVect(getAvalonFpSize());
    AvalonTools::getAvalonFP(*mol, *res, getAvalonFpSize(), isQuery, true,
                             bitFlags);
  } catch (...) {
    elog(ERROR, "makeAvalonBFP: Unknown exception");
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
#else
  elog(ERROR, "Avalon support not enabled");
  return NULL;
#endif
}

/* chemical reactions */

extern "C" void freeChemReaction(CChemicalReaction data) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  delete rxn;
}

extern "C" CChemicalReaction constructChemReact(Reaction *data) {
  ChemicalReaction *rxn = new ChemicalReaction();

  try {
    ByteA b(data);
    ReactionPickler::reactionFromPickle(b, rxn);
  } catch (ReactionPicklerException &e) {
    elog(ERROR, "reactionFromPickle: %s", e.message());
  } catch (...) {
    elog(ERROR, "constructChemReact: Unknown exception");
  }

  return (CChemicalReaction)rxn;
}

extern "C" Reaction *deconstructChemReact(CChemicalReaction data) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  ByteA b;

  try {
    ReactionPickler::pickleReaction(rxn, b);
  } catch (ReactionPicklerException &e) {
    elog(ERROR, "pickleReaction: %s", e.message());
  } catch (...) {
    elog(ERROR, "deconstructChemReact: Unknown exception");
  }

  return (Reaction *)b.toByteA();
}

extern "C" CChemicalReaction parseChemReactText(char *data, bool asSmarts,
                                                bool warnOnFail) {
  ChemicalReaction *rxn = NULL;

  try {
    if (asSmarts) {
      rxn = RxnSmartsToChemicalReaction(data);
    } else {
      rxn = RxnSmartsToChemicalReaction(data, 0, true);
    }
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    rxn = NULL;
  }
  if (rxn == NULL) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create chemical reaction from SMILES '%s'",
                      data)));
    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create chemical reaction  from SMILES '%s'",
                      data)));
    }
  }

  return (CChemicalReaction)rxn;
}

extern "C" CChemicalReaction parseChemReactBlob(char *data, int len) {
  ChemicalReaction *rxn = NULL;

  try {
    string binStr(data, len);
    rxn = new ChemicalReaction(binStr);
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("problem generating chemical reaction from blob data")));
  }
  if (rxn == NULL) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("blob data could not be parsed")));
  }

  return (CChemicalReaction)rxn;
}

extern "C" char *makeChemReactText(CChemicalReaction data, int *len,
				   bool asSmarts) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;

  try {
    if (!asSmarts) {
      StringData = ChemicalReactionToRxnSmiles(*rxn);
    } else {
      StringData = ChemicalReactionToRxnSmarts(*rxn);
    }
  } catch (...) {
    ereport(WARNING, (errcode(ERRCODE_WARNING),
                      errmsg("makeChemReactText: problems converting chemical "
                             "reaction  to SMILES/SMARTS")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" char *makeChemReactBlob(CChemicalReaction data, int *len) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  StringData.clear();
  try {
    ReactionPickler::pickleReaction(*rxn, StringData);
  } catch (...) {
    elog(ERROR, "makeChemReactBlob: Unknown exception");
  }

  *len = StringData.size();
  return (char *)StringData.data();
}

extern "C" CChemicalReaction parseChemReactCTAB(char *data, bool warnOnFail) {
  ChemicalReaction *rxn = NULL;

  try {
    rxn = RxnBlockToChemicalReaction(data);
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    rxn = NULL;
  }
  if (rxn == NULL) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create reaction from CTAB '%s'", data)));

    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create reaction from CTAB '%s'", data)));
    }
  }

  return (CChemicalReaction)rxn;
}

extern "C" char *makeCTABChemReact(CChemicalReaction data, int *len) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;

  try {
    StringData = ChemicalReactionToRxnBlock(*rxn);
  } catch (...) {
    ereport(
        WARNING,
        (errcode(ERRCODE_WARNING),
         errmsg("makeCTABChemReact: problems converting reaction to CTAB")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" int ChemReactNumReactants(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumReactantTemplates();
}

extern "C" int ChemReactNumProducts(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumProductTemplates();
}

extern "C" int ChemReactNumAgents(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumAgentTemplates();
}

extern "C" bytea *makeReactionSign(CChemicalReaction data) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  ExplicitBitVect *res = NULL;
  bytea *ret = NULL;

  try {
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(getReactionSubstructFpType());
    params.fpSize = getReactionSubstructFpSize();
    params.includeAgents = (!getIgnoreReactionAgents());
    params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();
    res = RDKit::StructuralFingerprintChemReaction(*rxn, params);

    if (res) {
      std::string sres = BitVectToBinaryText(*res);
      
      unsigned int varsize = VARHDRSZ + sres.size();
      ret = (bytea *) palloc0(varsize);
      memcpy(VARDATA(ret), sres.data(), sres.size());
      SET_VARSIZE(ret, varsize);
      
      delete res;
      res = 0;
    }
  } catch (...) {
    elog(ERROR, "makeReactionSign: Unknown exception");
    if (res) delete res;
  }
  return ret;
}

extern "C" int ReactionSubstruct(CChemicalReaction rxn,
                                 CChemicalReaction rxn2) {
  ChemicalReaction *rxnm = (ChemicalReaction *)rxn;
  ChemicalReaction *rxn2m = (ChemicalReaction *)rxn2;

  /* Reaction search */
  if (rxn2m->getNumReactantTemplates() != 0 &&
      rxn2m->getNumProductTemplates() != 0) {
    return hasReactionSubstructMatch(*rxnm, *rxn2m,
                                     (!getIgnoreReactionAgents()));
  }
  /* Product search */
  if (rxn2m->getNumReactantTemplates() == 0 &&
      rxn2m->getNumProductTemplates() != 0) {
    if (rxn2m->getNumAgentTemplates() != 0 && !getIgnoreReactionAgents()) {
      return (hasProductTemplateSubstructMatch(*rxnm, *rxn2m) &&
              hasAgentTemplateSubstructMatch(*rxnm, *rxn2m));
    }
    return hasProductTemplateSubstructMatch(*rxnm, *rxn2m);
  }
  /* Reactant search */
  if (rxn2m->getNumReactantTemplates() != 0 &&
      rxn2m->getNumProductTemplates() == 0) {
    if (rxn2m->getNumAgentTemplates() != 0 && !getIgnoreReactionAgents()) {
      return (hasReactantTemplateSubstructMatch(*rxnm, *rxn2m) &&
              hasAgentTemplateSubstructMatch(*rxnm, *rxn2m));
    }
    return hasReactantTemplateSubstructMatch(*rxnm, *rxn2m);
  }
  /* Agent search */
  if (rxn2m->getNumReactantTemplates() == 0 &&
      rxn2m->getNumProductTemplates() == 0 &&
      rxn2m->getNumAgentTemplates() != 0) {
    return hasAgentTemplateSubstructMatch(*rxnm, *rxn2m);
  }

  return false;
}

extern "C" int ReactionSubstructFP(CChemicalReaction rxn,
                                   CChemicalReaction rxnquery) {
  ChemicalReaction *rxnm = (ChemicalReaction *)rxn;
  ChemicalReaction *rxnqm = (ChemicalReaction *)rxnquery;

  RDKit::ReactionFingerprintParams params;
  params.fpType = static_cast<FingerprintType>(getReactionSubstructFpType());
  params.fpSize = getReactionSubstructFpSize();
  params.includeAgents = (!getIgnoreReactionAgents());
  params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();

  ExplicitBitVect *fp1 = StructuralFingerprintChemReaction(*rxnm, params);
  ExplicitBitVect *fp2 = StructuralFingerprintChemReaction(*rxnqm, params);

  if (fp1->getNumOnBits() < fp2->getNumOnBits()) {
    return false;
  }
  for (unsigned i = 0; i < fp1->getNumBits(); i++) {
    if ((fp1->getBit(i) & fp2->getBit(i)) != fp2->getBit(i)) {
      return false;
    }
  }
  return true;
}

// some helper functions in anonymous namespace
namespace {

struct MoleculeDescriptors {
  MoleculeDescriptors() : nAtoms(0), nBonds(0), nRings(0), MW(0.0) {}
  unsigned nAtoms;
  unsigned nBonds;
  unsigned nRings;
  double MW;
};

MoleculeDescriptors *calcMolecularDescriptorsReaction(
    RDKit::ChemicalReaction *rxn, RDKit::ReactionMoleculeType t) {
  MoleculeDescriptors *des = new MoleculeDescriptors();
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(*rxn, t);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(*rxn, t);
  for (; begin != end; ++begin) {
    des->nAtoms += begin->get()->getNumHeavyAtoms();
    des->nBonds += begin->get()->getNumBonds(true);
    des->MW = RDKit::Descriptors::calcAMW(*begin->get(), true);
    if (!begin->get()->getRingInfo()->isInitialized()) {
      begin->get()->updatePropertyCache();
      RDKit::MolOps::findSSSR(*begin->get());
    }
    des->nRings += begin->get()->getRingInfo()->numRings();
  }
  return des;
}

int compareMolDescriptors(const MoleculeDescriptors &md1,
                          const MoleculeDescriptors &md2) {
  int res = md1.nAtoms - md2.nAtoms;
  if (res) {
    return res;
  }
  res = md1.nBonds - md2.nBonds;
  if (res) {
    return res;
  }
  res = md1.nRings - md2.nRings;
  if (res) {
    return res;
  }
  res = int(md1.MW - md2.MW);
  if (res) {
    return res;
  }
  return 0;
}
}

extern "C" int reactioncmp(CChemicalReaction rxn, CChemicalReaction rxn2) {
  ChemicalReaction *rxnm = (ChemicalReaction *)rxn;
  ChemicalReaction *rxn2m = (ChemicalReaction *)rxn2;

  if (!rxnm) {
    if (!rxn2m) return 0;
    return -1;
  }
  if (!rxn2m) return 1;

  int res = rxnm->getNumReactantTemplates() - rxn2m->getNumReactantTemplates();
  if (res) {
    return res;
  }
  res = rxnm->getNumProductTemplates() - rxn2m->getNumProductTemplates();
  if (res) {
    return res;
  }
  if (!getIgnoreReactionAgents()) {
    res = rxnm->getNumAgentTemplates() - rxn2m->getNumAgentTemplates();
    if (res) {
      return res;
    }
  }

  MoleculeDescriptors *rxn_react =
      calcMolecularDescriptorsReaction(rxnm, Reactant);
  MoleculeDescriptors *rxn2_react =
      calcMolecularDescriptorsReaction(rxn2m, Reactant);
  res = compareMolDescriptors(*rxn_react, *rxn2_react);
  delete (rxn_react);
  delete (rxn2_react);
  if (res) {
    return res;
  }
  MoleculeDescriptors *rxn_product =
      calcMolecularDescriptorsReaction(rxnm, Product);
  MoleculeDescriptors *rxn2_product =
      calcMolecularDescriptorsReaction(rxn2m, Product);
  res = compareMolDescriptors(*rxn_product, *rxn2_product);
  delete (rxn_product);
  delete (rxn2_product);
  if (res) {
    return res;
  }
  if (!getIgnoreReactionAgents()) {
    MoleculeDescriptors *rxn_agent =
        calcMolecularDescriptorsReaction(rxnm, Agent);
    MoleculeDescriptors *rxn2_agent =
        calcMolecularDescriptorsReaction(rxn2m, Agent);
    res = compareMolDescriptors(*rxn_agent, *rxn2_agent);
    delete (rxn_agent);
    delete (rxn2_agent);
    if (res) {
      return res;
    }
  }

  RDKit::MatchVectType matchVect;
  if (hasReactionSubstructMatch(*rxnm, *rxn2m, (!getIgnoreReactionAgents()))) {
    return 0;
  }
  return -1;
}

extern "C" CSfp makeReactionDifferenceSFP(
    CChemicalReaction data, int size, int fpType) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  SparseFP *res = NULL;

  try {
    if (fpType > 3 || fpType < 1) {
      elog(ERROR, "makeReactionDifferenceSFP: Unknown Fingerprint type");
    }
    FingerprintType fp = static_cast<RDKit::FingerprintType>(fpType);
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(fpType);
    params.fpSize = size;
    params.includeAgents = (!getIgnoreReactionAgents());
    params.agentWeight = getReactionDifferenceFPWeightAgents();
    params.nonAgentWeight = getReactionDifferenceFPWeightNonagents();
    res = (SparseFP *)RDKit::DifferenceFingerprintChemReaction(*rxn, params);
  } catch (...) {
    elog(ERROR, "makeReactionDifferenceSFP: Unknown exception");
  }
  return (CSfp)res;
}

extern "C" CBfp makeReactionBFP(CChemicalReaction data, int size, int fpType) {
  ChemicalReaction *rxn = (ChemicalReaction *)data;
  ExplicitBitVect *res = NULL;

  try {
    if (fpType > 5 || fpType < 1) {
      elog(ERROR, "makeReactionBFP: Unknown Fingerprint type");
    }
    FingerprintType fp = static_cast<RDKit::FingerprintType>(fpType);
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(fpType);
    params.fpSize = size;
    params.includeAgents = (!getIgnoreReactionAgents());
    params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();
    res = (ExplicitBitVect *)RDKit::StructuralFingerprintChemReaction(*rxn,
                                                                      params);
  } catch (...) {
    elog(ERROR, "makeReactionBFP: Unknown exception");
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return NULL;
  }
}

extern "C" char *computeMolHash(CROMol data, int *len) {
  ROMol &mol = *(ROMol *)data;
  static string text;
  text.clear();
  try {
    // FIX: once R/S values are stored on the atoms, this will no longer be
    // needed
    MolOps::assignStereochemistry(mol);
    text = RDKit::MolHash::generateMoleculeHashSet(mol);
  } catch (...) {
    ereport(WARNING,
            (errcode(ERRCODE_WARNING), errmsg("computeMolHash: failed")));
    text.clear();
  }
  *len = text.length();
  return (char *)text.c_str();
}

extern "C" char *  // TEMP
    Mol2Smiles(CROMol data) {
  const ROMol &mol = *(ROMol *)data;
  static string text;
  text.clear();
  try {
    text = RDKit::MolToSmiles(mol);
  } catch (...) {
    ereport(WARNING,
            (errcode(ERRCODE_WARNING), errmsg("Mol2Smiles(): failed")));
    text.clear();
  }
  return (char *)text.c_str();
}

extern "C" char *findMCSsmiles(char *smiles, char *params) {
  static string mcs;
  mcs.clear();

  char *str = smiles;
  char *s = str;
  int len, nmols = 0;
  std::vector<RDKit::ROMOL_SPTR> molecules;
  while (*s && *s <= ' ') s++;
  while (*s > ' ') {
    len = 0;
    while (s[len] > ' ') len++;
    s[len] = '\0';
    if (0 == strlen(s)) continue;
    molecules.push_back(RDKit::ROMOL_SPTR(RDKit::SmilesToMol(s)));
    // elog(WARNING, s);
    s += len;
    s++;  // do s++; while(*s && *s <= ' ');
  }

  RDKit::MCSParameters p;

  if (params && 0 != strlen(params)) {
    try {
      RDKit::parseMCSParametersJSON(params, &p);
    } catch (...) {
      ereport(WARNING, (errcode(ERRCODE_WARNING),
                        errmsg("findMCS: Invalid argument \'params\'")));
      return (char *)mcs.c_str();
    }
  }

  try {
    MCSResult res = RDKit::findMCS(molecules, &p);
    mcs = res.SmartsString;
    if (!res.isCompleted())
      ereport(WARNING, (errcode(ERRCODE_WARNING),
                        errmsg("findMCS timed out, result is not maximal")));
  } catch (...) {
    ereport(WARNING, (errcode(ERRCODE_WARNING), errmsg("findMCS: failed")));
    mcs.clear();
  }
  return mcs.empty() ? (char *)"" : (char *)mcs.c_str();
}

extern "C" void *addMol2list(void *lst, Mol *mol) {
  try {
    if (!lst) {
      // elog(WARNING, "addMol2list: allocate new list");
      lst = new std::vector<RDKit::ROMOL_SPTR>;
    }
    std::vector<RDKit::ROMOL_SPTR> &mlst =
        *(std::vector<RDKit::ROMOL_SPTR> *)lst;
    // elog(WARNING, "addMol2list: create a copy of new mol");
    ROMol *m = (ROMol *)constructROMol(
        mol);  // new ROMol(*(const ROMol*)mol, false); // create a copy
    // elog(WARNING, "addMol2list: append new mol into list");
    mlst.push_back(RDKit::ROMOL_SPTR(m));
    // elog(WARNING, "addMol2list: finished");
  } catch (...) {
    // elog(WARNING, "addMol2list: ERROR");
    ereport(WARNING, (errcode(ERRCODE_WARNING), errmsg("addMol2list: failed")));
  }
  return lst;
}

extern "C" char *findMCS(void *vmols, char *params) {
  static string mcs;
  mcs.clear();
  std::vector<RDKit::ROMOL_SPTR> *molecules =
      (std::vector<RDKit::ROMOL_SPTR> *)vmols;
  // char t[256];
  // sprintf(t,"findMCS(): lst=%p, size=%u", molecules, molecules->size());
  // elog(WARNING, t);

  RDKit::MCSParameters p;

  if (params && 0 != strlen(params)) {
    try {
      RDKit::parseMCSParametersJSON(params, &p);
    } catch (...) {
      // mcs = params; //DEBUG
      ereport(WARNING, (errcode(ERRCODE_WARNING),
                        errmsg("findMCS: Invalid argument \'params\'")));
      return (char *)mcs.c_str();
    }
  }

  try {
    MCSResult res = RDKit::findMCS(*molecules, &p);
    if (!res.isCompleted())
      ereport(WARNING, (errcode(ERRCODE_WARNING),
                        errmsg("findMCS timed out, result is not maximal")));
    mcs = res.SmartsString;
  } catch (...) {
    ereport(WARNING, (errcode(ERRCODE_WARNING), errmsg("findMCS: failed")));
    mcs.clear();
  }
  // sprintf(t,"findMCS(): MCS='%s'", mcs.c_str());
  // elog(WARNING, t);
  delete molecules;
  // elog(WARNING, "findMCS(): molecules deleted. FINISHED.");
  return (char *)mcs.c_str();
}
