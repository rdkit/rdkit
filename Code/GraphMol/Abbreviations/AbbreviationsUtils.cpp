//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Abbreviations.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitQueries.h>
#include <boost/tokenizer.hpp>

using tokenizer = boost::tokenizer<boost::char_separator<char>>;

namespace RDKit {

namespace Abbreviations {

namespace common_properties {
const std::string numDummies = "_numDummies";
const std::string origAtomMapping = "_origAtomMapping";
const std::string origBondMapping = "_origBondMapping";
}  // namespace common_properties

namespace Utils {
namespace data {
/*
Translations of superatom labels to SMILES.

First atom of SMILES string should be the one connected to the rest of
the molecule.

ADAPTED FROM: https://github.com/openbabel/superatoms/blob/master/superatom.txt

Originally from http://cactus.nci.nih.gov/osra/

The left-aligned form is the one recognized in MDL alias lines;
the right-aligned form may be used in 2D depiction.

label smiles display_label display_label_w
*/
const std::string defaultAbbreviations =
    R"ABBREVS(CO2Et    C(=O)OCC  CO<sub>2</sub>Et    EtO<sub>2</sub>C
COOEt    C(=O)OCC  CO<sub>2</sub>Et    EtO<sub>2</sub>C
OiBu     OCC(C)C   OiBu     iBuO
nDec     CCCCCCCCCC  nDec
nNon     CCCCCCCCC   nNon
nOct     CCCCCCCC    nOct
nHept    CCCCCCC     nHept
nHex     CCCCCC      nHex
nPent    CCCCC       nPent
iPent    C(C)CCC     iPent
tBu      C(C)(C)C    tBu
iBu      C(C)CC      iBu
nBu      CCCC        nBu
iPr      C(C)C       iPr
nPr      CCC         nPr
Et       CC          Et
NCF3     NC(F)(F)F NCF<sub>3</sub>     F<sub>3</sub>CN
CF3      C(F)(F)F  CF<sub>3</sub>      F<sub>3</sub>C
CCl3     C(Cl)(Cl)Cl CCl<sub>3</sub>     Cl<sub>3</sub>C
CN       C#N       CN       NC
NC       [N+]#[C-] NC       CN 
N(OH)CH3 N([OH])C    N(OH)CH<sub>3</sub> CH<sub>3</sub>(OH)N
NO2      [N+](=O)[O-]  NO<sub>2</sub>      O<sub>2</sub>N
NO       N=O   NO       ON
SO3H     S(=O)(=O)[OH] SO<sub>3</sub>H     HO<sub>3</sub>S
CO2H     C(=O)[OH] CO<sub>2</sub>H     HO<sub>2</sub>C
COOH     C(=O)[OH] COOH     HOOC
OEt      OCC   OEt      EtO
OAc      OC(=O)C   OAc      AcO
NHAc     NC(=O)C   NHAc     AcNH
Ac       C(=O)C    Ac
CHO      C=O   CHO      OHC
NMe      NC    NMe      MeN
SMe      SC    SMe      MeS
OMe      OC    OMe      MeO
CO2-     C(=O)[O-]   COO<sup>-</sup>     <sup>-</sup>OOC
COO-     C(=O)[O-]   COO<sup>-</sup>     <sup>-</sup>OOC)ABBREVS";

/*
Translations of linker superatom labels to SMILES.

First atom of SMILES string should be a dummy connected to the rest of
the molecule. The other linker dummy/dummies show the other attachments

*/
const std::string defaultLinkers =
    R"ABBREVS(PEG6  *OCCOCCOCCOCCOCCOCC* PEG6
PEG5  *OCCOCCOCCOCCOCC* PEG5
PEG4  *OCCOCCOCCOCC* PEG4
PEG3  *OCCOCCOCC* PEG3
Dec   *CCCCCCCCCC*
Non   *CCCCCCCCC*
Oct   *CCCCCCCC*
Hept  *CCCCCCC*)ABBREVS";
// other possible abbreviations that might be useful:
/*
PEG6  *OCCOCCOCCOCCOCC* PEG6
PEG5  *OCCOCCOCCOCCOCC* PEG5
PEG4  *OCCOCCOCCOCC* PEG4
PEG3  *OCCOCCOCC* PEG3
Dec   *CCCCCCCCCC*
Non   *CCCCCCCCC*
Oct   *CCCCCCCC*
Hept  *CCCCCCC*
Hex   *CCCCCC*
Pent  *CCCCC*
Cy   *C1CCC(*)CC1  Cy
ala *N[C@@H](C)C(=O)* ala
arg *N[C@@H](CCCNC(N)=[NH])C(=O)* arg
asn *N[C@@H](CC(N)=O)C(=O)* asn
asp *N[C@@H](CC(O)=O)C(=O)* asp
cys *N[C@@H](CS)C(=O)* cys
gln *N[C@@H](CCC(N)=O)C(=O)* gln
glu *N[C@@H](CCC(O)=O)C(=O)* glu
gly *NCC(=O)* gly
his *N[C@@H](Cc1c[nH]cn1)C(=O)* his
ile *N[C@@H](C(C)CC)C(=O)* ile
leu *N[C@@H](CC(C)C)C(=O)* leu
lys *N[C@@H](CCCCN)C(=O)* lys
met *N[C@@H](CCSC)C(=O)* met
phe *N[C@@H](Cc1ccccc1)C(=O)* phe
pro *N1[C@@H](CCC1)C(=O)* pro
ser *N[C@@H](CO)C(=O)* ser
thr *N[C@@H](C(O)C)C(=O)* thr
trp *N[C@@H](Cc1c[nH]c2ccccc21)C(=O)* trp
tyr *N[C@@H](Cc1ccc(O)cc1)C(=O)* tyr
val *N[C@@H](C(C)C)C(=O)* val
*/
}  // namespace data

namespace detail {
ROMol *createAbbreviationMol(const std::string &txt, bool removeExtraDummies,
                             bool allowConnectionToDummies) {
  std::string smarts;
  if (txt[0] != '*') {
    smarts = "*" + txt;
  } else {
    smarts = txt;
  }
  RWMol *q = SmartsToMol(smarts);
  if (!q) {
    return q;
  }
  if (q->getNumAtoms() < 2) {
    BOOST_LOG(rdErrorLog) << "abbreviation with <2 atoms ignored" << std::endl;
    delete q;
    return nullptr;
  }
  MolOps::AdjustQueryParameters ps;
  ps.adjustDegree = true;
  ps.adjustDegreeFlags = MolOps::AdjustQueryWhichFlags::ADJUST_IGNOREDUMMIES;
  ps.adjustRingCount = true;
  ps.adjustRingCountFlags = MolOps::AdjustQueryWhichFlags::ADJUST_IGNOREDUMMIES;
  MolOps::adjustQueryProperties(*q, &ps);
  if (!allowConnectionToDummies) {
    auto qry = makeAtomNumQuery(0);
    qry->setNegation(true);
    q->getAtomWithIdx(0)->expandQuery(qry);
  }
  unsigned int nDummies = std::count_if(smarts.begin(), smarts.end(),
                                        [](char c) { return c == '*'; });
  if (removeExtraDummies) {
    q->beginBatchEdit();
    for (const auto at : q->atoms()) {
      if (!at->getIdx()) {
        // skip the first atom
        continue;
      }
      if (at->hasQuery() && at->getQuery()->getDescription() == "AtomNull") {
        q->removeAtom(at->getIdx());
        --nDummies;
      }
    }
    q->commitBatchEdit();
  }
  q->setProp(common_properties::numDummies, nDummies);
  return q;
}
}  // namespace detail

std::vector<AbbreviationDefinition> parseAbbreviations(
    const std::string &text, bool removeExtraDummies,
    bool allowConnectionToDummies) {
  std::vector<AbbreviationDefinition> res;
  boost::char_separator<char> lineSep("\n");
  tokenizer lines(text, lineSep);
  boost::char_separator<char> fieldSep(" \t");
  for (const auto &line : lines) {
    AbbreviationDefinition defn;
    tokenizer fields(line, fieldSep);
    tokenizer::iterator field = fields.begin();
    defn.label = *field;
    ++field;
    defn.smarts = *field;
    ++field;
    if (field != fields.end()) {
      defn.displayLabel = *field;
      ++field;
      if (field != fields.end()) {
        defn.displayLabelW = *field;
      }
    }
    defn.mol.reset(detail::createAbbreviationMol(
        defn.smarts, removeExtraDummies, allowConnectionToDummies));
    if (defn.mol) {
      res.push_back(defn);
    }
  }

  return res;
}
std::vector<AbbreviationDefinition> getDefaultAbbreviations() {
  static auto defs = parseAbbreviations(data::defaultAbbreviations);
  return defs;
}
std::vector<AbbreviationDefinition> getDefaultLinkers() {
  static auto defs = parseAbbreviations(data::defaultLinkers, true, true);
  return defs;
}
}  // namespace Utils

}  // namespace Abbreviations
}  // namespace RDKit
