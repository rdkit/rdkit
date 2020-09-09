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
}

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

left    right    SMILES
*/
const std::string defaultAbbreviations =
    R"ABBREVS(CO2Et    EtO2C    C(=O)OCC
COOEt    EtOOC    C(=O)OCC
OiBu     iBuO     OCC(C)C
tBu      tBu      C(C)(C)C
nBu      nBu      CCCC
iPr      iPr      C(C)C
nPr      nPr      CCC
Et       Et       CC
NCF3     F3CN     NC(F)(F)F
CF3      F3C      C(F)(F)F
CCl3     Cl3C     C(Cl)(Cl)Cl
CN       NC       C#N
NC       CN       [N+]#[C-]
N(OH)CH3 CH3(OH)N N([OH])C
NO2      O2N      [N+](=O)[O-]
NO       ON       N=O
SO3H     HO3S     S(=O)(=O)[OH]
CO2H     HOOC     C(=O)[OH]
COOH     HOOC     C(=O)[OH]
OEt      EtO      OCC
OAc      AcO      OC(=O)C
NHAc     AcNH     NC(=O)C
Ac       Ac       C(=O)C
CHO      OHC      C=O
NMe      MeN      NC
SMe      MeS      SC
OMe      MeO      OC
CO2-     -OOC     C(=O)[O-]
COO-     -OOC     C(=O)[O-])ABBREVS";

/*
Translations of linker superatom labels to SMILES.

First atom of SMILES string should be a dummy connected to the rest of
the molecule. The other linker dummy/dummies show the other attachments

*/
const std::string defaultLinkers = R"ABBREVS(PEG4  PEG4    *OCCOCCOCCOCC*
PEG3  PEG3    *OCCOCCOCC*
pentyl  pentyl  *CCCCC*
cyhex cyhex   *C1CCC(*)CC1
ala ala *N[C@@H](C)C(=O)*
arg arg *N[C@@H](CCCNC(N)=[NH])C(=O)*
asn asn *N[C@@H](CC(N)=O)C(=O)*
asp asp *N[C@@H](CC(O)=O)C(=O)*
cys cys *N[C@@H](CS)C(=O)*
gln gln *N[C@@H](CCC(N)=O)C(=O)*
glu glu *N[C@@H](CCC(O)=O)C(=O)*
gly gly *NCC(=O)*
his his *N[C@@H](Cc1c[nH]cn1)C(=O)*
ile ile *N[C@@H](C(C)CC)C(=O)*
leu leu *N[C@@H](CC(C)C)C(=O)*
lys lys *N[C@@H](CCCCN)C(=O)*
met met *N[C@@H](CCSC)C(=O)*
phe phe *N[C@@H](Cc1ccccc1)C(=O)*
pro pro *N1[C@@H](CCC1)C(=O)*
ser ser *N[C@@H](CO)C(=O)*
thr thr *N[C@@H](C(O)C)C(=O)*
trp trp *N[C@@H](Cc1c[nH]c2ccccc21)C(=O)*
tyr tyr *N[C@@H](Cc1ccc(O)cc1)C(=O)*
val val *N[C@@H](C(C)C)C(=O)*)ABBREVS";
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
    for (unsigned int i = q->getNumAtoms() - 1; i > 0; --i) {
      auto at = q->getAtomWithIdx(i);
      if (at->hasQuery() && at->getQuery()->getDescription() == "AtomNull") {
        q->removeAtom(i);
        --nDummies;
      }
    }
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
  for (const auto line : lines) {
    AbbreviationDefinition defn;
    tokenizer fields(line, fieldSep);
    tokenizer::iterator field = fields.begin();
    defn.llabel = *field;
    ++field;
    defn.rlabel = *field;
    ++field;
    defn.smarts = *field;
    defn.mol.reset(detail::createAbbreviationMol(
        defn.smarts, removeExtraDummies, allowConnectionToDummies));
    if (defn.mol) {
      res.push_back(defn);
    }
  }

  return res;
}
std::vector<AbbreviationDefinition> getDefaultAbbreviations() {
  return parseAbbreviations(data::defaultAbbreviations);
}
std::vector<AbbreviationDefinition> getDefaultLinkers() {
  return parseAbbreviations(data::defaultLinkers, true, true);
}
}  // namespace Utils

}  // namespace Abbreviations
}  // namespace RDKit
