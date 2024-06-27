//
//  Copyright (c) 2008-2023, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//   All rights reserved.
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
//       products derived from this software without specific prior
//       written permission.
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
//  created by Nik Stiefl May 2008
//  this file is heavily based on glandrum's MolFileParser
//

#include "FileParsers.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>
//
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>

namespace RDKit {

namespace {
void fixNitroSubstructureAndCharge(RWMol *res, unsigned int atIdx) {
  unsigned int noODblNeighbors = 0;
  ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
  unsigned int toModIdx = 0;
  boost::tie(nbrIdxIt, nbrEndIdxIt) =
      res->getAtomNeighbors(res->getAtomWithIdx(atIdx));
  while (nbrIdxIt != nbrEndIdxIt) {
    Bond *curBond = res->getBondBetweenAtoms(atIdx, *nbrIdxIt);
    if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
        curBond->getBondType() == Bond::DOUBLE) {
      ++noODblNeighbors;
      toModIdx = *nbrIdxIt;
    }
    ++nbrIdxIt;
  }
  if (noODblNeighbors == 2) {
    res->getBondBetweenAtoms(atIdx, toModIdx)->setBondType(Bond::SINGLE);
    res->getAtomWithIdx(atIdx)->setFormalCharge(1);
    res->getAtomWithIdx(toModIdx)->setFormalCharge(-1);
  }
}

void readFormalChargesFromAttr(std::istream *inStream, RWMol *res) {
  PRECONDITION(inStream, "inStream not valid");
  PRECONDITION(!inStream->eof(), "inStream is at eof");
  PRECONDITION(res, "RWMol not valid");
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(" \t\n");
  bool readNextAtomAttribs = true;
  unsigned int atomIdx = 0, noAtomAttr = 0;

  // std::streampos stPos = inStream->tellg();
  std::string tempStr = getLine(inStream);
  // there needs to be at least one entry
  if (inStream->eof()) {
    throw FileParseException("premature EOF in readFormalCharges");
  }
  while (readNextAtomAttribs) {
    tokenizer tokens(tempStr, sep);
    tokenizer::const_iterator itemIt = tokens.begin();
    try {
      atomIdx = boost::lexical_cast<unsigned int>(*itemIt);
      ++itemIt;
      noAtomAttr = boost::lexical_cast<unsigned int>(*itemIt);
    } catch (boost::bad_lexical_cast &) {
      throw FileParseException("Cannot process mol2 UnityAtomAttr.");
    }
    for (unsigned int i = 0; i < noAtomAttr; ++i) {
      std::string tempStr = getLine(inStream);
      if (inStream->eof()) {
        throw FileParseException("premature EOF in readFormalCharges");
      }
      tokenizer atmTokens(tempStr, sep);
      itemIt = atmTokens.begin();
      if (*itemIt == "AtomExpr") {
        int formCharge = 0;
        ++itemIt;
        // FIX:what if an atom has multiple properties? Seems like they might be
        // separated
        // with ";" ... need to look at that in more detail!
        if ((*itemIt).find("=") == std::string::npos) {
          try {
            formCharge = boost::lexical_cast<int>(*itemIt);
          } catch (boost::bad_lexical_cast &) {
            readNextAtomAttribs = false;
            throw FileParseException("Cannot process mol2 formal charge.");
          }
          // assign the charge
          res->getAtomWithIdx(atomIdx - 1)->setFormalCharge(formCharge);
        }
      }
    }  // endfor
    // check if we have finished reading UNITY_ATOM_ATTR
    if (inStream->eof()) {
      readNextAtomAttribs = false;
    } else {
      tempStr = getLine(inStream);
      if (tempStr == "" || tempStr[0] == '@' || tempStr[0] == '#') {
        readNextAtomAttribs = false;
      }
    }
  }
}

void guessFormalCharges(RWMol *res) {
  // FIX: this whole thing has problems with positively charged pyridines et al.
  for (RWMol::AtomIterator atomIt = res->beginAtoms();
       atomIt != res->endAtoms(); ++atomIt) {
    Atom *at = (*atomIt);
    // assign only if no formal charge set on that atom and atom is not carbon
    // (the latter
    // might needs changing later on - let's see) and not for query atoms (dummy
    // etc.)
    // since this happens only during cleanup of bad substructures
    // FIX: check nitro compounds et al.

    if (at->getFormalCharge() == 0 && at->getSymbol() != "C" &&
        !(at->hasQuery())) {
      // have to calculate the explicit valence without call to
      // calcExplicitValence since
      // this will barf when I have e.g. an uncharged 4 valent nitrogen ...
      int noAromBonds = 0;
      double accum = 0;  // FIX: could this give non int values ?
      for (const auto bnd : res->atomBonds(at)) {
        accum += bnd->getValenceContrib(at);
        if (bnd->getBondType() == Bond::AROMATIC) {
          ++noAromBonds;
        }
      }
      // Assumption: if there is an aromatic bridge atom the accum will be 4.5
      //(three aromatic bonds), e.g. naphthalenes. However those are not charged
      // so we
      // can stop here
      // FIX: a structure such as c12cccc[n+]2cccn1 will not work if we just
      // continue
      // for noAromBonds>2
      // special case checking - if atom c then go on
      if (noAromBonds > 2 && at->getSymbol() == "C") {
        continue;
      }

      // dbtranslate problems - for mols with no UNITY_ATOM_ATTR set - we won't
      // be able to guess
      // if this mol needs to be guessed or not ... hence, we don't guess on
      // atoms with
      // ar specification and not atom type X.ar and in 5-membered ring (there
      // is stuff like
      // c1cc[o+]cc1 that should be charged
      // e.g. 5ring with N.pl3 as NH atom or other atoms without ar
      // specification in aromatic ring
      // FIX: do we need make sure this only happens for atoms in ring?
      if (!res->getRingInfo()->isSssrOrBetter()) {
        MolOps::findSSSR(*res);
      }
      auto tATT = at->getProp<std::string>(common_properties::_TriposAtomType);
      if (at->getIsAromatic() && tATT.find("ar") == std::string::npos &&
          res->getRingInfo()->isAtomInRingOfSize(at->getIdx(), 5)) {
        continue;
      }

      // for dbtranslate a problem will also occur for N.ar with 3 aromatic
      // bonds and no charge
      // assigned for these we don't assign charges (corina will create
      // kekulized input for this
      //(at least in most cases) - anyway, throw a warning!
      if (noAromBonds == 3 && tATT == "N.ar") {
        auto nm = res->getProp<std::string>(common_properties::_Name);
        BOOST_LOG(rdWarningLog)
            << nm
            << ": warning - aromatic N with 3 aromatic bonds - "
               "skipping charge guess for this atom"
            << std::endl;
        continue;
      }

      // sometimes things like benzimidazoles can have only one bond of the
      // imidazole ring as aromatic and the other one as a single bond ...
      // catch that this way - see also the trick from GL
      auto expVal = static_cast<int>(std::round(accum + 0.1));
      const auto &valens =
          PeriodicTable::getTable()->getValenceList(at->getAtomicNum());

      // check default valence and compare to expVal - chg
      // the hypothesis is that we prefer positively charged atoms over
      // negatively charged ones
      // for multi default valence atoms (e.g. CS(O)(O) should end up being
      // C[S+]([O-])[O-] rather
      // than C[S-][O-][O-] but that might change based no different examples
      int nElectrons =
          PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
      int assignChg;
      if (nElectrons >= 4) {
        assignChg = expVal - valens.front();
      } else {
        assignChg = valens.front() - expVal;
      }
      if (assignChg > 0 && nElectrons >= 4) {
        for (auto vi : valens) {
          // Since we do this only for nocharged atoms we can get away without
          // including the charge into this apart from that we do not assign
          // charges higher than +/- 1 for atoms with multiple valence states
          // otherwise the early break would have to go away which in turn would
          // result in things like [S+4] for sulfonamides
          assignChg = expVal - vi;
          if (vi <= expVal && assignChg < 2) {
            break;
          }
        }
      }
      if (assignChg) {
        // no aromatic atom will get aa abs(charge) > 1
        if (at->getIsAromatic() && abs(assignChg) > 1) {
          at->setFormalCharge((assignChg > 0) -
                              (assignChg < 0));  // this results in -1 or +1
        } else {
          at->setFormalCharge(assignChg);
        }
        // corina will create strange nitro groups which look like N(=O)(=O)
        // which in turn will result in
        //[N+2](=O)(=O) this needs to be fixed at this stage since otherwise we
        // would have to check on all
        // N.pl3 during the cleanup substructures step
        if (assignChg == 2 && expVal == 5 && at->getSymbol() == "N") {
          fixNitroSubstructureAndCharge(res, at->getIdx());
        }
        // what if expVal != at->calcExplicitValence();?
        // cannot imagine a case where that will happen now.
      }
    }
  }
}

unsigned int chkNoHNeighbNOx(RWMol *res, ROMol::ADJ_ITER atIdxIt,
                             int &toModIdx) {
  Atom *at = res->getAtomWithIdx(*atIdxIt);
  unsigned int noHNbrs = 0;
  ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
  boost::tie(nbrIdxIt, nbrEndIdxIt) = res->getAtomNeighbors(at);
  while (nbrIdxIt != nbrEndIdxIt) {
    if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 1) {
      ++noHNbrs;
    } else if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
               res->getAtomDegree(res->getAtomWithIdx(*nbrIdxIt)) == 1) {
      // this is a N in an N-oxide constellation
      // we can do the above if clause since mol2 have explicit hydrogens
      toModIdx = *atIdxIt;
    }
    ++nbrIdxIt;
  }
  return noHNbrs;
}

bool cleanUpMol2Substructures(RWMol *res) {
  // NOTE: check the nitro fix in guess formal charges!
  boost::dynamic_bitset<> isFixed(res->getNumAtoms());
  for (auto at : res->atoms()) {
    unsigned int idx = at->getIdx();
    // make sure we haven't finished this atom already
    if (isFixed[idx]) {
      continue;
    }
    auto tAT = at->getProp<std::string>(common_properties::_TriposAtomType);

    if (tAT == "N.4") {
      at->setFormalCharge(1);
    } else if (tAT == "O.co2") {
      // negatively charged carboxylates with O.co2
      // according to Tripos, those should only appear in carboxylates and
      // phosphates,
      if (at->getDegree() != 1) {
        BOOST_LOG(rdWarningLog)
            << "Warning - O.co2 with degree >1." << std::endl;
        return false;
      }
      auto nbrs = res->atomNeighbors(at);
      // this should return only the C.2
      auto nbr = *nbrs.begin();
      auto tATT = nbr->getProp<std::string>(common_properties::_TriposAtomType);
      if (tATT == "P.3") {
        // special case for phosphates
        // we keep the first bond to O.co2 as double and make the rest single
        Bond *b = res->getBondBetweenAtoms(idx, nbr->getIdx());
        b->setBondType(Bond::DOUBLE);
        b->setIsAromatic(false);
        at->setIsAromatic(false);
        isFixed[idx] = 1;
        for (auto onbr : res->atomNeighbors(nbr)) {
          if (onbr->getAtomicNum() == 8 && !isFixed[onbr->getIdx()] &&
              onbr->getProp<std::string>(common_properties::_TriposAtomType) ==
                  "O.co2") {
            Bond *ob = res->getBondBetweenAtoms(nbr->getIdx(), onbr->getIdx());
            ob->setBondType(Bond::SINGLE);
            ob->setIsAromatic(false);
            onbr->setFormalCharge(-1);
            onbr->setIsAromatic(false);
            isFixed[onbr->getIdx()] = 1;
          }
        }
        nbr->setIsAromatic(false);
        isFixed[nbr->getIdx()] = 1;

      } else if (tATT == "C.2" || tATT == "S.o2") {
        // carboxylates and sulfonates
        // this should return only the bond between C.2 and O.co2
        Bond *b = res->getBondBetweenAtoms(idx, nbr->getIdx());
        if (!isFixed[nbr->getIdx()]) {
          // the first occurrence is negatively charged and has a single bond
          b->setBondType(Bond::SINGLE);
          b->setIsAromatic(false);
          at->setFormalCharge(-1);
          at->setIsAromatic(false);
          nbr->setIsAromatic(false);
          isFixed[idx] = 1;
          isFixed[nbr->getIdx()] = 1;
        } else {
          // the other occurrences are not charged and have a double bond
          b->setBondType(Bond::DOUBLE);
          b->setIsAromatic(false);
          at->setIsAromatic(false);
          isFixed[idx] = 1;
        }
      } else {
        std::string nm;
        res->getProp(common_properties::_Name, nm);
        BOOST_LOG(rdWarningLog)
            << nm << ": warning - O.co2 with non C.2 or S.o2 neighbor."
            << std::endl;
        return false;
      }
    } else if (tAT == "C.cat") {
      // positively charged guanidinium groups with C.cat
      // according to Tripos these should only appear in guanidinium groups
      // for the structural fix - the last nitrogen with the least number of
      // heavy atoms will get the double bond and the positive charge.
      // remember : this is not canonical!
      // first - set the C.cat as fixed
      isFixed[idx] = 1;
      ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt, tmpIdxIt;
      unsigned int lowestDeg = 100;
      boost::tie(nbrIdxIt, endNbrsIdxIt) = res->getAtomNeighbors(at);
      // one problem of programs like Corina is, that they will create also
      // C.cat
      // for groups that are not guanidinium. We cannot fix all, but the charged
      // amidine
      // in a ring is taken care of too.
      tmpIdxIt = nbrIdxIt;
      // declare and initialise toModIdx
      int toModIdx = -1;
      unsigned int noNNeighbors = 0;
      while (tmpIdxIt != endNbrsIdxIt) {
        if (res->getAtomWithIdx(*tmpIdxIt)->getSymbol() == "N") {
          ++noNNeighbors;
        }
        ++tmpIdxIt;
      }
      if (noNNeighbors < 2 || noNNeighbors > 3) {
        std::string nm;
        res->getProp(common_properties::_Name, nm);
        BOOST_LOG(rdWarningLog)
            << nm << ": Error - C.Cat with bad number of N neighbors."
            << std::endl;
        return false;
      } else if (noNNeighbors == 2) {
        // the idea is that we assign the positive charge according to the
        // following precedence:
        // 1. is part of N-oxide
        // 2. atom with highest number of hydrogen atoms
        // 3. atom in ring
        // 4. random
        // first we identify the N atoms
        ROMol::ADJ_ITER idxIt1 = nbrIdxIt, idxIt2 = nbrIdxIt;
        bool firstIdent = false;
        while (nbrIdxIt != endNbrsIdxIt) {
          if (res->getAtomWithIdx(*nbrIdxIt)->getSymbol() == "N") {
            // fix the bond to one - only the modified N will have a double bond
            // to C.cat
            res->getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(Bond::SINGLE);
            res->getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
            res->getAtomWithIdx(*nbrIdxIt)->setIsAromatic(false);
            // FIX: what is happening if we hit an atom that was fixed before -
            // probably nothing.
            // since I cannot think of a case where this is a problem - throw a
            // warning
            if (isFixed[*nbrIdxIt]) {
              std::string nm;
              res->getProp(common_properties::_Name, nm);
              BOOST_LOG(rdWarningLog)
                  << nm << ": warning - charged amidine and isFixed atom."
                  << std::endl;
            }
            isFixed[*nbrIdxIt] = 1;
            if (firstIdent) {
              idxIt2 = nbrIdxIt;
            } else {
              idxIt1 = nbrIdxIt;
              firstIdent = true;
            }
          }
          ++nbrIdxIt;
        }
        // now that we know which are the relevant atoms we check the above
        // features
        // is part of N-oxide?
        // number of hydrogens on each neighbour
        unsigned int noHNbrs1 = chkNoHNeighbNOx(res, idxIt1, toModIdx);
        unsigned int noHNbrs2 = chkNoHNeighbNOx(res, idxIt2, toModIdx);
        if (toModIdx < 0) {
          // no N-oxide
          if (noHNbrs1 != noHNbrs2) {
            if (noHNbrs1 > noHNbrs2) {
              toModIdx = *idxIt1;
            } else {
              toModIdx = *idxIt2;  // this is random if both have the same
                                   // number of atoms
            }
          } else {
            // perceive the rings
            if (!res->getRingInfo()->isSssrOrBetter()) {
              MolOps::findSSSR(*res);
            }
            // then we check if both atoms are in a ring
            unsigned int rIdx1 = res->getRingInfo()->numAtomRings((*idxIt1));
            unsigned int rIdx2 = res->getRingInfo()->numAtomRings((*idxIt2));
            if (rIdx1 > rIdx2) {
              toModIdx = *idxIt1;
            } else {
              toModIdx = *idxIt2;
            }
          }
        }
        res->getBondBetweenAtoms(idx, toModIdx)->setBondType(Bond::DOUBLE);
        res->getBondBetweenAtoms(idx, toModIdx)->setIsAromatic(false);
        res->getAtomWithIdx(toModIdx)->setFormalCharge(1);
        at->setIsAromatic(false);
      } else {
        while (nbrIdxIt != endNbrsIdxIt) {
          if (!isFixed[*nbrIdxIt]) {
            // we get in here if this N.pl3 was not seen / fixed before
            Atom *nbr = res->getAtomWithIdx(*nbrIdxIt);
            // get the number of heavy atoms connected to this atom
            ROMol::ADJ_ITER nbrNbrIdxIt, nbrEndNbrsIdxIt;
            unsigned int hvyAtDeg = 0;
            boost::tie(nbrNbrIdxIt, nbrEndNbrsIdxIt) =
                res->getAtomNeighbors(nbr);
            while (nbrNbrIdxIt != nbrEndNbrsIdxIt) {
              if (res->getAtomWithIdx(*nbrNbrIdxIt)->getAtomicNum() > 1) {
                std::string nbrAT;
                res->getAtomWithIdx(*nbrNbrIdxIt)
                    ->getProp(common_properties::_TriposAtomType, nbrAT);
                if (nbrAT == "C.cat") {
                  hvyAtDeg += 2;  // that way we reduce the risk of ionising the
                                  // N attached to another C.cat ...
                } else {
                  ++hvyAtDeg;
                }
              }
              ++nbrNbrIdxIt;
            }
            // now check for lowest heavy atom degree
            if (hvyAtDeg < lowestDeg) {
              toModIdx = *nbrIdxIt;
              lowestDeg = hvyAtDeg;
            }
            // modify the bond between C.Cat and the N.pl3
            Bond *b = res->getBondBetweenAtoms(idx, *nbrIdxIt);
            b->setBondType(Bond::SINGLE);
            b->setIsAromatic(false);
            nbr->setIsAromatic(false);
            // set N.pl3 as fixed
            isFixed[*nbrIdxIt] = 1;
          } else {
            // the N is already fixed - since we don't touch this atom make the
            // bond to single
            // FIX: check on 3-way symmetric guanidinium mol -
            //     this could produce a only single bonded C.cat for bad H mols
            res->getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(Bond::SINGLE);
            res->getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
          }
          ++nbrIdxIt;
        }
        // now modify the respective N and the C.cat
        Bond *b = res->getBondBetweenAtoms(idx, toModIdx);
        b->setBondType(Bond::DOUBLE);
        b->setIsAromatic(false);
        res->getAtomWithIdx(toModIdx)->setFormalCharge(1);
        at->setIsAromatic(false);
      }
    }
    idx++;
  }
  return true;
}

Atom *ParseMol2FileAtomLine(const std::string atomLine, RDGeom::Point3D &pos) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(" \t\n");
  std::string tAN, tAT;
  tokenizer tokens(atomLine, sep);
  tokenizer::const_iterator itemIt = tokens.begin();
  if (itemIt == tokens.end()) {
    throw FileParseException("no info in mol2 atom line");
  }

  auto *res = new Atom();

  // skip TriposAtomId
  ++itemIt;
  if (itemIt == tokens.end()) {
    delete res;
    throw FileParseException("premature end of mol2 atom line");
  }
  // the sybyl atom name does not necessarily make sense - into atom property
  tAN = *itemIt;
  ++itemIt;
  if (itemIt == tokens.end()) {
    delete res;
    throw FileParseException("premature end of mol2 atom line");
  }
  try {
    pos.x = boost::lexical_cast<double>(*itemIt);
    ++itemIt;
    if (itemIt == tokens.end()) {
      delete res;
      throw FileParseException("premature end of mol2 atom line");
    }
    pos.y = boost::lexical_cast<double>(*itemIt);
    ++itemIt;
    if (itemIt == tokens.end()) {
      delete res;
      throw FileParseException("premature end of mol2 atom line");
    }
    pos.z = boost::lexical_cast<double>(*itemIt);
    ++itemIt;
    if (itemIt == tokens.end()) {
      delete res;
      throw FileParseException("premature end of mol2 atom line");
    }
  } catch (boost::bad_lexical_cast &) {
    delete res;
    throw FileParseException("Cannot process mol2 coordinates.");
  }
  // now it becomes interesting - this is the SYBYL atom type. I put this into
  // an atom property and deduce the
  // symbol from that (everything before the dot)
  tAT = *itemIt;
  std::string symb = (*itemIt).substr(0, (*itemIt).find('.'));

  // bad symbols:
  // LP is not an atom so remove it ...
  if (symb == "LP") {
    delete res;
    return nullptr;
  } else if (symb == "ANY" || symb == "Du") {
    // queryAtoms
    // according to the SYBYL spec, these match anything
    auto *query = new QueryAtom(0);
    query->setQuery(makeAtomNullQuery());
    delete res;
    res = query;
  } else if (symb == "HEV") {
    auto *query = new QueryAtom(1);
    query->getQuery()->setNegation(true);
    delete res;
    res = query;
  } else if (symb == "HET") {
    // Tripos: N,O,P,S
    auto *query = new QueryAtom(7);
    query->expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR, true);
    query->expandQuery(makeAtomNumQuery(15), Queries::COMPOSITE_OR, true);
    query->expandQuery(makeAtomNumQuery(16), Queries::COMPOSITE_OR, true);
    delete res;
    res = query;
  } else if (symb == "HAL") {
    // Tripos: F,Cl,Br,I
    auto *query = new QueryAtom(9);
    query->expandQuery(makeAtomNumQuery(17), Queries::COMPOSITE_OR, true);
    query->expandQuery(makeAtomNumQuery(35), Queries::COMPOSITE_OR, true);
    query->expandQuery(makeAtomNumQuery(53), Queries::COMPOSITE_OR, true);
    delete res;
    res = query;
  } else {
    res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
  }

  // now assign the properties
  res->setProp("_TriposAtomName", tAN);  // maybe remove that since it's
                                         // useless?
  res->setProp(common_properties::_TriposAtomType, tAT);
  // no implicit hydrogens for mol2 files
  res->setNoImplicit(true);

  // up to here the fields must be written - check if we have more
  // next comes the subst_id, subst_name - skip those
  ++itemIt;
  if (itemIt == tokens.end()) {
    return res;
  }
  ++itemIt;
  if (itemIt == tokens.end()) {
    return res;
  }
  ++itemIt;
  // the Partial charge in the file
  if (itemIt != tokens.end()) {
    res->setProp("_TriposPartialCharge", *itemIt);
  }
  // we skip the status bit ...

  return res;
}

Bond *ParseMol2FileBondLine(const std::string bondLine,
                            const INT_VECT &idxCorresp) {
  unsigned int idx1, idx2;

  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(" \t\n");

  tokenizer tokens(bondLine, sep);
  tokenizer::const_iterator itemIt = tokens.begin();
  if (itemIt == tokens.end()) {
    throw FileParseException("no info in mol2 bond line");
  }

  try {
    // tripos bond id skip
    ++itemIt;
    if (itemIt == tokens.end()) {
      throw FileParseException("no info in mol2 bond line");
    }
    idx1 = boost::lexical_cast<unsigned int>(*itemIt);
    ++itemIt;
    if (itemIt == tokens.end()) {
      throw FileParseException("no info in mol2 bond line");
    }
    idx2 = boost::lexical_cast<unsigned int>(*itemIt);
    ++itemIt;
    if (itemIt == tokens.end()) {
      throw FileParseException("no info in mol2 bond line");
    }
  } catch (boost::bad_lexical_cast &) {
    throw FileParseException("Cannot process mol2 bonds.");
  }

  // adjust the numbering
  idx1--;
  idx2--;

  // if either of both ends of the bond is not an atom in the mol - return NULL
  if (!(idx1 < idxCorresp.size() || idx2 < idxCorresp.size())) {
    throw FileParseException("index mismatch");
  }

  if (idxCorresp[idx1] < 0 || idxCorresp[idx2] < 0) {
    return nullptr;
  }

  // lexical casts of bond types is not really useful in this case since we
  // could also have strings as values ...
  // use if/else if/end statements for now - maybe change to mapping later - no
  // idea if it's worth it ...?
  Bond::BondType type;
  std::string tBType = *itemIt;
  if (tBType == "1" || tBType == "am") {
    type = Bond::SINGLE;
  } else if (tBType == "2") {
    type = Bond::DOUBLE;
  } else if (tBType == "3") {
    type = Bond::TRIPLE;
  } else if (tBType == "ar") {
    type = Bond::AROMATIC;
  } else if (tBType == "du" || tBType == "un") {
    type = Bond::UNSPECIFIED;  // have to come back here - see if there is any
                               // comment in th documentation ...
  } else {
    // this happens only if some weird thing is in the file or if we encounter a
    // "nc" - not connected ...
    // but why would anyone specify a bond which is not connected?
    BOOST_LOG(rdWarningLog) << "Warning - unsupported bond type: " << tBType
                            << " ignored!" << std::endl;
    return nullptr;
  }
  auto *res = new Bond(type);
  res->setBeginAtomIdx(idxCorresp[idx1]);
  res->setEndAtomIdx(idxCorresp[idx2]);
  return res;
}

void ParseMol2AtomBlock(std::istream *inStream, RWMol *res, unsigned int nAtoms,
                        INT_VECT &idxCorresp) {
  PRECONDITION(inStream, "inStream not valid");
  PRECONDITION(!inStream->eof(), "inStream is at eof");
  PRECONDITION(res, "RWMol not valid");
  PRECONDITION(idxCorresp.size() == nAtoms, "vector size mismatch");
  unsigned int nLP = 0;
  bool hasHAtoms = false;

  std::vector<RDGeom::Point3D> threeDPs;
  threeDPs.reserve(nAtoms);

  for (unsigned int i = 0; i < nAtoms; ++i) {
    std::string tempStr = getLine(inStream);
    if (inStream->eof()) {
      throw FileParseException("premature EOF");
    }
    RDGeom::Point3D pos;
    Atom *atom = ParseMol2FileAtomLine(tempStr, pos);

    // if atom is NULL then we hit LP
    if (atom) {
      int aid = res->addAtom(atom, false, true);
      idxCorresp[i] = aid;
      threeDPs.push_back(pos);
      if (atom->getSymbol() == "H") {
        hasHAtoms = true;
      }
    } else {
      ++nLP;
    }
  }
  // mol2 files need to have hydrogen atoms otherwise formal charge estimation
  // will be problematic
  if (!hasHAtoms) {
    std::string nm;
    res->getProp(common_properties::_Name, nm);
    BOOST_LOG(rdWarningLog) << nm
                            << ": Warning - no explicit hydrogens in "
                               "mol2 file but needed for formal charge "
                               "estimation."
                            << std::endl;
  }
  // create conformer based on 3DPoints and add to RWMol
  auto *conf = new Conformer(nAtoms - nLP);
  std::vector<RDGeom::Point3D>::const_iterator threeDPsIt = threeDPs.begin();
  for (unsigned int i = 0; i < threeDPs.size(); ++i) {
    conf->setAtomPos(i, *threeDPsIt);
    ++threeDPsIt;
  }
  res->addConformer(conf, true);

  POSTCONDITION(res->getNumAtoms() == (nAtoms - nLP),
                "Wrong number of atoms in molecule");
}

void ParseMol2BondBlock(std::istream *inStream, RWMol *res, unsigned int nBonds,
                        const INT_VECT &idxCorresp) {
  PRECONDITION(inStream, "inStream not valid");
  PRECONDITION(!inStream->eof(), "inStream is at eof");
  PRECONDITION(res, "RWMol not valid");
  unsigned int nBadBonds = 0;

  for (unsigned int i = 0; i < nBonds; ++i) {
    std::string tempStr = getLine(inStream);
    if (inStream->eof()) {
      throw FileParseException("premature EOF");
    }
    Bond *bond = ParseMol2FileBondLine(tempStr, idxCorresp);
    // if something weird happened there will be no bond for that line
    if (bond) {
      // if we got an aromatic bond set the flag on the bond and the connected
      // atoms
      if (bond->getBondType() == Bond::AROMATIC) {
        bond->setIsAromatic(true);
        res->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
        res->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
      }
      res->addBond(bond, true);
    } else {
      nBadBonds++;
    }
  }
  POSTCONDITION(res->getNumBonds() == (nBonds - nBadBonds),
                "Wrong number of atoms in molecule");
}

};  // end of anonymous namespace

namespace v2 {
namespace FileParsers {

//------------------------------------------------
//
//  Read a molecule from a stream
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMol2DataStream(std::istream &inStream,
                                             const Mol2ParserParams &params) {
  std::string tempStr, lineBeg;
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(" \t\n");
  Utils::LocaleSwitcher ls;

  // all molecules start with a @<TRIPOS>MOLECULE! There is no other way to
  // define an end of
  // molecule than to find a new one or an eof. Hence I have to read until I
  // find one of the two ...
  std::streampos molStart = 0, atomStart = 0, bondStart = 0, chargeStart = 0;
  while (!inStream.eof() && !inStream.fail()) {
    tempStr = getLine(inStream);
    if (inStream.eof()) {
      break;
    }

    if (tempStr != "" && tempStr[0] == '@') {
      tokenizer tokens(tempStr, sep);
      std::string firstToken = *tokens.begin();
      if (firstToken == "@<TRIPOS>MOLECULE") {
        if (!molStart) {
          // we reach that point in a multimol2 when we have not seen a
          // @<TRIPOS>MOLECULE
          // and set all important flags
          molStart = inStream.tellg();
        } else {
          break;
        }
      } else if (firstToken == "@<TRIPOS>ATOM") {
        atomStart = inStream.tellg();
      } else if (firstToken == "@<TRIPOS>BOND") {
        bondStart = inStream.tellg();
      } else if (firstToken == "@<TRIPOS>UNITY_ATOM_ATTR") {
        // tripos dbtranslate will write not Tripos conform atom types for
        // various aromatic atoms
        // this results in problems with the formal charge guesser. Hence, if we
        // find the
        // UNITY_ATOM_ATTR that contains formal charges we will read those and
        // skip the charger and
        // the substructure cleanup
        chargeStart = inStream.tellg();
      }  // end if seenMolBefore
    }    // end if @
  }      // end while

  // we should reach this point with at least the molStart and atomStart set
  if (!molStart) {
    throw FileParseException("No MOLECULE block found in Mol2 data");
  }
  if (!atomStart) {
    throw FileParseException("No ATOM block found in Mol2 data");
  }

  if (inStream.eof()) {
    inStream.clear();
  }

  inStream.seekg(molStart, std::ios::beg);
  tempStr = getLine(inStream);
  auto res = std::make_unique<RWMol>();
  boost::trim_right(tempStr);
  res->setProp(common_properties::_Name, tempStr);

  tempStr = getLine(inStream);
  tokenizer tokens(tempStr, sep);
  if (tokens.begin() == tokens.end()) {
    throw FileParseException("Empty counts line");
  }

  unsigned int nAtoms = 0, nBonds = 0;
  tokenizer::const_iterator itemIt = tokens.begin();
  // counts line, this is where we really get started
  try {
    nAtoms = boost::lexical_cast<unsigned int>(*itemIt);

    ++itemIt;
    if (itemIt != tokens.end()) {
      nBonds = boost::lexical_cast<unsigned int>(*itemIt);
    }
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert " << *itemIt << " to unsigned int";
    throw FileParseException(errout.str());
  }

  if (nAtoms == 0) {
    throw FileParseException("molecule has no atoms");
  }
  tempStr = getLine(inStream);  // mol_type - ignore
  tempStr = getLine(inStream);
  boost::trim(tempStr);
  res->setProp("_TriposChargeType", tempStr);
  // stop here since we don't support anything else from the MOLECULE block
  INT_VECT idxCorresp(nAtoms, -1);
  inStream.seekg(atomStart, std::ios::beg);
  ParseMol2AtomBlock(&inStream, res.get(), nAtoms, idxCorresp);
  if (nBonds) {
    // stop here since we don't support anything else from the MOLECULE block
    inStream.seekg(bondStart, std::ios::beg);
    ParseMol2BondBlock(&inStream, res.get(), nBonds, idxCorresp);
  }

  if (!chargeStart) {
    bool molFixed;
    if (params.cleanupSubstructures) {
      molFixed = cleanUpMol2Substructures(res.get());
    } else {
      molFixed = true;
    }

    if (!molFixed) {
      return nullptr;
    }

    // mol2 format does not support formal charge information, hence we need to
    // guess it based on default and explicit valences
    guessFormalCharges(res.get());
  } else {
    inStream.seekg(chargeStart, std::ios::beg);
    readFormalChargesFromAttr(&inStream, res.get());
  }

  // set chirality prior to sanitization since it happens from 3D and it's not
  // possible anymore once the hydrogens are removed
  // FIX: for now this is only for the first conformer - need to be changed once
  // we use multiconformer files
  MolOps::assignChiralTypesFrom3D(*res);

  if (res && params.sanitize) {
    MolOps::cleanUp(*res);

    try {
      if (params.removeHs) {
        // Bond stereo detection must happen before H removal, or
        // else we might be removing stereogenic H atoms in double
        // bonds (e.g. imines). But before we run stereo detection,
        // we need to run mol cleanup so don't have trouble with
        // e.g. nitro groups. Sadly, this a;; means we will find
        // run both cleanup and ring finding twice (a fast find
        // rings in bond stereo detection, and another in
        // sanitization's SSSR symmetrization).
        unsigned int failedOp = 0;
        MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
        MolOps::detectBondStereochemistry(*res);
        MolOps::removeHs(*res, false, false);
      } else {
        MolOps::sanitizeMol(*res);
        MolOps::detectBondStereochemistry(*res);
      }

    } catch (MolSanitizeException &se) {
      BOOST_LOG(rdWarningLog) << "sanitize ";
      std::string molName;
      res->getProp(common_properties::_Name, molName);
      BOOST_LOG(rdWarningLog) << molName << ": ";
      throw se;
    }

    res->updatePropertyCache(false);
    MolOps::assignStereochemistry(*res, true, true);
  }

  return res;
};

//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMol2Block(const std::string &molBlock,
                                        const Mol2ParserParams &params) {
  std::istringstream inStream(molBlock);
  return MolFromMol2DataStream(inStream, params);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMol2File(const std::string &fName,
                                       const Mol2ParserParams &params) {
  // FIX: this binary mode of opening file is here because of a bug in VC++ 6.0
  // the function "tellg" does not work correctly if we do not open it this way
  //   Jan 2009: Confirmed that this is still the case in visual studio 2008
  std::ifstream inStream(fName.c_str(), std::ios_base::binary);
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    return MolFromMol2DataStream(inStream, params);
  } else {
    return nullptr;
  }
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
