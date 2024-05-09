//
//  Copyright (c) 2011-2022 Novartis Institutes for BioMedical Research Inc. and
//  other RDkit contributors
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

//
// Known issues:
//
// - Allene stereochemistry is not processed
//
// - advanced InChI features - such as fixed-H layer - have not been tested
//
// - InChI-write issues on broken molecules (e.g. PubChem Compound 42622894,
// 42622893, 42620342, 42621905, 42622374, 42617647), because RDKit and standard
// InChI binary will "fix" them differently. However, if the molecules have been
// preprocessed by RDKit, then most have no write issue.
//
// - InChI-read issues on molecules with metals.
//
// - For molecules with large ring and no coordinates, RDKit does not provide
// sufficient ring stereochemistry required by InChI and will result in a
// warning about undefined stereo. InChI requires all single bond in a ring with
// 8 or more bonds to have E/Z parity assigned. If coordinates are provided,
// then InChI will infer stereochemistry from them.
//
// - Radical electrons messed up by InChI are not repaired. One example molecule
// is PubChem Compound 10784031, 10784032
//
#include <string>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <inchi_api.h>
#include <cstring>
#include <vector>
#include <stack>
#include <set>
#include <queue>
#include "inchi.h"
#include <algorithm>

#include <RDGeneral/BoostStartInclude.h>
#include <tuple>
#include <RDGeneral/BoostEndInclude.h>

//#define DEBUG 1
namespace RDKit {
namespace {
/* assignBondDirs
 * assign bond direction for neighboring bonds of stereo double bonds based
 * on two sets of constraints: zBondPairs gives the pairs of bonds that must
 * have the same direction and eBondPairs gives the pairs of bonds that must
 * have different directions
 *
 * return true on success and false when it is not doable
 */
typedef std::pair<int, int> INT_PAIR;
typedef std::vector<INT_PAIR> INT_PAIR_VECT;
bool assignBondDirs(RWMol& mol, INT_PAIR_VECT& zBondPairs,
                    INT_PAIR_VECT& eBondPairs) {
  // bonds to assign
  std::set<int> pending;
  for (const auto& pair : zBondPairs) {
    pending.insert(pair.first);
    pending.insert(pair.second);
  }
  for (const auto& pair : eBondPairs) {
    pending.insert(pair.first);
    pending.insert(pair.second);
  }
  // a queue for pending assignments
  typedef std::queue<std::pair<int, Bond::BondDir>> ASSIGNMENTQTYPE;
  ASSIGNMENTQTYPE queue;
  // in a loop, modify one bond at a time, until all bonds are assigned
  while (!pending.empty() || !queue.empty()) {
    if (queue.empty()) {
      // pumping one bond from pending to queue
      queue.push(std::make_pair(*(pending.begin()), Bond::ENDUPRIGHT));
    } else {
      // pop one entry from queue and do the actual assignment
      int curBondIdx;
      Bond::BondDir dir;
      boost::tie(curBondIdx, dir) = queue.front();
      queue.pop();
      Bond* bond = mol.getBondWithIdx(curBondIdx);
      // is it assigned already?
      if (bond->getBondDir() != Bond::NONE) {
        // assigned. then check conflict
        if (bond->getBondDir() != dir) {
          // not doable
          return false;
        }
      } else {
        // assign since it's not assigned yet
        bond->setBondDir(dir);
        auto searchItr = pending.find(curBondIdx);
        if (searchItr != pending.end()) {
          pending.erase(searchItr);
        }
        // find all affecting bonds and add to queue by going thru all rules
        Bond::BondDir otherDir =
            dir == Bond::ENDUPRIGHT ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT;
        // same routine for zBondPairs and eBondPairs
        // use a switch _ to go through both by setting _ to 0 and then 1
        for (int _ = 0; _ < 2; _++) {
          INT_PAIR_VECT* _rules = _ == 0 ? &zBondPairs : &eBondPairs;
          Bond::BondDir _dir = _ == 0 ? dir : otherDir;
          for (const auto& pair : *_rules) {
            int other = -1;
            if (pair.first == curBondIdx) {
              other = pair.second;
            } else if (pair.second == curBondIdx) {
              other = pair.first;
            }
            // a match?
            if (other != curBondIdx && other != -1) {
              Bond* otherBond = mol.getBondWithIdx(other);
              // check if it is assigned
              if (otherBond->getBondDir() != Bond::NONE) {
                // assigned. check conflict
                if (otherBond->getBondDir() != _dir) {
                  // not doable
                  return false;
                }
              } else {
                // not assigned, then add to queue
                queue.push(std::make_pair(otherBond->getIdx(), _dir));
              }  // end if otherBond's bond direction check
            }    // end if there is a match
          }      // end loop over pairs in _rules
        }        // end for _ to go thru rule sets
      }          // end if this bond is assigned
    }            // end if queue is empty
  }              // end while on pending set and queue
  return true;
}

/* findAlternatingBonds
 *
 * This is a modified DFS that returns the shortest path consisting of
 * alternating bonds from the current node to a node with desired atomic
 * number.
 *
 * The DFS uses a static variable to remember which nodes have already been
 * visited and therefore is not threadsafe.
 *
 * The traversal is done recursively. At any point of the traversal, one
 * single copy of <path> is maintained. If the desired atom has not been
 * found, <path> is empty. If it has been found for once, <path> maintains
 * the path between the desired atom to the lowest common ancestor between
 * the desired atom and the current node being visited. If it is found for a
 * second time, the shortest path will survive. <path> always maintain the
 * suffix of the final search result with member bonds in reverse order.
 * This is doable because the call stack implicitly keeps track of the path,
 * and we just reproduce the path through tracing back the call stack.
 *
 * The return value of the function is a pointer value that is either NULL
 * if it could find a better path or points to a target atom if it is able
 * to do better than the best-so-far result before it is called. At each
 * point of the traversal of the search tree, the function asks the subtree
 * rooted at the current node whether they could enhance the current best
 * <path>. If a subtree answers yes (and returns a non-NULL pointer), then
 * the <path> value has been updated, and this call should push itself to
 * the <path> and return the non-NULL pointer to its caller. Otherwise, it
 * should signal its caller that it cannot enhance <path> by returning a
 * NULL pointer.
 *
 * If maxPathLength is larger than 2, than we are looking for a path with
 * alternating single and double bond. If maxPathLength is 2, then it's
 * basically a path with desiredNextBondType and then a
 * desiredEndingBondType. If maxPathLength is 1, you are looking at
 * immediate neighbor and desiredNextBondType and desiredEndingBondType must
 * be the same.
 */
Atom* findAlternatingBonds(
    ROMol& mol, Atom* current, int desiredAtomicNumber, int desiredAtomCharge,
    Bond::BondType desiredNextBondType, Bond::BondType desiredEndingBondType,
    unsigned int currentPathLength, unsigned int maxPathLength, Bond* lastBond,
    /*OUT*/ std::stack<Bond*>& path, std::set<int>& _visited) {
  // memory for what has been visited
  if (lastBond == nullptr) {
    _visited.clear();
    while (!path.empty()) {
      path.pop();
    }
  }
  _visited.insert(current->getIdx());

  // for (int i = 0; i < currentPathLength; i ++)
  //  std::cerr << ".";
  // std::cerr << (int) current->getIdx() << "("
  //  << (int) current->getAtomicNum()
  //  << ")" << std::endl;

  // is this atom the desired one?
  if (lastBond && current->getAtomicNum() == desiredAtomicNumber &&
      lastBond->getBondType() == desiredEndingBondType &&
      current->getFormalCharge() == desiredAtomCharge) {
    // Yes! But am I better than the existing one - if one exists?
    if (path.size() == 0 || path.size() > currentPathLength) {
      // Yes! clear the path and repopulate it
      while (!path.empty()) {
        path.pop();
      }
      // add myself to the path
      path.push(lastBond);
      return current;
    } else {
      // I am no better than the existing one. This will also cause the
      // path search to not continue down
      return nullptr;
    }
  }

  // searching too far?
  if (maxPathLength <= currentPathLength) {
    return nullptr;
  }

  // continue searching down
  RWMol::ADJ_ITER nid, end;
  Atom *target = nullptr, *temp;
  for (boost::tie(nid, end) = mol.getAtomNeighbors(current); nid != end;
       nid++) {
    if (_visited.find(*nid) != _visited.end()) {
      continue;
    }
    // check whether bond is valid for search to go down through it
    Bond* bond = mol.getBondBetweenAtoms(current->getIdx(), *nid);
    if (bond->getBondType() == desiredNextBondType) {
      // recursive call: for all ways to extend the path, ask each to try
      // enhancing the current best path (stored in <path>)
      // by setting SINGLE as the default, we allow a very special case to
      // be supported: a TRIPLE bond followed by a SINGLE bond
      // This is used in _Valence5NCleanUp2
      Bond::BondType nextBondType = Bond::SINGLE;
      if (desiredNextBondType == Bond::SINGLE) {
        nextBondType = Bond::DOUBLE;
      }
      if ((temp = findAlternatingBonds(
               mol, mol.getAtomWithIdx(*nid), desiredAtomicNumber,
               desiredAtomCharge, nextBondType, desiredEndingBondType,
               currentPathLength + 1, maxPathLength, bond, path, _visited)) !=
          nullptr) {
        target = temp;
      }
    } else if (desiredEndingBondType != Bond::SINGLE &&
               desiredEndingBondType != Bond::DOUBLE &&
               bond->getBondType() == desiredEndingBondType) {
      // try recursive call limited to one level down to see whether
      // this can serve as the last leg of the path. This is done only if
      // the desiredEndingBondType is not part of the alternating bonds
      if ((temp = findAlternatingBonds(
               mol, mol.getAtomWithIdx(*nid), desiredAtomicNumber,
               desiredAtomCharge, Bond::UNSPECIFIED, /* no next */
               desiredEndingBondType, currentPathLength + 1,
               0, /* this limits the recursion */
               bond, path, _visited))) {
        target = temp;
      }
    }
  }

  // about the return
  if (target != nullptr) {
    if (lastBond) {
      path.push(lastBond);
    }
    return target;
  }
  return nullptr;
}

int getNumDoubleBondedNegativelyChargedNeighboringSi(ROMol& mol, Atom* a) {
  RWMol::ADJ_ITER nid1, end1;
  boost::tie(nid1, end1) = mol.getAtomNeighbors(a);
  int nSi = 0;
  int thisId = a->getIdx();
  while (nid1 != end1) {
    Atom* nbr = mol.getAtomWithIdx(*nid1);
    Bond* bond = mol.getBondBetweenAtoms(*nid1, thisId);
    if (nbr->getAtomicNum() == 14 && nbr->getFormalCharge() == -1 &&
        bond->getBondType() == Bond::DOUBLE) {
      nSi++;
    }
    nid1++;
  }
  return nSi;
}

// clean C1=NN=[N-]=N1
bool _Valence4NCleanUp1(RWMol& mol, Atom* atom) {
  // replace the N- with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != -1 ||
      atom->calcExplicitValence(false) != 4) {
    return false;
  }
  atom->setAtomicNum(50);
  atom->setFormalCharge(0);

  // substructure matching
  auto* query = new RWMol();
  query->addAtom(new Atom(6), false, true);   // 0
  query->addAtom(new Atom(7), false, true);   // 1
  query->addAtom(new Atom(50), false, true);  // 2
  query->addAtom(new Atom(7), false, true);   // 3
  query->addAtom(new Atom(7), false, true);   // 4
  query->addBond(0, 1, Bond::SINGLE);
  query->addBond(1, 2, Bond::DOUBLE);
  query->addBond(2, 3, Bond::DOUBLE);
  query->addBond(3, 4, Bond::SINGLE);
  query->addBond(4, 0, Bond::DOUBLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;
  // no action if none or more than one match was found
  if (fgpMatches.size() != 1) {
    atom->setAtomicNum(7);
    atom->setFormalCharge(-1);
    return false;
  }

  // collect matching atoms
  int map[5];
  MatchVectType match = fgpMatches[0];
  for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
       mi++) {
    map[mi->first] = mi->second;
  }
  // flip bonds
  mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[2], map[3])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[3], map[4])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[4], map[0])->setBondType(Bond::SINGLE);
  // change the problematic N-
  atom->setAtomicNum(7);
  atom->setFormalCharge(-1);
  return true;
}

// directly to a N via double bond
bool _Valence4NCleanUp2(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }

  stack.top()->setBondType(Bond::SINGLE);
  atom->setFormalCharge(0);
  target->setFormalCharge(-1);
  return true;
}

// try search for valence-5 N connected to a N+
bool _Valence5NCleanUp1(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 1, Bond::DOUBLE, Bond::DOUBLE, 0, 5,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }
  target->setFormalCharge(0);
  target->calcExplicitValence(false);
  while (!stack.empty()) {
    if (stack.top()->getBondType() == Bond::DOUBLE) {
      stack.top()->setBondType(Bond::SINGLE);
    } else {
      stack.top()->setBondType(Bond::DOUBLE);
    }
    stack.pop();
  }
  atom->setFormalCharge(1);
  return true;
}

// N connected to N- through a tiple then single bond
bool _Valence5NCleanUp2(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, -1, Bond::TRIPLE, Bond::SINGLE, 0, 2,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }

  Bond* bond = stack.top();
  bond->setBondType(Bond::SINGLE);
  if (bond->getBeginAtomIdx() == atom->getIdx()) {
    mol.getAtomWithIdx(bond->getEndAtomIdx())->setFormalCharge(-1);
  } else {
    mol.getAtomWithIdx(bond->getBeginAtomIdx())->setFormalCharge(-1);
  }
  stack.pop();
  stack.top()->setBondType(Bond::DOUBLE);
  target->setFormalCharge(0);
  target->calcExplicitValence(false);
  atom->calcExplicitValence(false);
  return true;
}

// directly to a N via double bond
bool _Valence5NCleanUp3(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }

  // we are double bonded to a neighboring N. Check to see if we are also
  // double bonded to an O. If so, we don't want to mess with the other N
  // this occurs because the InChI code produces this structure:
  //  CN(=O)=N(=O)C
  // and we don't want to mess with that.
  // this was github #1572

  std::stack<Bond*> stack2;
  std::set<int> _visited2;
  Atom* target2 =
      findAlternatingBonds(mol, atom, 8, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
                           nullptr, stack2, _visited2);
  if (target2 == nullptr) {
    target->setFormalCharge(-1);
    target->calcExplicitValence(false);
    stack.top()->setBondType(Bond::SINGLE);
    atom->setFormalCharge(1);
    atom->calcExplicitValence(false);
  }
  return true;
}

// N connected to two Si- via double bonds; also a positive charged S
// connected to a non-charged C. shift the charge to the C
bool _Valence5NCleanUp4(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  RWMol::ADJ_ITER nid1, end1;
  int nSi = 0;
  int thisId = atom->getIdx();
  Atom* nbrs[2];
  Bond* bonds[2];
  boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
  while (nid1 != end1) {
    Atom* nbr = mol.getAtomWithIdx(*nid1);
    Bond* bond = mol.getBondBetweenAtoms(*nid1, thisId);
    if (nbr->getAtomicNum() == 14 && nbr->getFormalCharge() == -1 &&
        bond->getBondType() == Bond::DOUBLE) {
      if (nSi >= 2) {
        return false;
      }
      nbrs[nSi] = nbr;
      bonds[nSi] = bond;
      nSi++;
    }
    ++nid1;
  }
  if (nSi != 2) {
    return false;
  }
  nbrs[0]->setFormalCharge(0);
  nbrs[1]->setFormalCharge(0);
  bonds[0]->setBondType(Bond::SINGLE);
  bonds[1]->setBondType(Bond::SINGLE);

#if 0
      // FIX
      // not clear why this is here, but it almost definitely shouldn't be
      Atom* s = NULL;
      Atom* c = NULL;
      Bond* sc_bond;
      ROMol::VERTEX_ITER atBegin,atEnd;
      boost::tie(atBegin,atEnd) = mol.getVertices();
      while (atBegin != atEnd) {
          ATOM_SPTR at2 = mol[*atBegin];
          if (at2->getAtomicNum() == 16 && at2->getFormalCharge() == 1) {
            boost::tie(nid1, end1) = mol.getAtomNeighbors(at2);
           while (nid1 != end1) {
             Atom* nbr = mol.getAtomWithIdx(*nid1);
             Bond* bond = mol.getBondBetweenAtoms(*nid1, at2->getIdx());
             if (nbr->getAtomicNum() == 6 && nbr->getFormalCharge() == 0 &&
                 bond->getBondType() == Bond::DOUBLE) {
               s = &(*at2);
               c = nbr;
               sc_bond = bond;
               break;
             }
             ++nid1
           }
          }
          ++atBegin;
      }

      if (s == NULL) return false;
      s->setFormalCharge(0);
      c->setFormalCharge(-1);
      sc_bond->setBondType(Bond::SINGLE);
      atom->setFormalCharge(0);
#endif
  return true;
}

bool _Valence5NCleanUp5(RWMol& mol, Atom* atom, int atomicNum) {
  PRECONDITION(
      atomicNum == 8 || atomicNum == 16 || atomicNum == 9 || atomicNum == 17,
      "this cleanup looks for O or S or Cl or F");
  std::stack<Bond*> stackCharged, stackUncharged, *stack;
  // try search for valence-5 N connected to O or S, determined by the
  // <atomicNum> parameter with alternating
  // bonds if there is a charged Oxygen and an uncharged one both
  // connected to our N through alternating bonds, strip the charge
  // and hydrogen from the charged one, and the use the uncharged
  // one in our procedure
  // see InChI for PubChem compound 10775236:
  //   CC(C1=CC=CC=N1=C2C(OC)=O)CC2=[OH+]
  // is converted into
  //   COC(O)=C1[n+]2ccccc2C(C)CC1=O
  Atom *unchargedOxygen, *chargedOxygen;
  std::set<int> _visited;
  unchargedOxygen =
      findAlternatingBonds(mol, atom, atomicNum, 0, Bond::DOUBLE, Bond::DOUBLE,
                           0, 7, nullptr, stackUncharged, _visited);
  chargedOxygen =
      findAlternatingBonds(mol, atom, atomicNum, 1, Bond::DOUBLE, Bond::DOUBLE,
                           0, 7, nullptr, stackCharged, _visited);
  if (unchargedOxygen == nullptr && chargedOxygen == nullptr) {
    return false;
  }

  stack = &stackUncharged;
  if (unchargedOxygen == nullptr) {
    stack = &stackCharged;
  }
  if (unchargedOxygen && chargedOxygen) {
    // both exists. fix the charged oxygen now by set it to neutral
    // with its hydrogen taken and moved later to the uncharged one
    CHECK_INVARIANT(chargedOxygen->getFormalCharge() == 1,
                    "expecting +1 charge");
    chargedOxygen->setFormalCharge(0);
    chargedOxygen->setNumExplicitHs(0);  // this hydrogen will be
                                         // added to the uncharged
                                         // oxygen later
  }
  if (unchargedOxygen || chargedOxygen) {
    // set charge on N
    atom->setFormalCharge(1);
    // switch all bonds
    Bond* b;
    while (!stack->empty()) {
      b = stack->top();
      if (b->getBondType() == Bond::DOUBLE) {
        b->setBondType(Bond::SINGLE);
      } else {
        b->setBondType(Bond::DOUBLE);
      }
      stack->pop();
    }
    if (unchargedOxygen && chargedOxygen) {
      // both charged and uncharged oxygen are found, the uncharged
      // remains uncharged and take the hydrogen from the charged
      // one
      unchargedOxygen->setNumExplicitHs(1);
    } else if (unchargedOxygen) {
      // if only uncharged oxygen is found, not the oxygen has -1
      // charge
      unchargedOxygen->setFormalCharge(-1);
    } else {
      // if only charged oxygen is found, it's neutral now (and
      // keeps its hydrogen)
      chargedOxygen->setFormalCharge(0);
    }
    if (chargedOxygen) {
      chargedOxygen->calcExplicitValence(false);
    }
    if (unchargedOxygen) {
      unchargedOxygen->calcExplicitValence(false);
    }
  }
  return true;
}

// clean CN1=CCN=CC=1
// example: PubChem 10781979
bool _Valence5NCleanUp6(RWMol& mol, Atom* atom) {
  // replace the N with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
      atom->calcExplicitValence(false) != 5) {
    return false;
  }
  atom->setAtomicNum(50);

  // substructure matching
  auto* query = new RWMol();
  query->addAtom(new Atom(6), false, true);   // 0
  query->addAtom(new Atom(6), false, true);   // 1
  query->addAtom(new Atom(50), false, true);  // 2
  query->addAtom(new Atom(6), false, true);   // 3
  query->addAtom(new Atom(6), false, true);   // 4
  query->addAtom(new Atom(7), false, true);   // 5
  query->addAtom(new Atom(6), false, true);   // 6
  query->addBond(0, 1, Bond::SINGLE);
  query->addBond(1, 2, Bond::DOUBLE);
  query->addBond(2, 3, Bond::DOUBLE);
  query->addBond(3, 4, Bond::UNSPECIFIED);
  query->addBond(4, 5, Bond::SINGLE);
  query->addBond(5, 0, Bond::DOUBLE);
  query->addBond(2, 6, Bond::SINGLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;
  // no action if none or more than one match was found
  if (fgpMatches.size() != 1) {
    atom->setAtomicNum(7);
    return false;
  }

  // collect matching atoms
  int map[7];
  MatchVectType match = fgpMatches[0];
  for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
       mi++) {
    map[mi->first] = mi->second;
  }
  // flip bonds
  mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[4], map[5])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
  // change the problematic N
  atom->setAtomicNum(7);
  atom->setFormalCharge(1);
  return true;
}

// clean CN1=NCOCC=1 that is connected via alternating bonds to O
// example: PubChem 10781979
bool _Valence5NCleanUp7(RWMol& mol, Atom* atom) {
  // is it connected to O via alternating bonds?
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 8, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 5,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }
  // replace the N with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
      atom->calcExplicitValence(false) != 5) {
    return false;
  }
  atom->setAtomicNum(50);

  // substructure matching
  auto* query = new RWMol();
  query->addAtom(new Atom(6), false, true);   // 0
  query->addAtom(new Atom(6), false, true);   // 1
  query->addAtom(new Atom(50), false, true);  // 2
  query->addAtom(new Atom(7), false, true);   // 3
  query->addAtom(new Atom(6), false, true);   // 4
  query->addAtom(new Atom(8), false, true);   // 5
  query->addAtom(new Atom(6), false, true);   // 6
  query->addBond(0, 1, Bond::UNSPECIFIED);
  query->addBond(1, 2, Bond::DOUBLE);
  query->addBond(2, 3, Bond::DOUBLE);
  query->addBond(3, 4, Bond::SINGLE);
  query->addBond(4, 5, Bond::SINGLE);
  query->addBond(5, 0, Bond::SINGLE);
  query->addBond(2, 6, Bond::SINGLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;
  // no action if none or more than one match was found
  if (fgpMatches.size() != 1) {
    atom->setAtomicNum(7);
    return false;
  }

  // collect matching atoms
  int map[7];
  MatchVectType match = fgpMatches[0];
  for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
       mi++) {
    map[mi->first] = mi->second;
  }
  // flip bonds
  mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
  Bond* b;
  while (!stack.empty()) {
    b = stack.top();
    if (b->getBondType() == Bond::DOUBLE) {
      b->setBondType(Bond::SINGLE);
    } else {
      b->setBondType(Bond::DOUBLE);
    }
    stack.pop();
  }
  // set charge on oxygen
  target->setFormalCharge(-1);
  // change the problematic N
  atom->setAtomicNum(7);
  return true;
}

// clean [N]=C1N=CN=N1
// example: PubChem 10782655
bool _Valence5NCleanUp8(RWMol& mol, Atom* atom) {
  // replace the N with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
      atom->calcExplicitValence(false) != 5) {
    return false;
  }
  atom->setAtomicNum(50);

  // substructure matching
  auto* query = new RWMol();
  query->addAtom(new Atom(6), false, true);   // 0
  query->addAtom(new Atom(7), false, true);   // 1
  query->addAtom(new Atom(6), false, true);   // 2
  query->addAtom(new Atom(7), false, true);   // 3
  query->addAtom(new Atom(7), false, true);   // 4
  query->addAtom(new Atom(50), false, true);  // 5
  query->addBond(0, 1, Bond::SINGLE);
  query->addBond(1, 2, Bond::DOUBLE);
  query->addBond(2, 3, Bond::SINGLE);
  query->addBond(3, 4, Bond::DOUBLE);
  query->addBond(4, 0, Bond::SINGLE);
  query->addBond(5, 0, Bond::DOUBLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;

  if (fgpMatches.size() != 1) {
    atom->setAtomicNum(7);
    return false;
  }

  // collect matching atoms
  int map[6];
  MatchVectType match = fgpMatches[0];
  for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
       mi++) {
    map[mi->first] = mi->second;
  }
  // flip bonds
  mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[2], map[3])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[3], map[4])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[4], map[0])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
  mol.getAtomWithIdx(map[1])->setFormalCharge(-1);
  // change the problematic N
  atom->setAtomicNum(7);
  atom->setFormalCharge(1);
  return true;
}

// clean [N]=C1C=CN=N1
// example: PubChem 10785993
bool _Valence5NCleanUp9(RWMol& mol, Atom* atom) {
  // replace the N with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
      atom->calcExplicitValence(false) != 5) {
    return false;
  }
  atom->setAtomicNum(50);

  // substructure matching
  auto* query = new RWMol();
  query->addAtom(new Atom(6), false, true);   // 0
  query->addAtom(new Atom(7), false, true);   // 1
  query->addAtom(new Atom(7), false, true);   // 2
  query->addAtom(new Atom(6), false, true);   // 3
  query->addAtom(new Atom(6), false, true);   // 4
  query->addAtom(new Atom(50), false, true);  // 5
  query->addBond(0, 1, Bond::SINGLE);
  query->addBond(1, 2, Bond::DOUBLE);
  query->addBond(2, 3, Bond::SINGLE);
  query->addBond(3, 4, Bond::DOUBLE);
  query->addBond(4, 0, Bond::SINGLE);
  query->addBond(5, 0, Bond::DOUBLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;

  if (fgpMatches.size() != 1) {
    atom->setAtomicNum(7);
    return false;
  }

  // collect matching atoms
  int map[6];
  MatchVectType match = fgpMatches[0];
  for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
       mi++) {
    map[mi->first] = mi->second;
  }
  // flip bonds
  mol.getBondBetweenAtoms(map[0], map[1])->setBondType(Bond::DOUBLE);
  mol.getBondBetweenAtoms(map[1], map[2])->setBondType(Bond::SINGLE);
  mol.getBondBetweenAtoms(map[5], map[0])->setBondType(Bond::SINGLE);
  mol.getAtomWithIdx(map[2])->setFormalCharge(-1);
  // change the problematic N
  atom->setAtomicNum(7);
  atom->setFormalCharge(1);
  return true;
}

// N connected via alternating bonds to N=N
bool _Valence5NCleanUpA(RWMol& mol, Atom* atom) {
  // replace the N with Sn
  if (atom->getAtomicNum() != 7 || atom->getFormalCharge() != 0 ||
      atom->calcExplicitValence(false) != 5) {
    return false;
  }
  // first find the N=N
  auto* query = new RWMol();
  query->addAtom(new Atom(7), false, true);  // 0
  query->addAtom(new Atom(7), false, true);  // 1
  query->addBond(0, 1, Bond::DOUBLE);

  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *query, fgpMatches);
  delete query;

  if (fgpMatches.size() == 0) {
    return false;
  }

  std::stack<Bond*> bestPath;
  for (const auto& match : fgpMatches) {
    // does the match contains the current atom?
    if (match[0].second == static_cast<int>(atom->getIdx()) ||
        match[1].second == static_cast<int>(atom->getIdx())) {
      continue;
    }
    // set both matched N to Sn
    mol.getAtomWithIdx(match[0].second)->setAtomicNum(50);
    mol.getAtomWithIdx(match[1].second)->setAtomicNum(50);
    // now search the path from current atom to these atoms
    std::stack<Bond*> stack;
    std::set<int> _visited;
    Atom* target =
        findAlternatingBonds(mol, atom, 50, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 9,
                             nullptr, stack, _visited);
    if (target && (bestPath.empty() || stack.size() < bestPath.size())) {
      bestPath = stack;
    }
    mol.getAtomWithIdx(match[0].second)->setAtomicNum(7);
    mol.getAtomWithIdx(match[1].second)->setAtomicNum(7);
  }

  if (!bestPath.empty()) {
    while (!bestPath.empty()) {
      Bond* bond = bestPath.top();
      if (bond->getBondType() == Bond::SINGLE) {
        bond->setBondType(Bond::DOUBLE);
      } else {
        bond->setBondType(Bond::SINGLE);
      }
      bestPath.pop();
    }
    atom->setFormalCharge(1);
    atom->calcExplicitValence(false);
    return true;
  }
  return false;
}
//
// directly to a C via double bond; this is last resort
bool _Valence5NCleanUpB(RWMol& mol, Atom* atom) {
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 6, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
                           nullptr, stack, _visited);
  if (target == nullptr) {
    return false;
  }

  target->setFormalCharge(-1);
  target->calcExplicitValence(false);
  stack.top()->setBondType(Bond::SINGLE);
  atom->setFormalCharge(1);
  atom->calcExplicitValence(false);
  return true;
}

//   C([S-](=O)(=O)=O)
// to:
//   C(S([O-])(=O)=O)
// for instance:
//   CC(C)(C1=CC([S-](=O)(=O)=O)=[N+](F)C=C1)C
// is converted to:
//   CC(C)(C1=CC(S[O-](=O)=O)=[N+](F)C=C1)C
bool _Valence7SCleanUp1(RWMol& mol, Atom* atom) {
  if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
      atom->calcExplicitValence(false) != 7) {
    return false;
  }
  int aid = atom->getIdx();
  int neighborsC = 0;
  int neighborsO = 0;
  RWMol::ADJ_ITER nid, nid1, end1;
  boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
  nid = end1;
  while (nid1 != end1) {
    Atom* otherAtom = mol.getAtomWithIdx(*nid1);
    if (otherAtom->getAtomicNum() == 8) {
      if (mol.getBondBetweenAtoms(*nid1, aid)->getBondType() != Bond::DOUBLE) {
        neighborsO = 100;
        break;
      } else {
        nid = nid1;
        neighborsO++;
      }
    } else if (otherAtom->getAtomicNum() == 6) {
      if (mol.getBondBetweenAtoms(*nid1, aid)->getBondType() != Bond::SINGLE) {
        neighborsC = 100;
        break;
      } else {
        neighborsC++;
      }
    } else {
      neighborsC = 100;
      break;
    }
    nid1++;
  }
  if (nid != end1 && (neighborsC == 1 || neighborsO == 3)) {
    mol.getBondBetweenAtoms(*nid, aid)->setBondType(Bond::SINGLE);
    Atom* otherAtom = mol.getAtomWithIdx(*nid);
    otherAtom->setFormalCharge(-1);
    atom->setFormalCharge(0);
    otherAtom->calcExplicitValence(false);
    atom->calcExplicitValence(false);
    return true;
  } else {
    return false;
  }
}

// [S-]=CC#N
bool _Valence7SCleanUp2(RWMol& mol, Atom* atom) {
  if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
      atom->calcExplicitValence(false) != 7) {
    return false;
  }

  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::TRIPLE, 0, 3,
                           nullptr, stack, _visited);
  if (target) {
    while (!stack.empty()) {
      Bond* bond = stack.top();
      if (bond->getBondType() == Bond::SINGLE) {
        bond->setBondType(Bond::DOUBLE);
      } else if (bond->getBondType() == Bond::DOUBLE) {
        bond->setBondType(Bond::SINGLE);
      } else if (bond->getBondType() == Bond::TRIPLE) {
        bond->setBondType(Bond::DOUBLE);
      }
      stack.pop();
    }
    atom->setFormalCharge(0);
    atom->calcExplicitValence(false);
    return true;
  } else {
    return false;
  }
}

// S- connected to a N via double bond
bool _Valence7SCleanUp3(RWMol& mol, Atom* atom) {
  if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
      atom->calcExplicitValence(false) != 7) {
    return false;
  }

  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 1,
                           nullptr, stack, _visited);
  if (target) {
    stack.top()->setBondType(Bond::SINGLE);
    target->setFormalCharge(-1);
    atom->setFormalCharge(0);
    atom->calcExplicitValence(false);
    return true;
  } else {
    return false;
  }
}

// S- connected to a N via alternating bond
bool _Valence8SCleanUp1(RWMol& mol, Atom* atom) {
  if (atom->getAtomicNum() != 16 || atom->getFormalCharge() != -1 ||
      atom->calcExplicitValence(false) != 7) {
    return false;
  }

  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 7, 0, Bond::DOUBLE, Bond::DOUBLE, 0, 9,
                           nullptr, stack, _visited);

  if (!target) {
    return false;
  }

  while (!stack.empty()) {
    if (stack.top()->getBondType() == Bond::DOUBLE) {
      stack.top()->setBondType(Bond::SINGLE);
    } else {
      stack.top()->setBondType(Bond::DOUBLE);
    }
    stack.pop();
  }
  target->setFormalCharge(-1);
  target->calcExplicitValence(false);
  target->setNumExplicitHs(0);
  atom->setFormalCharge(0);
  atom->calcExplicitValence(false);
  return true;
}

//    [Cl-](=O)(=O)(=O)(=O)
// to:
//    [Cl+3]([O-])([O-])([O-])[O-]
bool _Valence8ClCleanUp1(RWMol& mol, Atom* atom) {
  if (atom->calcExplicitValence(false) != 8 || atom->getFormalCharge() != -1) {
    return false;
  }
  int aid = atom->getIdx();
  bool neighborsAllO = true;
  RWMol::ADJ_ITER nid1, end1;
  boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
  while (nid1 != end1) {
    if (mol.getAtomWithIdx(*nid1)->getAtomicNum() != 8) {
      neighborsAllO = false;
      break;
    }
    nid1++;
  }
  if (neighborsAllO) {
    atom->setFormalCharge(3);
    boost::tie(nid1, end1) = mol.getAtomNeighbors(atom);
    while (nid1 != end1) {
      Bond* b = mol.getBondBetweenAtoms(aid, *nid1);
      if (b->getBondType() == Bond::DOUBLE) {
        b->setBondType(Bond::SINGLE);
        Atom* otherAtom = mol.getAtomWithIdx(*nid1);
        otherAtom->setFormalCharge(-1);
        otherAtom->calcExplicitValence(false);
      }
      nid1++;
    }
    atom->calcExplicitValence(false);
    return true;
  }
  return false;
}

// [Cl+][O-] to Cl=O
bool _Valence5ClCleanUp1(RWMol& mol, Atom* atom) {
  if (atom->calcExplicitValence(false) != 6 || atom->getFormalCharge() != 1) {
    return false;
  }
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 8, -1, Bond::SINGLE, Bond::SINGLE, 0, 1,
                           nullptr, stack, _visited);
  if (!target) {
    return false;
  }
  stack.top()->setBondType(Bond::DOUBLE);
  atom->setFormalCharge(0);
  target->setFormalCharge(0);
  atom->calcExplicitValence(false);
  return true;
}
//
// Cl#S to ClS
bool _Valence3ClCleanUp1(RWMol& mol, Atom* atom) {
  if (atom->calcExplicitValence(false) != 3 || atom->getFormalCharge() != 0) {
    return false;
  }
  std::stack<Bond*> stack;
  std::set<int> _visited;
  Atom* target =
      findAlternatingBonds(mol, atom, 16, 0, Bond::TRIPLE, Bond::TRIPLE, 0, 1,
                           nullptr, stack, _visited);
  if (!target) {
    return false;
  }
  stack.top()->setBondType(Bond::SINGLE);
  atom->calcExplicitValence(false);
  return true;
}

void cleanUp(RWMol& mol) {
  ROMol::AtomIterator ai;
  bool aromHolder;
  for (ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
    switch ((*ai)->getAtomicNum()) {
      case 7:
        if ((*ai)->calcExplicitValence(false) == 4) {
          if (_Valence4NCleanUp1(mol, *ai)) {
            continue;
          }
          if ((*ai)->getFormalCharge() == -1) {
            if (_Valence4NCleanUp2(mol, *ai)) {
              continue;
            }
          }
          continue;
        }

        if ((*ai)->getFormalCharge()) {
          continue;
        }
        aromHolder = (*ai)->getIsAromatic();
        (*ai)->setIsAromatic(0);

        if ((*ai)->calcExplicitValence(false) == 5) {
          // rings CN1=CCN=CC=1, CN1=NCOCC=1, [N]=C1N=CN=N1, [N]=C1C=CN=N1
          (_Valence5NCleanUp6(mol, *ai)) || (_Valence5NCleanUp7(mol, *ai)) ||
              (_Valence5NCleanUp8(mol, *ai)) ||
              (_Valence5NCleanUp9(mol, *ai)) ||
              (_Valence5NCleanUpA(mol, *ai)) ||
              // try search for valence-5 N connected to a N+
              (_Valence5NCleanUp1(mol, *ai)) ||
              // connected to N- through a tiple then single bond
              (_Valence5NCleanUp2(mol, *ai)) ||
              // directly to a N
              (_Valence5NCleanUp3(mol, *ai)) ||
              // to two Si- via double bonds
              (_Valence5NCleanUp4(mol, *ai)) ||
              // alternating bonds to O
              (_Valence5NCleanUp5(mol, *ai, 8)) ||
              // alternating bonds to S
              (_Valence5NCleanUp5(mol, *ai, 16)) ||
              // alternating bonds to S
              (_Valence5NCleanUp5(mol, *ai, 9)) ||
              // alternating bonds to S
              (_Valence5NCleanUp5(mol, *ai, 17)) ||
              // last resort
              (_Valence5NCleanUpB(mol, *ai));
        }
        if (aromHolder) {
          (*ai)->setIsAromatic(1);
        }
        break;
      case 17:
        if ((*ai)->calcExplicitValence(false) == 8 &&
            _Valence8ClCleanUp1(mol, *ai)) {
          continue;
        }
        if ((*ai)->calcExplicitValence(false) == 5 &&
            _Valence5ClCleanUp1(mol, *ai)) {
          continue;
        }
        if ((*ai)->calcExplicitValence(false) == 3 &&
            _Valence3ClCleanUp1(mol, *ai)) {
          continue;
        }
        break;
      case 16:
        if ((*ai)->calcExplicitValence(false) == 7) {
          if (_Valence7SCleanUp1(mol, *ai)) {
            continue;
          }
          if (_Valence7SCleanUp2(mol, *ai)) {
            continue;
          }
          if (_Valence7SCleanUp3(mol, *ai)) {
            continue;
          }
          _Valence8SCleanUp1(mol, *ai);
        } else if ((*ai)->calcExplicitValence(false) == 8) {
          _Valence8SCleanUp1(mol, *ai);
        }
        break;
      case 35:
        if ((*ai)->calcExplicitValence(false) == 3 &&
            (*ai)->getFormalCharge() == 0) {
          // connected to Se. Example: PubChem 10787526
          if ((*ai)->getDegree() == 1) {
            RWMol::ADJ_ITER nid, end;
            boost::tie(nid, end) = mol.getAtomNeighbors(*ai);
            if (mol.getAtomWithIdx(*nid)->getAtomicNum() == 34) {
              mol.getBondBetweenAtoms((*ai)->getIdx(), *nid)
                  ->setBondType(Bond::SINGLE);
            }
          }
        }
        break;

    }  // end the switch block
  }    // end the for loop that iterates over atoms
}  // end cleanUp
}  // namespace

RWMol* InchiToMol(const std::string& inchi, ExtraInchiReturnValues& rv,
                  bool sanitize, bool removeHs) {
  // input
  char* _inchi = new char[inchi.size() + 1];
  char options[1] = "";
  strcpy(_inchi, inchi.c_str());
  inchi_InputINCHI inchiInput;
  inchiInput.szInChI = _inchi;
  inchiInput.szOptions = options;

  // creating RWMol for return
  RWMol* m = nullptr;
  {
    // output structure
    inchi_OutputStruct inchiOutput;
    // DLL call
    int retcode = GetStructFromINCHI(&inchiInput, &inchiOutput);

    // prepare output
    rv.returnCode = retcode;
    if (inchiOutput.szMessage) {
      rv.messagePtr = std::string(inchiOutput.szMessage);
    }
    if (inchiOutput.szLog) {
      rv.logPtr = std::string(inchiOutput.szLog);
    }

    // for isotopes of H
    typedef std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        ISOTOPES_t;
    ISOTOPES_t isotopes;
    if (retcode == inchi_Ret_OKAY || retcode == inchi_Ret_WARNING) {
      m = new RWMol;
      std::vector<unsigned int> indexToAtomIndexMapping;
      PeriodicTable* periodicTable = PeriodicTable::getTable();
      unsigned int nAtoms = inchiOutput.num_atoms;
      for (unsigned int i = 0; i < nAtoms; i++) {
        inchi_Atom* inchiAtom = &(inchiOutput.atom[i]);
        // use element name to set atomic number
        int atomicNumber = periodicTable->getAtomicNumber(inchiAtom->elname);
        Atom* atom = new Atom(atomicNumber);
        double averageWeight = atom->getMass();
        int refWeight = static_cast<int>(averageWeight + 0.5);
        int isotope = 0;
        if (inchiAtom->isotopic_mass) {
          isotope = inchiAtom->isotopic_mass - ISOTOPIC_SHIFT_FLAG;
        }
        if (isotope) {
          atom->setIsotope(isotope + refWeight);
        }
        // set charge
        atom->setFormalCharge(inchiAtom->charge);
        // set radical
        if (inchiAtom->radical) {
          if (inchiAtom->radical != 3 && inchiAtom->radical != 2) {
            BOOST_LOG(rdWarningLog)
                << "expect radical to be either 2 or 3 while getting "
                << inchiAtom->radical << ". Ignore radical." << std::endl;
          } else {
            atom->setNumRadicalElectrons(inchiAtom->radical - 1);
          }
        }
        // number of hydrogens
        atom->setNumExplicitHs(inchiAtom->num_iso_H[0]);
        if (inchiAtom->num_iso_H[1]) {
          isotopes.push_back(std::make_tuple(1, i, inchiAtom->num_iso_H[1]));
        } else if (inchiAtom->num_iso_H[2]) {
          isotopes.push_back(std::make_tuple(2, i, inchiAtom->num_iso_H[2]));
        } else if (inchiAtom->num_iso_H[3]) {
          isotopes.push_back(std::make_tuple(3, i, inchiAtom->num_iso_H[3]));
        }
        // at this point the molecule has all Hs it should have. Set the
        // noImplicit flag so
        // we don't end up with extras later (this was github #562):
        atom->setNoImplicit(true);
        // add atom to molecule
        unsigned int aid = m->addAtom(atom, false, true);
        indexToAtomIndexMapping.push_back(aid);
#ifdef DEBUG
        BOOST_LOG(rdWarningLog)
            << "adding " << aid << ":" << atom->getAtomicNum() << ":"
            << (int)inchiAtom->num_iso_H[0]
            << " charge: " << (int)inchiAtom->charge << std::endl;
#endif
      }

      // adding bonds
      std::set<std::pair<unsigned int, unsigned int>> bondRegister;
      for (unsigned int i = 0; i < nAtoms; i++) {
        inchi_Atom* inchiAtom = &(inchiOutput.atom[i]);
        unsigned int nBonds = inchiAtom->num_bonds;
        for (unsigned int b = 0; b < nBonds; b++) {
          unsigned int nbr = inchiAtom->neighbor[b];
          // check register to avoid duplication
          if (bondRegister.find(std::make_pair(i, nbr)) != bondRegister.end() ||
              bondRegister.find(std::make_pair(nbr, i)) != bondRegister.end()) {
            continue;
          }
          bondRegister.insert(std::make_pair(i, nbr));
          Bond* bond = nullptr;
          // bond type
          if ((unsigned int)inchiAtom->bond_type[b] <= INCHI_BOND_TYPE_TRIPLE) {
            bond = new Bond((Bond::BondType)inchiAtom->bond_type[b]);
          } else if ((unsigned int)inchiAtom->bond_type[b] ==
                     INCHI_BOND_TYPE_ALTERN) {
            BOOST_LOG(rdWarningLog)
                << "receive ALTERN bond type which should be avoided. "
                << "This is treated as aromatic." << std::endl;
            bond = new Bond(Bond::AROMATIC);
            bond->setIsAromatic(true);
          } else {
            BOOST_LOG(rdErrorLog) << "illegal bond type ("
                                  << (unsigned int)inchiAtom->bond_type[b]
                                  << ") in InChI" << std::endl;
            delete m;
            return nullptr;
          }
          // bond ends
          bond->setBeginAtomIdx(indexToAtomIndexMapping[i]);
          bond->setEndAtomIdx(indexToAtomIndexMapping[nbr]);
          // bond stereo
          switch (inchiAtom->bond_stereo[b]) {
            case INCHI_BOND_STEREO_NONE:
              break;
            case INCHI_BOND_STEREO_SINGLE_1UP:
            case INCHI_BOND_STEREO_SINGLE_2DOWN:
              bond->setBondDir(Bond::BEGINWEDGE);
              break;
            case INCHI_BOND_STEREO_SINGLE_1DOWN:
            case INCHI_BOND_STEREO_SINGLE_2UP:
              bond->setBondDir(Bond::BEGINDASH);
              break;
            case INCHI_BOND_STEREO_SINGLE_1EITHER:
              bond->setBondDir(Bond::UNKNOWN);
              break;
            case INCHI_BOND_STEREO_DOUBLE_EITHER:
              bond->setBondDir(Bond::EITHERDOUBLE);
              break;
          }
          // add bond
          m->addBond(bond, true);
#ifdef DEBUG
          BOOST_LOG(rdWarningLog)
              << "adding " << (int)bond->getBeginAtomIdx() << "("
              << m->getAtomWithIdx(bond->getBeginAtomIdx())->getAtomicNum()
              << ")"
              << "-" << (int)bond->getEndAtomIdx() << "("
              << m->getAtomWithIdx(bond->getEndAtomIdx())->getAtomicNum() << ")"
              << "[" << (int)bond->getBondType() << "]" << std::endl;
#endif
        }
      }

      // adding isotopes at the end
      for (auto& ii : isotopes) {
        auto [isotope, aid, repeat] = ii;
        aid = indexToAtomIndexMapping[aid];
        for (unsigned int i = 0; i < repeat; i++) {
          // create atom
          Atom* atom = new Atom;
          atom->setAtomicNum(1);
          // set mass
          atom->setIsotope(isotope);
          int j = m->addAtom(atom, false, true);
          // add bond
          Bond* bond = new Bond(Bond::SINGLE);
          bond->setEndAtomIdx(aid);
          bond->setBeginAtomIdx(j);
          m->addBond(bond, true);
        }
      }

      // basic topological structure is ready. calculate valence
      m->updatePropertyCache(false);

      // 0Dstereo
      INT_PAIR_VECT eBondPairs;
      INT_PAIR_VECT zBondPairs;
      unsigned int numStereo0D = inchiOutput.num_stereo0D;
      if (numStereo0D) {
        // calculate CIPCode as they might be used
        UINT_VECT ranks;
        Chirality::assignAtomCIPRanks(*m, ranks);
        for (unsigned int i = 0; i < numStereo0D; i++) {
          inchi_Stereo0D* stereo0DPtr = inchiOutput.stereo0D + i;
          if (stereo0DPtr->parity == INCHI_PARITY_NONE ||
              stereo0DPtr->parity == INCHI_PARITY_UNDEFINED) {
            continue;
          }
          switch (stereo0DPtr->type) {
            case INCHI_StereoType_None:
              break;
            case INCHI_StereoType_DoubleBond: {
              // find the bond
              unsigned left = indexToAtomIndexMapping[stereo0DPtr->neighbor[1]];
              unsigned right =
                  indexToAtomIndexMapping[stereo0DPtr->neighbor[2]];
              int originalLeftNbr =
                  indexToAtomIndexMapping[stereo0DPtr->neighbor[0]];
              int originalRightNbr =
                  indexToAtomIndexMapping[stereo0DPtr->neighbor[3]];

              Bond* bond = m->getBondBetweenAtoms(left, right);
              if (!bond) {
                // Likely to be allene stereochemistry, which we don't handle.
                BOOST_LOG(rdWarningLog)
                    << "Extended double-bond stereochemistry (e.g. C=C=C=C) "
                       "ignored"
                    << std::endl;
                continue;
              }
              // also find neighboring atoms. Note we cannot use what InChI
              // returned in stereo0DPtr->neighbor as there can be hydrogen in
              // it, which is later removed and is therefore not reliable. Plus,
              // InChI seems to use lower CIPRank-neighbors rather than
              // higher-CIPRank ones (hence the use of hydrogen neighbor).
              // However, if the neighbors we selected differ from what are in
              // stereo0DPtr->neighbor, we might also need to switch E and Z

              auto findNbrAtoms = [&m, &ranks](unsigned ref) {
                int nbr = -1;
                int extraNbr = -1;
                int cip = -1;
                int _cip = -1;
                for (auto bond : m->atomBonds(m->getAtomWithIdx(ref))) {
                  if (bond->getBondType() != Bond::SINGLE &&
                      bond->getBondType() != Bond::AROMATIC) {
                    continue;
                  }
                  auto atom = bond->getOtherAtomIdx(ref);
                  if ((_cip = ranks[atom]) > cip) {
                    if (nbr >= 0) {
                      extraNbr = nbr;
                    }
                    nbr = atom;
                    cip = _cip;
                  } else {
                    extraNbr = atom;
                  }
                }
                return std::make_pair(nbr, extraNbr);
              };
              auto [leftNbr, extraLeftNbr] = findNbrAtoms(left);
              auto [rightNbr, extraRightNbr] = findNbrAtoms(right);

              if (leftNbr < 0 || rightNbr < 0) {
                BOOST_LOG(rdWarningLog)
                    << "Ignoring stereochemistry on double-bond without appropriate neighbors"
                    << std::endl;
                continue;
              }

              bool switchEZ = false;
              if ((originalLeftNbr == leftNbr &&
                   originalRightNbr != rightNbr) ||
                  (originalLeftNbr != leftNbr &&
                   originalRightNbr == rightNbr)) {
                switchEZ = true;
              }

              char parity = stereo0DPtr->parity;
              if (parity == INCHI_PARITY_ODD && switchEZ) {
                parity = INCHI_PARITY_EVEN;
              } else if (parity == INCHI_PARITY_EVEN && switchEZ) {
                parity = INCHI_PARITY_ODD;
              }

              auto findBondPairs = [&m, &zBondPairs, &eBondPairs](
                                       unsigned ref, int nbr, int extraNbr) {
                auto bond = m->getBondBetweenAtoms(ref, nbr);
                if (extraNbr >= 0) {
                  // modifier to track whether bond is reversed
                  int modifier = -1;
                  if (bond->getBeginAtomIdx() != ref) {
                    modifier *= -1;
                  }
                  auto extraBond = m->getBondBetweenAtoms(ref, extraNbr);
                  if (extraBond->getBeginAtomIdx() != ref) {
                    modifier *= -1;
                  }
                  if (modifier == 1) {
                    zBondPairs.push_back(
                        std::make_pair(bond->getIdx(), extraBond->getIdx()));
                  } else {
                    eBondPairs.push_back(
                        std::make_pair(bond->getIdx(), extraBond->getIdx()));
                  }
                }
                return bond;
              };

              auto leftBond = findBondPairs(left, leftNbr, extraLeftNbr);
              auto rightBond = findBondPairs(right, rightNbr, extraRightNbr);

              int modifier = -1;  // modifier to track whether bond is reversed
              if (leftBond->getBeginAtomIdx() != left) {
                modifier *= -1;
              }
              if (rightBond->getBeginAtomIdx() != right) {
                modifier *= -1;
              }

              if (parity == INCHI_PARITY_ODD) {
                bond->setStereo(Bond::STEREOZ);
                if (modifier == 1) {
                  eBondPairs.push_back(
                      std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
                } else {
                  zBondPairs.push_back(
                      std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
                }
              } else if (parity == INCHI_PARITY_EVEN) {
                bond->setStereo(Bond::STEREOE);
                if (modifier == 1) {
                  zBondPairs.push_back(
                      std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
                } else {
                  eBondPairs.push_back(
                      std::make_pair(leftBond->getIdx(), rightBond->getIdx()));
                }
              } else if (parity == INCHI_PARITY_NONE) {
                bond->setStereo(Bond::STEREONONE);
              } else {
                bond->setStereo(Bond::STEREOANY);
              }
              // set the stereo atoms for the double bond
              bond->getStereoAtoms().push_back(leftNbr);
              bond->getStereoAtoms().push_back(rightNbr);
              break;
            }
            case INCHI_StereoType_Tetrahedral: {
              unsigned int c =
                  indexToAtomIndexMapping[stereo0DPtr->central_atom];
              Atom* atom = m->getAtomWithIdx(c);
              // find number of swaps for the members
              int nSwaps = 0;
              unsigned int nid = 0;
              if (stereo0DPtr->neighbor[0] == stereo0DPtr->central_atom) {
                // 3-neighbor case
                nid = 1;
                if (atom->getDegree() == 3) {
                  // this happens with chiral three-coordinate S
                  nSwaps = 1;
                }
              }
              // if (atom->getTotalNumHs(true) == 1)
              //  nSwaps = 1;
              // std::cerr<<"build atom: "<<c<<" "<<atom->getTotalNumHs(true);
              std::list<int> neighbors;
              for (; nid < 4; nid++) {
                unsigned end =
                    indexToAtomIndexMapping[stereo0DPtr->neighbor[nid]];
                Bond* bond = m->getBondBetweenAtoms(c, end);
                neighbors.push_back(bond->getIdx());
                // std::cerr<<" "<<end<<"("<<bond->getIdx()<<")";
              }
              nSwaps += atom->getPerturbationOrder(neighbors);
              // std::cerr<<" swaps: "<<nSwaps<<" parity: "<<
              //  (stereo0DPtr->parity==INCHI_PARITY_EVEN?"even":"odd")<<std::endl;
              if (stereo0DPtr->parity == INCHI_PARITY_ODD) {
                atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
              } else {
                atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
              }
              if (nSwaps % 2) {
                atom->invertChirality();
              }
              break;
            }
            case INCHI_StereoType_Allene:
              BOOST_LOG(rdWarningLog) << "Allene-style stereochemistry is not "
                                         "supported yet and will be ignored."
                                      << std::endl;
              break;
            default:
              BOOST_LOG(rdWarningLog)
                  << "Unrecognized stereo0D type (" << (int)stereo0DPtr->type
                  << ") is ignored!" << std::endl;
          }  // end switch stereotype
        }    // end for loop over all stereo0D entries
        // set the bond directions
        if (!assignBondDirs(*m, zBondPairs, eBondPairs)) {
          BOOST_LOG(rdWarningLog)
              << "Cannot assign bond directions!" << std::endl;
          ;
        }
      }  // end if (if stereo0D presents)
    }    // end if (if return code is success)

    // clean up
    delete[] _inchi;
    FreeStructFromINCHI(&inchiOutput);
  }

  // clean up the molecule to be acceptable to RDKit
  if (m) {
    cleanUp(*m);
    try {
      if (sanitize) {
        if (removeHs) {
          MolOps::removeHs(*m, false, false);
        } else {
          MolOps::sanitizeMol(*m);
        }
      }
    } catch (const MolSanitizeException&) {
      delete m;
      throw;
    }
    // call assignStereochemistry just to be safe; otherwise, MolToSmiles may
    // overwrite E/Z and/or bond direction on double bonds.
    MolOps::assignStereochemistry(*m, true, true);
  }

  return m;
}

void fixOptionSymbol(const char* in, char* out) {
  unsigned int i;
  for (i = 0; i < strlen(in); i++) {
#ifdef _WIN32
    if (in[i] == '-') {
      out[i] = '/';

#else
    if (in[i] == '/') {
      out[i] = '-';

#endif
    } else {
      out[i] = in[i];
    }
  }
  out[i] = '\0';
}

/*! "reverse" clean up: prepare a molecule to be used with InChI sdk */
void rCleanUp(RWMol& mol) {
  RWMol* q = SmilesToMol("[O-][Cl+3]([O-])([O-])O");
  std::vector<MatchVectType> fgpMatches;
  SubstructMatch(mol, *q, fgpMatches);
  delete q;
  // replace all matches
  for (auto match : fgpMatches) {
    // collect matching atoms
    int map[5];
    for (MatchVectType::const_iterator mi = match.begin(); mi != match.end();
         mi++) {
      map[mi->first] = mi->second;
    }
    // check charges
    if (mol.getAtomWithIdx(map[1])->getFormalCharge() != 3) {
      return;
    }
    int unchargedFound = -1;
    for (int i = 0; i < 5; i++) {
      if (i == 1) {
        continue;
      }
      Atom* o = mol.getAtomWithIdx(map[i]);
      if (o->getFormalCharge() == 0) {
        if (unchargedFound != -1) {
          return;  // too many uncharged oxygen
        } else {
          unchargedFound = i;
        }
      }
    }

    // flip bonds and remove charges
    for (int i = 0; i < 5; i++) {
      if (i == 1) {
        continue;
      }
      if (i == unchargedFound) {
        continue;
      }
      if (unchargedFound == -1 && i == 0) {
        mol.getBondBetweenAtoms(map[1], map[i])->setBondType(Bond::SINGLE);
        mol.getAtomWithIdx(map[i])->setFormalCharge(-1);
        continue;
      }
      mol.getBondBetweenAtoms(map[1], map[i])->setBondType(Bond::DOUBLE);
      mol.getAtomWithIdx(map[i])->setFormalCharge(0);
    }
    mol.getAtomWithIdx(map[1])->setFormalCharge(0);
  }
  return;
}

std::string MolToInchi(const ROMol& mol, ExtraInchiReturnValues& rv,
                       const char* options) {
  std::unique_ptr<RWMol> m{new RWMol(mol)};

  // assign stereochem:
  if (mol.needsUpdatePropertyCache()) {
    m->updatePropertyCache(false);
  }
  // kekulize
  MolOps::Kekulize(*m, false);

  // "reverse" cleanup: undo some clean up done by RDKit
  rCleanUp(*m);

  unsigned int nAtoms = m->getNumAtoms();
  unsigned int nBonds = m->getNumBonds();

  // Make array of inchi_atom (storage space)
  auto* inchiAtoms = new inchi_Atom[nAtoms];
  // and a vector for stereo0D
  std::vector<inchi_Stereo0D> stereo0DEntries;

  PeriodicTable* periodicTable = PeriodicTable::getTable();
  // Fill inchi_Atom's by atoms in RWMol
  for (unsigned int i = 0; i < nAtoms; i++) {
    Atom* atom = m->getAtomWithIdx(i);
    inchiAtoms[i].num_bonds = 0;

    // coordinates
    if (!m->getNumConformers()) {
      inchiAtoms[i].x = 0;
      inchiAtoms[i].y = 0;
      inchiAtoms[i].z = 0;
    } else {
      auto conformerIter = m->beginConformers();
      RDGeom::Point3D coord = (*conformerIter)->getAtomPos(i);
      inchiAtoms[i].x = coord[0];
      inchiAtoms[i].y = coord[1];
      inchiAtoms[i].z = coord[2];
    }

    // element name
    unsigned int atomicNumber = atom->getAtomicNum();
    std::string elementName = periodicTable->getElementSymbol(atomicNumber);
    strcpy(inchiAtoms[i].elname, elementName.c_str());

    // isotopes
    int isotope = atom->getIsotope();
    if (isotope) {
      inchiAtoms[i].isotopic_mass =
          ISOTOPIC_SHIFT_FLAG + isotope -
          static_cast<int>(periodicTable->getAtomicWeight(atomicNumber) + 0.5);
    } else {
      // check explicit iso property. If this is set, we have a 0 offset
      // Example: CHEMBL220875
      // if (atom->getIsotope()){
      //  inchiAtoms[i].isotopic_mass = ISOTOPIC_SHIFT_FLAG + 0;
      //} else {
      inchiAtoms[i].isotopic_mass = 0;
      //}
    }

    // charge
    inchiAtoms[i].charge = atom->getFormalCharge();

    // number of iso H
    int nHs = -1;
    switch (atom->getAtomicNum()) {
      case 6:
      case 7:
      case 8:
      case 9:
      case 17:
      case 35:
      case 53:
        nHs = -1;
        break;
      default:
        nHs = atom->getTotalNumHs();
    }
    inchiAtoms[i].num_iso_H[0] = nHs;
    inchiAtoms[i].num_iso_H[1] = 0;
    inchiAtoms[i].num_iso_H[2] = 0;
    inchiAtoms[i].num_iso_H[3] = 0;

    // radical
    inchiAtoms[i].radical = 0;
    if (atom->getNumRadicalElectrons()) {
      // the direct specification of radicals in InChI is tricky since they use
      // the MDL representation (singlet, double, triplet) and we just have the
      // number of unpaired electrons. Instead we set the number of implicit Hs
      // here, that together with the atom identity and charge should be
      // sufficient
      inchiAtoms[i].num_iso_H[0] = atom->getTotalNumHs();
    } else {
    }

    // convert tetrahedral chirality info to Stereo0D
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED ||
        atom->hasProp("molParity")) {
      // we ignore the molParity if the number of neighbors are below 3
      atom->calcImplicitValence();
      if (atom->getNumImplicitHs() + atom->getDegree() < 3) {
        continue;
      }
      inchi_Stereo0D stereo0D;
      stereo0D.central_atom = i;
      stereo0D.type = INCHI_StereoType_Tetrahedral;
      ROMol::ADJ_ITER nbrIter, endNbrIter;
      boost::tie(nbrIter, endNbrIter) = m->getAtomNeighbors(atom);
      std::vector<std::pair<unsigned int, unsigned int>> neighbors;
      while (nbrIter != endNbrIter) {
        int cip = 0;
        // if (m->getAtomWithIdx(*nbrIter)->hasProp("_CIPRank"))
        //   m->getAtomWithIdx(*nbrIter)->getProp("_CIPRank", cip);
        neighbors.emplace_back(cip, *nbrIter);
        ++nbrIter;
      }
      // std::sort(neighbors.begin(), neighbors.end());
      unsigned char nid = 0;
      // std::cerr<<" at: "<<atom->getIdx();
      for (const auto& p : neighbors) {
        stereo0D.neighbor[nid++] = p.second;
        // std::cerr<<" "<<p.second;
      }
      if (nid == 3) {
        // std::cerr<<" nid==3, reorder";
        // std::cerr<<" "<<i;
        for (; nid > 0; nid--) {
          stereo0D.neighbor[nid] = stereo0D.neighbor[nid - 1];
          // std::cerr<<" "<<stereo0D.neighbor[nid];
        }
        stereo0D.neighbor[0] = i;
      }
      // std::cerr<<std::endl;
      Atom::ChiralType chiralTag;
      if ((chiralTag = atom->getChiralTag()) != Atom::CHI_UNSPECIFIED) {
        bool pushIt = false;
        if (atom->getDegree() == 4) {
          if (chiralTag == Atom::CHI_TETRAHEDRAL_CW) {
            stereo0D.parity = INCHI_PARITY_EVEN;
            pushIt = true;
          } else {
            stereo0D.parity = INCHI_PARITY_ODD;
            pushIt = true;
          }
        } else {
          // std::cerr<<"tag: "<<chiralTag<<std::endl;
          if (chiralTag == Atom::CHI_TETRAHEDRAL_CCW) {
            stereo0D.parity = INCHI_PARITY_EVEN;
            pushIt = true;
          } else if (chiralTag == Atom::CHI_TETRAHEDRAL_CW) {
            stereo0D.parity = INCHI_PARITY_ODD;
            pushIt = true;
          } else {
            BOOST_LOG(rdWarningLog)
                << "unrecognized chirality tag (" << chiralTag << ") on atom "
                << i << " is ignored." << std::endl;
          }
        }
        if (pushIt) {
          // this was github #296
          // with molecules like C[S@@](=O)C(C)(C)C the stereochem of the sulfur
          // from
          // the inchi comes back reversed if we don't have wedged bonds. There
          // must
          // be something with the way S stereochem is being handled that I'm
          // not
          // getting.
          // There's something of an explanation at around line 258 of
          // inchi_api.h
          // but that didn't help that much.
          // For want of a better idea, detect this pattern
          // and flip the stereochem:
          // if(atom->getAtomicNum()==16 &&
          //    atom->getDegree()==3 && atom->getExplicitValence()==4){
          //   if(stereo0D.parity==INCHI_PARITY_EVEN){
          //     stereo0D.parity=INCHI_PARITY_ODD;
          //   } else if(stereo0D.parity==INCHI_PARITY_ODD){
          //     stereo0D.parity=INCHI_PARITY_EVEN;
          //   }
          // }
          stereo0DEntries.push_back(stereo0D);
        }

      } else {
        // std::string molParity;
        // atom->getProp("molParity", molParity);
        // if (molParity == "2") {
        //  stereo0D.parity = INCHI_PARITY_EVEN;
        //  stereo0DEntries.push_back(stereo0D);
        //} else if (molParity == "1") {
        //  stereo0D.parity = INCHI_PARITY_ODD;
        //  stereo0DEntries.push_back(stereo0D);
        //} else if (molParity == "0") {
        //  stereo0D.parity = INCHI_PARITY_NONE;
        //  stereo0DEntries.push_back(stereo0D);
        //} else if (molParity == "3") {
        //  stereo0D.parity = INCHI_PARITY_UNKNOWN;
        //  stereo0DEntries.push_back(stereo0D);
        //} else {
        //  BOOST_LOG(rdWarningLog) << "unrecognized parity on atom "
        //    << molParity << " is ignored." << std::endl;
        //}
      }
    }
  }

  // read bond info
  for (unsigned int i = 0; i < nBonds; i++) {
    Bond* bond = m->getBondWithIdx(i);
    unsigned int atomIndex1 = bond->getBeginAtomIdx();
    unsigned int atomIndex2 = bond->getEndAtomIdx();
    int bondDirectionModifier = 1;
    // update only for the atom having smaller index
    if (atomIndex1 > atomIndex2) {
      std::swap(atomIndex1, atomIndex2);
      bondDirectionModifier = -1;
    }

    // neighbor
    unsigned int idx = inchiAtoms[atomIndex1].num_bonds;
    inchiAtoms[atomIndex1].neighbor[idx] = atomIndex2;

    // bond type
    Bond::BondType bondType = bond->getBondType();
    if (bondType > Bond::TRIPLE) {
      BOOST_LOG(rdWarningLog) << "bond type above 3 (" << bondType
                              << ") is treated as unspecified!" << std::endl;
      bondType = Bond::UNSPECIFIED;
    }
    inchiAtoms[atomIndex1].bond_type[idx] = bondType;

    // stereo
    Bond::BondDir bondDirection = bond->getBondDir();
    switch (bondDirection) {
      case Bond::BEGINWEDGE:
        inchiAtoms[atomIndex1].bond_stereo[idx] =
            bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1UP;
        break;
      case Bond::BEGINDASH:
        inchiAtoms[atomIndex1].bond_stereo[idx] =
            bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1DOWN;
        break;
      case Bond::EITHERDOUBLE:
        inchiAtoms[atomIndex1].bond_stereo[idx] =
            INCHI_BOND_STEREO_DOUBLE_EITHER;
        break;
      case Bond::UNKNOWN:
        inchiAtoms[atomIndex1].bond_stereo[idx] =
            bondDirectionModifier * INCHI_BOND_STEREO_SINGLE_1EITHER;
        break;
      case Bond::NONE:
      default:
        inchiAtoms[atomIndex1].bond_stereo[idx] = INCHI_BOND_STEREO_NONE;
    }

    // double bond stereochemistry
    // single bond in the big ring will get E/Z assigned as well. Though rdkit
    // will eventually remove it, I added it any way
    if (  // bondType == Bond::DOUBLE and
        bond->getStereo() > Bond::STEREOANY &&
        bond->getStereoAtoms().size() >= 2) {
      inchi_Stereo0D stereo0D;
      if (bond->getStereo() == Bond::STEREOZ ||
          bond->getStereo() == Bond::STEREOCIS) {
        stereo0D.parity = INCHI_PARITY_ODD;
      } else {
        stereo0D.parity = INCHI_PARITY_EVEN;
      }
      stereo0D.neighbor[0] = bond->getStereoAtoms()[0];
      stereo0D.neighbor[3] = bond->getStereoAtoms()[1];
      stereo0D.neighbor[1] = atomIndex1;
      stereo0D.neighbor[2] = atomIndex2;
      if (!m->getBondBetweenAtoms(stereo0D.neighbor[0], stereo0D.neighbor[1])) {
        std::swap(stereo0D.neighbor[0], stereo0D.neighbor[3]);
      }
      stereo0D.central_atom = NO_ATOM;
      stereo0D.type = INCHI_StereoType_DoubleBond;
      stereo0DEntries.push_back(stereo0D);
    } else if (bond->getStereo() == Bond::STEREOANY) {
      // have to treat STEREOANY separately because RDKit will clear out
      // StereoAtoms information.
      // Here we just change the coordinates of the two end atoms - to bring
      // them really close - so that InChI will not try to infer stereobond
      // info from coordinates.
      inchiAtoms[atomIndex1].x = inchiAtoms[atomIndex2].x;
      inchiAtoms[atomIndex1].y = inchiAtoms[atomIndex2].y;
      inchiAtoms[atomIndex1].z = inchiAtoms[atomIndex2].z;
    }

    // number of bonds
    inchiAtoms[atomIndex1].num_bonds++;
  }

  // create stereo0D
  inchi_Stereo0D* stereo0Ds;
  if (stereo0DEntries.size()) {
    stereo0Ds = new inchi_Stereo0D[stereo0DEntries.size()];
    for (unsigned int i = 0; i < stereo0DEntries.size(); i++) {
      stereo0Ds[i] = stereo0DEntries[i];
    }
  } else {
    stereo0Ds = nullptr;
  }

  // create input
  inchi_Input input;
  input.atom = inchiAtoms;
  input.stereo0D = stereo0Ds;
  if (options) {
    char* _options = new char[strlen(options) + 1];
    fixOptionSymbol(options, _options);
    input.szOptions = _options;
  } else {
    input.szOptions = nullptr;
  }
  input.num_atoms = nAtoms;
  input.num_stereo0D = stereo0DEntries.size();

  // create output
  inchi_Output output;

  // call DLL
  std::string inchi;
  {
    int retcode = GetINCHI(&input, &output);

    // generate output
    rv.returnCode = retcode;
    if (output.szInChI) {
      inchi = std::string(output.szInChI);
    }
    if (output.szMessage) {
      rv.messagePtr = std::string(output.szMessage);
    }
    if (output.szLog) {
      rv.logPtr = std::string(output.szLog);
    }
    if (output.szAuxInfo) {
      rv.auxInfoPtr = std::string(output.szAuxInfo);
    }

    // clean up
    FreeINCHI(&output);
  }
  if (input.szOptions) {
    delete[] input.szOptions;
  }

  delete[] inchiAtoms;
  if (stereo0Ds) {
    delete[] stereo0Ds;
  }

  return inchi;
}

std::string MolBlockToInchi(const std::string& molBlock,
                            ExtraInchiReturnValues& rv, const char* options) {
  // create output
  inchi_Output output;
  memset((void*)&output, 0, sizeof(output));
  // call DLL
  std::string inchi;
  {
    char* _options = nullptr;
    if (options) {
      _options = new char[strlen(options) + 1];
      fixOptionSymbol(options, _options);
      options = _options;
    }
    int retcode =
        MakeINCHIFromMolfileText(molBlock.c_str(), (char*)options, &output);

    // generate output
    rv.returnCode = retcode;
    if (output.szInChI) {
      inchi = std::string(output.szInChI);
    }
    if (output.szMessage) {
      rv.messagePtr = std::string(output.szMessage);
    }
    if (output.szLog) {
      rv.logPtr = std::string(output.szLog);
    }
    if (output.szAuxInfo) {
      rv.auxInfoPtr = std::string(output.szAuxInfo);
    }

    // clean up
    FreeINCHI(&output);
    delete[] _options;
  }
  return inchi;
}

std::string InchiToInchiKey(const std::string& inchi) {
  char inchiKey[29];
  char xtra1[65], xtra2[65];
  int ret = 0;
  { ret = GetINCHIKeyFromINCHI(inchi.c_str(), 0, 0, inchiKey, xtra1, xtra2); }
  std::string error;
  switch (ret) {
    case INCHIKEY_OK:
      return std::string(inchiKey);
    case INCHIKEY_UNKNOWN_ERROR:
      error = "Unknown error";
      break;
    case INCHIKEY_EMPTY_INPUT:
      error = "Empty input";
      break;
    case INCHIKEY_INVALID_INCHI_PREFIX:
      error = "Invalid InChI prefix";
      break;
    case INCHIKEY_NOT_ENOUGH_MEMORY:
      error = "Not enough memory";
      break;
    case INCHIKEY_INVALID_INCHI:
      error = "Invalid input InChI string";
      break;
    case INCHIKEY_INVALID_STD_INCHI:
      error = "Invalid standard InChI string";
      break;
  }
  BOOST_LOG(rdErrorLog) << error << " in generating InChI Key" << std::endl;
  return std::string();
}
}  // namespace RDKit
