//
//  Copyright (C) 2018-2020 Susan H. Leung and Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Tautomer.h"
#include "Fragment.h"
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <limits>

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <utility>

using namespace RDKit;

namespace RDKit {

namespace MolStandardize {

namespace detail {
std::vector<std::pair<unsigned int, unsigned int>> pairwise(
    const std::vector<int> vect) {
  std::vector<std::pair<unsigned int, unsigned int>> pvect;
  for (size_t i = 0; i < vect.size() - 1; ++i) {
    std::pair<unsigned int, unsigned int> p =
        std::pair<unsigned int, unsigned int>(vect[i], vect[i + 1]);
    pvect.push_back(p);
  }
  return pvect;
}

}  // namespace detail

namespace TautomerScoringFunctions {
int scoreRings(const ROMol &mol) {
  int score = 0;
  auto ringInfo = mol.getRingInfo();
  std::unique_ptr<ROMol> cp;
  if (!ringInfo->isInitialized()) {
    cp.reset(new ROMol(mol));
    MolOps::symmetrizeSSSR(*cp);
    ringInfo = cp->getRingInfo();
  }
  boost::dynamic_bitset<> isArom(mol.getNumBonds());
  boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
  for (const auto &bnd : mol.bonds()) {
    if (bnd->getIsAromatic()) {
      isArom.set(bnd->getIdx());
      if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
          bnd->getEndAtom()->getAtomicNum() == 6) {
        bothCarbon.set(bnd->getIdx());
      }
    }
  }
  for (const auto &bring : ringInfo->bondRings()) {
    bool allC = true;
    bool allAromatic = true;
    for (const auto bidx : bring) {
      if (!isArom[bidx]) {
        allAromatic = false;
        break;
      }
      if (!bothCarbon[bidx]) {
        allC = false;
      }
    }
    if (allAromatic) {
      score += 100;
      if (allC) {
        score += 150;
      }
    }
  }
  return score;
};

struct smarts_mol_holder {
  std::string d_smarts;
  std::shared_ptr<ROMol> dp_mol;
  smarts_mol_holder(const std::string &smarts) : d_smarts(smarts) {
    dp_mol.reset(SmartsToMol(smarts));
  }
};

typedef boost::flyweight<
    boost::flyweights::key_value<std::string, smarts_mol_holder>,
    boost::flyweights::no_tracking>
    smarts_mol_flyweight;

struct SubstructTerm {
  std::string name;
  std::string smarts;
  int score;
  std::shared_ptr<ROMol> matcher;
  SubstructTerm(std::string aname, std::string asmarts, int ascore)
      : name(std::move(aname)), smarts(std::move(asmarts)), score(ascore) {
    matcher = smarts_mol_flyweight(smarts).get().dp_mol;
  };
};

int scoreSubstructs(const ROMol &mol) {
  // a note on efficiency here: we'll construct the SubstructTerm objects here
  // repeatedly, but the SMARTS parsing for each entry will only be done once
  // since we're using the boost::flyweights above to cache them
  const std::vector<SubstructTerm> substructureTerms{
      {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]",
       25},
      {"oxim", "[#6]=[N][OH]", 4},
      {"C=O", "[#6]=,:[#8]", 2},
      {"N=O", "[#7]=,:[#8]", 2},
      {"P=O", "[#15]=,:[#8]", 2},
      {"C=hetero", "[#6]=[!#1;!#6]", 1},
      {"methyl", "[CX4H3]", 1},
      {"guanidine terminal=N", "[#7][#6](=[NR0])[#7H0]", 1},
      {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2},
      {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4}};
  int score = 0;
  for (const auto &term : substructureTerms) {
    if (!term.matcher) {
      BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
                            << " is invalid, ignoring it." << std::endl;
      continue;
    }
    SubstructMatchParameters params;
    const auto matches = SubstructMatch(mol, *term.matcher, params);
    // if (!matches.empty()) {
    //   std::cerr << " " << matches.size() << " matches to " << term.name
    //             << std::endl;
    // }
    score += matches.size() * term.score;
  }
  return score;
};
int scoreHeteroHs(const ROMol &mol) {
  int score = 0;
  for (const auto &at : mol.atoms()) {
    int anum = at->getAtomicNum();
    if (anum == 15 || anum == 16 || anum == 34 || anum == 52) {
      score -= at->getTotalNumHs();
    }
  }
  return score;

  return 1;
};
}  // namespace TautomerScoringFunctions

unsigned int MAX_TAUTOMERS = 1000;

ROMol *TautomerEnumerator::pickCanonical(
    const std::vector<ROMOL_SPTR> &tautomers,
    boost::function<int(const ROMol &mol)> scoreFunc) const {
  PRECONDITION(scoreFunc, "no scoring function");
  if (tautomers.size() == 1) {
    return new ROMol(*tautomers[0]);
  }
  // Calculate score for each tautomer
  int bestScore = std::numeric_limits<int>::min();
  std::string bestSmiles = "";
  ROMOL_SPTR bestMol;
  for (const auto t : tautomers) {
    auto score = scoreFunc(*t);
    if (score > bestScore) {
      bestScore = score;
      bestSmiles = MolToSmiles(*t);
      bestMol = t;
    } else if (score == bestScore) {
      auto smiles = MolToSmiles(*t);
      if (smiles < bestSmiles) {
        bestSmiles = smiles;
        bestMol = t;
      }
    }
  }
  return new ROMol(*bestMol);
}

std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(const ROMol &mol) const {
  // std::cout << "**********************************" << std::endl;

  PRECONDITION(dp_catalog, "no catalog!");
  const TautomerCatalogParams *tautparams = dp_catalog->getCatalogParams();

  PRECONDITION(tautparams, "");
  const std::vector<TautomerTransform> &transforms =
      tautparams->getTransforms();

  // Enumerate all possible tautomers and return them as a list.
  std::string smi = MolToSmiles(mol, true);
  boost::shared_ptr<ROMol> taut(new ROMol(mol));
  std::map<std::string, boost::shared_ptr<ROMol>> tautomers = {{smi, taut}};
  // Create a kekulized form of the molecule to match the SMARTS against
  boost::shared_ptr<RWMol> kekulized(new RWMol(mol));
  MolOps::Kekulize(*kekulized, false);
  std::map<std::string, boost::shared_ptr<ROMol>> kekulized_mols = {
      {smi, kekulized}};
  std::vector<std::string> done;
  bool broken = false;

  while (tautomers.size() < MAX_TAUTOMERS) {
    // std::map automatically sorts tautomers into alphabetical order (SMILES)
    for (const auto &tautomer : tautomers) {
      // std::cout << "Done : " << std::endl;
      // for (const auto d : done) {
      //   std::cout << d << std::endl;
      // }
      // std::cout << "Looking at tautomer: " << tautomer.first << std::endl;
      std::string tsmiles;
      if (std::find(done.begin(), done.end(), tautomer.first) != done.end()) {
        continue;
      } else {
        // done does not contain tautomer
        for (const auto &transform : transforms) {
          // find kekulized_mol in kekulized_mols with same smiles as taut
          auto kmol = kekulized_mols.find(tautomer.first);
          //					if (search !=
          // kekulized_mols.end() 					for
          // (const auto &mol : kekulized_mols) { if (mol.first ==
          // tautomer.first) { std::cout << mol.first << std::endl;
          //							}
          //					std::cout <<
          // MolToSmiles(*transform.Mol)
          //<< std::endl;
          std::vector<MatchVectType> matches;
          unsigned int matched =
              SubstructMatch(*(kmol->second), *(transform.Mol), matches);
          std::string name;
          (transform.Mol)->getProp(common_properties::_Name, name);

          if (!matched) {
            continue;
          } else {
            // std::cout << "kmol: " << kmol->first << std::endl;
            // std::cout << MolToSmiles(*(kmol->second)) << std::endl;
            // std::cout << "transform mol: " << MolToSmarts(*(transform.Mol))
            //           << std::endl;

            // std::cout << "Matched: " << name << std::endl;
          }
          for (const auto &match : matches) {
            std::vector<int> idx_matches;
            for (const auto &pair : match) {
              idx_matches.push_back(pair.second);
            }
            // Create a copy of in the input molecule so we can modify it
            // Use kekule form so bonds are explicitly single/double instead of
            // aromatic
            boost::shared_ptr<ROMol> product(new ROMol(*(kmol->second)));
            // Remove a hydrogen from the first matched atom and add one to the
            // last
            Atom *first = product->getAtomWithIdx(idx_matches[0]);
            Atom *last = product->getAtomWithIdx(idx_matches.back());
            first->setNumExplicitHs(
                std::max((unsigned int)0, first->getTotalNumHs() - 1));
            last->setNumExplicitHs(last->getTotalNumHs() + 1);
            // Remove any implicit hydrogens from the first and last atoms
            // now we have set the count explicitly
            first->setNoImplicit(true);
            last->setNoImplicit(true);
            // Adjust bond orders
            unsigned int bi = 0;
            std::vector<std::pair<unsigned int, unsigned int>> pvect =
                detail::pairwise(idx_matches);
            for (const auto &pair : pvect) {
              Bond *bond =
                  product->getBondBetweenAtoms(pair.first, pair.second);
              // check if bonds is specified in tatuomer.in file
              if (!transform.BondTypes.empty()) {
                bond->setBondType(transform.BondTypes[bi]);
                ++bi;
              } else {
                Bond::BondType bondtype = bond->getBondType();
                //								std::cout
                //<< "Bond as double: " << bond->getBondTypeAsDouble() <<
                // std::endl;
                // std::cout
                // << bondtype << std::endl;
                if (bondtype == 1) {
                  bond->setBondType(Bond::DOUBLE);
                  //									std::cout
                  //<< "Set bond to double" << std::endl;
                }
                if (bondtype == 2) {
                  bond->setBondType(Bond::SINGLE);
                  //									std::cout
                  //<< "Set bond to single" << std::endl;
                }
              }
            }
            // TODO adjust charges
            if (!transform.Charges.empty()) {
              unsigned int ci = 0;
              for (const auto idx : idx_matches) {
                Atom *atom = product->getAtomWithIdx(idx);
                atom->setFormalCharge(atom->getFormalCharge() +
                                      transform.Charges[ci]);
                ++ci;
              }
            }

            boost::shared_ptr<RWMol> wproduct(new RWMol(*product));
            // wproduct->updatePropertyCache(false);
            // std::cout << "pre-sanitization: "
            //           << MolToSmiles(*wproduct, true, true) << std::endl;
            MolOps::sanitizeMol(*wproduct);
            //						MolOps::sanitizeMol(*static_cast<RWMol*>(product.get()));
            tsmiles = MolToSmiles(*wproduct, true);
            //						std::string name;
            //						(transform.Mol)->getProp(common_properties::_Name,
            // name);
            // std::cout << "Applied rule: " << name << " to " << tautomer.first
            //           << std::endl;
            const bool is_in = tautomers.find(tsmiles) != tautomers.end();
            if (!is_in) {
              // std::cout << "New tautomer produced: " << tsmiles << std::endl;
              boost::shared_ptr<RWMol> kekulized_product(new RWMol(*wproduct));
              tautomers[tsmiles] = wproduct;
              MolOps::Kekulize(*kekulized_product, false);
              kekulized_mols[tsmiles] = kekulized_product;

              // std::cout << "Now completed: " << std::endl;
              // for (const auto &tautomer : tautomers) {
              //   std::cout << tautomer.first << std::endl;
              // }

            } else {
              // std::cout << "Previous tautomer produced again: " << tsmiles
              //           << std::endl;
            }
          }
        }
      }
      done.push_back(tautomer.first);
    }
    if (tautomers.size() == done.size()) {
      broken = true;
      break;
    }
  }  // while
  if (!broken) {
    BOOST_LOG(rdWarningLog) << "Tautomer enumeration stopped at maximum "
                            << MAX_TAUTOMERS << std::endl;
  }

  // Clean up stereochemistry
  for (auto &tautomer : tautomers) {
    auto &tmp = tautomer.second;
    MolOps::assignStereochemistry(*tmp, true, true);
    //		for (auto &bond : (tmp)->getBonds()) {
    for (size_t i = 0; i < tmp->getNumBonds(); ++i) {
      Bond *bond = (tmp)->getBondWithIdx(i);
      if (bond->getBondType() == 2 &&
          bond->getStereo() > Bond::BondStereo::STEREOANY) {
        unsigned int begin = bond->getBeginAtomIdx();
        unsigned int end = bond->getEndAtomIdx();
        for (auto &other_tautomer : tautomers) {
          auto &other_tmp = other_tautomer.second;
          if (!(other_tmp->getBondBetweenAtoms(begin, end)->getBondType() ==
                2)) {
            Atom *begin_at = tmp->getAtomWithIdx(begin);
            ROMol::OEDGE_ITER beg, end;
            boost::tie(beg, end) = tmp->getAtomBonds(begin_at);
            // std::cout << "BEG " << std::endl;
            // std::cout << *beg << std::endl;
            while (beg != end) {
              Bond::BondDir bonddir = (*tmp)[*beg]->getBondDir();
              if (bonddir == Bond::BondDir::ENDUPRIGHT ||
                  bonddir == Bond::BondDir::ENDDOWNRIGHT) {
                (*tmp)[*beg]->setBondDir(Bond::BondDir::NONE);
              }

              ++beg;
            }
            MolOps::assignStereochemistry(*tmp, true, true);
            // std::cout << "Removed stereochemistry from unfixed double bond"
            //           << std::endl;
            break;
          }
        }
      }
    }
  }

  // get vector of enumerated smiles
  std::vector<ROMOL_SPTR> res;
  for (const auto &tautomer : tautomers) {
    res.push_back(tautomer.second);
    // std::cout << MolToSmiles(*(tautomer.second)) << std::endl;
  }
  return res;
}

}  // namespace MolStandardize
}  // namespace RDKit
