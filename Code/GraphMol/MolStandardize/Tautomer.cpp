//
//  Copyright (C) 2018 Susan H. Leung
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
#include <algorithm>

using namespace RDKit;

namespace RDKit {

namespace MolStandardize {

unsigned int MAX_TAUTOMERS = 1000;

ROMol *TautomerCanonicalizer::canonicalize(const ROMol &mol,
                                           TautomerCatalog *tautcat) {
  PRECONDITION(tautcat, "tautcat not provided");
  // REVIEW: I think it's a good idea to raise an error if this unfinished code
  // is called.
  UNDER_CONSTRUCTION("Tautomer canonicalization not yet implemented");
  TautomerEnumerator tenum;
  std::vector<ROMOL_SPTR> tautomers = tenum.enumerate(mol, tautcat);
  if (tautomers.size() == 1) {
    return tautomers[0].get();
  }
  // Calculate score for each tautomer
  for (const auto t : tautomers) {
    std::string smiles = MolToSmiles(*t);
    std::cout << "Tautomer: " << smiles << std::endl;
    // unsigned int score = 0;
    // Add aromatic ring scores
    VECT_INT_VECT rings;
    MolOps::symmetrizeSSSR(*t, rings);
    // for (const auto ring : rings) {
    //   for (const auto pair : MolStandardize::pairwise(ring)) {
    //     			std::cout << pair.first << " " <<
    //     pair.second
    //     << std::endl;
    //     Bond::BondType btype =
    //         t->getBondBetweenAtoms(pair.first, pair.second)->getBondType();
    //     			std::cout << btype << std::endl;
    //     Stopping for the moment to do the Python wrap
    //   }
    // }
  }

  return new ROMol(mol);
}

std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(
    const ROMol &mol, TautomerCatalog *tautcat) {
  // std::cout << "**********************************" << std::endl;

  PRECONDITION(tautcat, "");
  const TautomerCatalogParams *tautparams = tautcat->getCatalogParams();

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
    // std::map automaticaly sorts tautomers into alphabetical order (SMILES)
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
            // std::cout << "transform mol: " << MolToSmiles(*(transform.Mol))
            //           << std::endl;
            //
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
                MolStandardize::pairwise(idx_matches);
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

}  // namespace MolStandardize
}  // namespace RDKit
