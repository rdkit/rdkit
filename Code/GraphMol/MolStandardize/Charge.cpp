//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Charge.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/range/adaptor/reversed.hpp>

namespace RDKit {
namespace MolStandardize {

// The default list of ChargeCorrections.
std::vector<ChargeCorrection> CHARGE_CORRECTIONS = {
    ChargeCorrection("[Li,Na,K]", "[Li,Na,K;X0+0]", 1),
    ChargeCorrection("[Mg,Ca]", "[Mg,Ca;X0+0]", 2),
    ChargeCorrection("[Cl]", "[Cl;X0+0]", -1)};

// constructor
Reionizer::Reionizer() {
  AcidBaseCatalogParams abparams(defaultCleanupParameters.acidbaseFile);
  this->d_abcat = new AcidBaseCatalog(&abparams);
  this->d_ccs = CHARGE_CORRECTIONS;
}

Reionizer::Reionizer(const std::string acidbaseFile) {
  AcidBaseCatalogParams abparams(acidbaseFile);
  this->d_abcat = new AcidBaseCatalog(&abparams);
  this->d_ccs = CHARGE_CORRECTIONS;
}

Reionizer::Reionizer(const std::string acidbaseFile,
                     const std::vector<ChargeCorrection> ccs) {
  AcidBaseCatalogParams abparams(acidbaseFile);
  this->d_abcat = new AcidBaseCatalog(&abparams);
  this->d_ccs = ccs;
}

Reionizer::~Reionizer() { delete d_abcat; }

// Reionizer::Reionizer(const AcidBaseCatalog *abcat, const
// std::vector<ChargeCorrection> ccs = CHARGE_CORRECTIONS) 	:
// d_abcat(abcat),
// d_css(css) {};

ROMol *Reionizer::reionize(const ROMol &mol) {
  PRECONDITION(this->d_abcat, "");
  const AcidBaseCatalogParams *abparams = this->d_abcat->getCatalogParams();

  PRECONDITION(abparams, "");
  const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> abpairs =
      abparams->getPairs();

  ROMOL_SPTR omol(new ROMol(mol));
  int start_charge = MolOps::getFormalCharge(*omol);

  for (const auto &cc : this->d_ccs) {
    std::vector<MatchVectType> res;
    ROMOL_SPTR ccmol(SmartsToMol(cc.Smarts));
    unsigned int matches = SubstructMatch(*omol, *ccmol, res);
    if (matches) {
      for (const auto &match : res) {
        for (const auto &pair : match) {
          auto idx = pair.second;
          Atom *atom = omol->getAtomWithIdx(idx);
          BOOST_LOG(rdInfoLog)
              << "Applying charge correction " << cc.Name << " "
              << atom->getSymbol() << " " << cc.Charge << "\n";
          atom->setFormalCharge(cc.Charge);
        }
      }
    }
  }
  int current_charge = MolOps::getFormalCharge(*omol);
  int charge_diff = current_charge - start_charge;
  // std::cout << "Current charge: " << current_charge << std::endl;
  // std::cout << "Charge diff: " << charge_diff << std::endl;

  // If molecule is now neutral, assume everything is now fixed
  // But otherwise, if charge has become more positive,
  // look for additional protonated acid groups to ionize

  if (current_charge != 0) {
    while (charge_diff > 0) {
      // returns the acid strength ranking (ppos)
      // and the substruct match (poccur) in a pair
      std::shared_ptr<std::pair<unsigned int, std::vector<unsigned int>>> res(
          this->strongestProtonated(mol, abpairs));
      if (res == nullptr) {
        break;
      } else {
        unsigned int ppos = res->first;
        std::vector<unsigned int> poccur = res->second;
        std::string abname;
        std::pair<ROMOL_SPTR, ROMOL_SPTR> abpair = abpairs[ppos];
        (abpair.first)->getProp(common_properties::_Name, abname);
        BOOST_LOG(rdInfoLog) << "Ionizing " << abname
                             << " to balance previous charge corrections\n";
        Atom *patom = omol->getAtomWithIdx(poccur.back());
        patom->setFormalCharge(patom->getFormalCharge() - 1);

        if (patom->getNumExplicitHs() > 0) {
          patom->setNumExplicitHs(patom->getNumExplicitHs() - 1);
        }

        patom->updatePropertyCache();
        --charge_diff;
      }
    }
  }

  // std::cout << MolToSmiles(*omol) << std::endl;
  // std::cout << "Charge diff: " << charge_diff << std::endl;

  std::set<std::vector<unsigned int>> already_moved;
  while (true) {
    std::shared_ptr<std::pair<unsigned int, std::vector<unsigned int>>> sp_res(
        this->strongestProtonated(*omol, abpairs));
    std::shared_ptr<std::pair<unsigned int, std::vector<unsigned int>>> wi_res(
        this->weakestIonized(*omol, abpairs));
    if (sp_res != nullptr && wi_res != nullptr) {
      unsigned int ppos = sp_res->first;
      unsigned int ipos = wi_res->first;
      std::vector<unsigned int> poccur = sp_res->second;
      std::vector<unsigned int> ioccur = wi_res->second;
      if (ppos < ipos) {
        if (poccur.back() == ioccur.back()) {
          // Bad! H wouldn't be moved, resulting in infinite loop.
          BOOST_LOG(rdInfoLog)
              << "Aborted reionization due to unexpected situation\n";
          break;
        }

        std::vector<unsigned int> key = {poccur.back(), ioccur.back()};
        std::sort(key.begin(), key.end());
        const bool is_in = already_moved.find(key) != already_moved.end();
        if (is_in) {
          BOOST_LOG(rdInfoLog)
              << "Aborting reionization to avoid infinite loop due \
								to it being ambiguous where to put a Hydrogen\n";
          break;
        }
        already_moved.insert(key);

        std::string prot_name, ionized_name;
        std::pair<ROMOL_SPTR, ROMOL_SPTR> prot_pair = abpairs[ppos];
        std::pair<ROMOL_SPTR, ROMOL_SPTR> ionized_pair = abpairs[ipos];
        (prot_pair.first)->getProp(common_properties::_Name, prot_name);
        (ionized_pair.first)->getProp(common_properties::_Name, ionized_name);

        BOOST_LOG(rdInfoLog) << "Moved proton from " << prot_name << " to "
                             << ionized_name << "\n";
        // Remove hydrogen from strongest protonated
        Atom *patom = omol->getAtomWithIdx(poccur.back());
        patom->setFormalCharge(patom->getFormalCharge() - 1);
        // If no implicit Hs to autoremove, and at least 1 explicit H to remove,
        //  reduce explicit count by 1
        if (patom->getNumImplicitHs() == 0 && patom->getNumExplicitHs() > 0) {
          patom->setNumExplicitHs(patom->getNumExplicitHs() - 1);
          // TODO: Remove any chiral label on patom?
        }
        patom->updatePropertyCache();

        // Add hydrogen to weakest ionized
        Atom *iatom = omol->getAtomWithIdx(ioccur.back());
        iatom->setFormalCharge(iatom->getFormalCharge() + 1);
        // Increase explicit H count if no implicit, or aromatic N or P,
        // or non default valence state
        const PeriodicTable *table = PeriodicTable::getTable();
        INT_VECT valence_list = table->getValenceList(iatom->getAtomicNum());
        bool found =
            (std::find(valence_list.begin(), valence_list.end(),
                       iatom->getTotalValence()) != valence_list.end());
        if (iatom->getNoImplicit() ||
            ((patom->getAtomicNum() == 7 || patom->getAtomicNum() == 15) &&
             patom->getIsAromatic()) ||
            !found) {
          iatom->setNumExplicitHs(iatom->getNumExplicitHs() + 1);
        }
        iatom->updatePropertyCache();
      } else {
        break;
      }
    } else {
      break;
    }
  }  // while loop

  RWMOL_SPTR wmol(new RWMol(*omol));
  MolOps::sanitizeMol(*wmol);

  return new ROMol(*wmol);
}

std::pair<unsigned int, std::vector<unsigned int>>
    *Reionizer::strongestProtonated(
        const ROMol &mol,
        const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs) {
  // position is the position in the acid list.
  unsigned int position = 0;
  for (const auto &abpair : abpairs) {
    RDKit::MatchVectType res;
    unsigned int matches = SubstructMatch(mol, *(abpair.first), res);
    if (matches > 0) {
      std::vector<unsigned int> occurence;
      for (const auto &pair : res) {
        occurence.push_back(pair.second);
      }
      return new std::pair<unsigned int, std::vector<unsigned int>>(position,
                                                                    occurence);
    }
    ++position;
  }
  return nullptr;
}

std::pair<unsigned int, std::vector<unsigned int>> *Reionizer::weakestIonized(
    const ROMol &mol,
    const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs) {
  // position is the position in the acid list.
  unsigned int position = 0;
  for (const auto &abpair : boost::adaptors::reverse(abpairs)) {
    RDKit::MatchVectType res;
    unsigned int matches = SubstructMatch(mol, *(abpair.second), res);
    if (matches > 0) {
      std::vector<unsigned int> occurence;
      for (const auto &pair : res) {
        occurence.push_back(pair.second);
      }
      return new std::pair<unsigned int, std::vector<unsigned int>>(
          (abpairs.size() - position - 1), occurence);
    }
    ++position;
  }
  return nullptr;
}

Uncharger::Uncharger()
    : pos_h(SmartsToMol("[+!H0!$(*~[-])]")),
      pos_quat(SmartsToMol("[+H0!$(*~[-])]")),
      neg(SmartsToMol("[-!$(*~[+H0])]")),
      neg_acid(SmartsToMol("[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]")){};

Uncharger::Uncharger(const Uncharger &other) {
  pos_h = other.pos_h;
  pos_quat = other.pos_quat;
  neg = other.neg;
  neg_acid = other.neg_acid;
};

Uncharger::~Uncharger(){};

ROMol *Uncharger::uncharge(const ROMol &mol) {
  BOOST_LOG(rdInfoLog) << "Running Uncharger\n";
  ROMol *omol = new ROMol(mol);

  std::vector<MatchVectType> p_matches;
  std::vector<MatchVectType> q_matches;
  std::vector<MatchVectType> n_matches;
  std::vector<MatchVectType> a_matches;

  // Get atom ids for matches
  SubstructMatch(*omol, *(this->pos_h), p_matches);
  unsigned int q_matched = SubstructMatch(*omol, *(this->pos_quat), q_matches);
  unsigned int n_matched = SubstructMatch(*omol, *(this->neg), n_matches);
  unsigned int a_matched = SubstructMatch(*omol, *(this->neg_acid), a_matches);

  // trying to understand how to use n_matches as a vector...
  //	std::cout << "Size " << n_matches.size() << std::endl;
  //	std::cout << n_matches[0][0].second << std::endl;
  //	std::cout << n_matches[1][0].second << std::endl;
  //	for (auto &i : n_matches) {
  //		for (auto &j : i) {
  //			std::cout << j.second << std::endl;
  //		}
  //	}
  //

  // Neutralize negative charges
  if (q_matched > 0) {
    // Surplus negative charges more than non-neutralizable positive charges
    int neg_surplus = n_matched - q_matched;
    if (a_matched > 0 && neg_surplus > 0) {
      // zwitterion with more negative charges than quaternary positive centres
      while (neg_surplus > 0 && a_matched > 0) {
        // Add hydrogen to first negative acid atom, increase formal charge
        // Until quaternary positive == negative total or no more negative acid
        Atom *atom = omol->getAtomWithIdx(a_matches[0][0].second);
        atom->setNoImplicit(true);
        a_matches.erase(a_matches.begin());
        atom->setNumExplicitHs(atom->getNumExplicitHs() + 1);
        atom->setFormalCharge(atom->getFormalCharge() + 1);
        --neg_surplus;
        BOOST_LOG(rdInfoLog) << "Removed negative charge.\n";
      }
    }
  } else {
    std::vector<unsigned int> n_idx_matches;
    for (const auto &match : n_matches) {
      for (const auto &pair : match) {
        n_idx_matches.push_back(pair.second);
      }
    }
    for (const auto &idx : n_idx_matches) {
      Atom *atom = omol->getAtomWithIdx(idx);
      atom->setNoImplicit(true);
      while (atom->getFormalCharge() < 0) {
        atom->setNumExplicitHs(atom->getNumExplicitHs() + 1);
        atom->setFormalCharge(atom->getFormalCharge() + 1);
        BOOST_LOG(rdInfoLog) << "Removed negative charge.\n";
      }
    }
  }
  // Neutralize positive charges
  std::vector<unsigned int> p_idx_matches;
  for (const auto &match : p_matches) {
    for (const auto &pair : match) {
      p_idx_matches.push_back(pair.second);
    }
  }
  for (const auto &idx : p_idx_matches) {
    Atom *atom = omol->getAtomWithIdx(idx);
    if (!atom->getNumExplicitHs()) {
      // atoms from places like Mol blocks are normally missing explicit Hs:
      atom->setNumExplicitHs(atom->getTotalNumHs());
    }
    atom->setNoImplicit(true);
    while (atom->getFormalCharge() > 0 && atom->getNumExplicitHs() > 0) {
      atom->setNumExplicitHs(atom->getTotalNumHs() - 1);
      atom->setFormalCharge(atom->getFormalCharge() - 1);
      BOOST_LOG(rdInfoLog) << "Removed positive charge.\n";
    }
  }
  return omol;
}

}  // namespace MolStandardize
}  // namespace RDKit
