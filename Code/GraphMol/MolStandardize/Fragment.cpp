//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Fragment.h"
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <RDGeneral/types.h>

namespace RDKit {
namespace MolStandardize {

//constructor
FragmentRemover::FragmentRemover(){
	BOOST_LOG(rdInfoLog) << "Initializing FragmentRemover\n" ;
  FragmentCatalogParams fparams(defaultCleanupParameters.fragmentFile);
//  unsigned int numfg = fparams->getNumFuncGroups();
//  TEST_ASSERT(fparams->getNumFuncGroups() == 61);
  this->d_fcat = new FragmentCatalog(&fparams);
	this->LEAVE_LAST = true; 
}

//overloaded constructor
FragmentRemover::FragmentRemover(const std::string fragmentFile, const bool leave_last){
  FragmentCatalogParams fparams(fragmentFile);
  this->d_fcat = new FragmentCatalog(&fparams);
	this->LEAVE_LAST = leave_last; 
}

//Destructor
FragmentRemover::~FragmentRemover(){
	delete d_fcat;
};

ROMol *FragmentRemover::remove(const ROMol &mol) {
	BOOST_LOG(rdInfoLog) << "Running FragmentRemover\n" ;
  PRECONDITION(this->d_fcat, "");
  const FragmentCatalogParams *fparams = this->d_fcat->getCatalogParams();

  PRECONDITION(fparams, "");

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  auto *removed = new ROMol(mol);

  for (auto &fgci : fgrps) {
    std::vector<boost::shared_ptr<ROMol>> frags = MolOps::getMolFrags(*removed);
    // If nothing is left or leave_last and only one fragment, end here
    if (removed->getNumAtoms() == 0 ||
        (this->LEAVE_LAST && frags.size() <= 1)) {
      break;
    }

    std::string fname;
    fgci->getProp(common_properties::_Name, fname);
    ROMol *tmp = RDKit::deleteSubstructs(*removed, *fgci, true);

    if (tmp->getNumAtoms() != removed->getNumAtoms()) {
      BOOST_LOG(rdInfoLog) << "Removed fragment: " << fname << "\n";
    }

    if (this->LEAVE_LAST && tmp->getNumAtoms() == 0) {
      // All the remaining fragments match this pattern - leave them all
			delete tmp;
      break;
    }
		delete removed;
    removed = tmp;
  }
  return removed;
}

bool isOrganic(const ROMol &frag) {
  // Returns true if fragment contains at least one carbon atom.
  for (const auto at : frag.atoms()) {
    if (at->getAtomicNum() == 6) {
      return true;
    }
  }
  return false;
}

LargestFragmentChooser::LargestFragmentChooser(
    const LargestFragmentChooser &other) {
	BOOST_LOG(rdInfoLog) << "Initializing LargestFragmentChooser\n";
  PREFER_ORGANIC = other.PREFER_ORGANIC;
}

ROMol *LargestFragmentChooser::choose(const ROMol &mol) {
	BOOST_LOG(rdInfoLog) << "Running LargestFragmentChooser\n";

  std::vector<boost::shared_ptr<ROMol>> frags = MolOps::getMolFrags(mol);
  LargestFragmentChooser::Largest l;

  for (const auto &frag : frags) {
    std::string smiles = MolToSmiles(*frag);
		BOOST_LOG(rdInfoLog) << "Fragment: " << smiles << "\n";
    bool organic = isOrganic(*frag);
    if (this->PREFER_ORGANIC) {
      // Skip this fragment if not organic and we already have an organic
      // fragment as the largest so far
      if (l.Fragment != nullptr && l.Organic && !organic) continue;
      // Reset largest if it wasn't organic and this fragment is organic
      // if largest and organic and not largest['organic']:
      if (l.Fragment != nullptr && organic && !l.Organic) {
        l.Fragment = nullptr;
      }
    }
    unsigned int numatoms = 0;
    for (const auto at : frag->atoms()) {
      numatoms += 1 + at->getTotalNumHs();
    }
    // Skip this fragment if fewer atoms than the largest
    if (l.Fragment != nullptr && (numatoms < l.NumAtoms)) continue;

    // Skip this fragment if equal number of atoms but weight is lower
    double weight = Descriptors::calcExactMW(*frag);
    if (l.Fragment != nullptr && (numatoms == l.NumAtoms) &&
        (weight < l.Weight))
      continue;

    // Skip this fragment if equal number of atoms and equal weight but smiles
    // comes last alphabetically
    if (l.Fragment != nullptr && (numatoms == l.NumAtoms) &&
        (weight == l.Weight) && (smiles > l.Smiles))
      continue;

		BOOST_LOG(rdInfoLog) << "New largest fragment: " << smiles << " (" <<
				numatoms << ")\n";
    // Otherwise this is the largest so far
    l.Smiles = smiles;
    l.Fragment = frag;
    l.NumAtoms = numatoms;
    l.Weight = weight;
    l.Organic = organic;
  }

  return new ROMol(*(l.Fragment));
}

LargestFragmentChooser::Largest::Largest()
    : Smiles(""), Fragment(nullptr), NumAtoms(0), Weight(0), Organic(false) {}

LargestFragmentChooser::Largest::Largest(
    std::string &smiles, const boost::shared_ptr<ROMol> &fragment,
    unsigned int &numatoms, double &weight, bool &organic)
    : Smiles(smiles),
      Fragment(fragment),
      NumAtoms(numatoms),
      Weight(weight),
      Organic(organic) {}

}  // namespace MolStandardize
}  // namespace RDKit
