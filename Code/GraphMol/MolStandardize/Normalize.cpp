//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Normalize.h"
#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

using namespace std;
using namespace RDKit;

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

typedef boost::flyweight<
    boost::flyweights::key_value<std::string, TransformCatalogParams>,
    boost::flyweights::no_tracking>
    param_flyweight;

// unsigned int MAX_RESTARTS = 200;

// constructor
Normalizer::Normalizer() {
  BOOST_LOG(rdInfoLog) << "Initializing Normalizer\n";
  const TransformCatalogParams *tparams =
      &(param_flyweight(defaultCleanupParameters.normalizations).get());
  //  unsigned int ntransforms = tparams->getNumTransformations();
  //  TEST_ASSERT(ntransforms == 22);
  this->d_tcat = new TransformCatalog(tparams);
  this->MAX_RESTARTS = 200;
}

// overloaded constructor
Normalizer::Normalizer(const std::string normalizeFile,
                       const unsigned int maxRestarts) {
  BOOST_LOG(rdInfoLog) << "Initializing Normalizer\n";
  const TransformCatalogParams *tparams =
      &(param_flyweight(normalizeFile).get());
  this->d_tcat = new TransformCatalog(tparams);
  this->MAX_RESTARTS = maxRestarts;
}

// overloaded constructor
Normalizer::Normalizer(std::istream &normalizeStream,
                       const unsigned int maxRestarts) {
  BOOST_LOG(rdInfoLog) << "Initializing Normalizer\n";
  TransformCatalogParams tparams(normalizeStream);
  this->d_tcat = new TransformCatalog(&tparams);
  this->MAX_RESTARTS = maxRestarts;
}

// destructor
Normalizer::~Normalizer() { delete d_tcat; }

ROMol *Normalizer::normalize(const ROMol &mol) {
  BOOST_LOG(rdInfoLog) << "Running Normalizer\n";
  PRECONDITION(this->d_tcat, "");
  const TransformCatalogParams *tparams = this->d_tcat->getCatalogParams();

  PRECONDITION(tparams, "");
  const std::vector<std::shared_ptr<ChemicalReaction>> &transforms =
      tparams->getTransformations();
  bool sanitizeFrags = false;
  std::vector<boost::shared_ptr<ROMol>> frags =
      MolOps::getMolFrags(mol, sanitizeFrags);
  std::vector<ROMOL_SPTR> nfrags;  //( frags.size() );
  for (const auto &frag : frags) {
    frag->updatePropertyCache(false);
    ROMOL_SPTR nfrag(this->normalizeFragment(*frag, transforms));
    nfrags.push_back(nfrag);
  }
  ROMol *outmol = new ROMol(*(nfrags.back()));
  nfrags.pop_back();
  for (const auto &nfrag : nfrags) {
    ROMol *tmol = combineMols(*outmol, *nfrag);
    delete outmol;
    outmol = tmol;
    //		delete nfrag;
  }
  return outmol;
}

boost::shared_ptr<ROMol> Normalizer::normalizeFragment(
    const ROMol &mol,
    const std::vector<std::shared_ptr<ChemicalReaction>> &transforms) {
  boost::shared_ptr<ROMol> nfrag(new ROMol(mol));
  MolOps::fastFindRings(
      *nfrag);  // this doesn't do anything if rings are already there
  for (unsigned int i = 0; i < MAX_RESTARTS; ++i) {
    bool loop_brake = false;
    // Iterate through Normalization transforms and apply each in order
    for (auto &transform : transforms) {
      boost::shared_ptr<ROMol> product =
          this->applyTransform(nfrag, *transform);
      if (product != nullptr) {
        BOOST_LOG(rdInfoLog)
            << "Rule applied: "
            << transform->getProp<std::string>(common_properties::_Name)
            << "\n";
        nfrag = product;
        loop_brake = true;
        break;
      }
    }
    // For loop finishes normally, all applicable transforms have been applied
    if (!loop_brake) {
      return nfrag;
    }
  }
  BOOST_LOG(rdInfoLog) << "Gave up normalization after " << MAX_RESTARTS
                       << " restarts.\n";
  return nfrag;
}

boost::shared_ptr<ROMol> Normalizer::applyTransform(
    const boost::shared_ptr<ROMol> mol, ChemicalReaction &transform) {
  // Repeatedly apply normalization transform to molecule until no changes
  // occur.
  //
  // It is possible for multiple products to be produced when a rule is applied.
  // The rule is applied repeatedly to each of the products, until no further
  // changes occur or after 20 attempts.
  //
  // If there are multiple unique products after the final application, the
  // first product (sorted alphabetically by SMILES) is chosen.

  MOL_SPTR_VECT mols;
  mols.push_back(mol);

  if (!transform.isInitialized()) transform.initReactantMatchers();
  // REVIEW: what's the source of the 20 in the next line?
  for (unsigned int i = 0; i < 20; ++i) {
    std::vector<Normalizer::Product> pdts;
    for (auto &m : mols) {
      std::vector<MOL_SPTR_VECT> products = transform.runReactants({m});
      for (auto &pdt : products) {
        // shared_ptr<ROMol> p0( new RWMol(*pdt[0]) );
        //				std::cout << MolToSmiles(*p0) <<
        // std::endl;
        unsigned int failed;
        try {
          RWMol *tmol = static_cast<RWMol *>(pdt[0].get());
          // we'll allow atoms with a valence that's too high to make it
          // through, but we should fail if we just created something that
          // can't, for example, be kekulized.
          unsigned int sanitizeOps = MolOps::SANITIZE_ALL ^
                                     MolOps::SANITIZE_CLEANUP ^
                                     MolOps::SANITIZE_PROPERTIES;
          MolOps::sanitizeMol(*tmol, failed, sanitizeOps);
          // REVIEW: is it actually important that we use canonical SMILES here?
          Normalizer::Product np(MolToSmiles(*tmol), pdt[0]);
          pdts.push_back(np);
        } catch (MolSanitizeException &) {
          BOOST_LOG(rdInfoLog) << "FAILED sanitizeMol.\n";
        }
      }
    }
    if (pdts.size() != 0) {
      std::sort(pdts.begin(), pdts.end());
      mols.clear();
      mols.push_back(pdts[0].Mol);
    } else {
      if (i > 0) {
        return mols[0];
      } else {
        return nullptr;
      }
    }
  }
  if (mols.size())
    return mols[0];
  else
    return nullptr;
}

}  // namespace MolStandardize
}  // namespace RDKit
