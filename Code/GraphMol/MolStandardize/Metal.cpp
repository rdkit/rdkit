//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Metal.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>

using namespace std;
using namespace RDKit;
namespace RDKit {
class RWMol;

class ROMol;

namespace MolStandardize {
MetalDisconnector::MetalDisconnector()
    : metal_nof(
          SmartsToMol("[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,"
                      "Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,"
                      "W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[#7,#8,F]")),
      metal_non(SmartsToMol(
          "[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,"
          "W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]")) {
  BOOST_LOG(rdInfoLog) << "Initializing MetalDisconnector\n";
};

MetalDisconnector::MetalDisconnector(const MetalDisconnector &other) {
  metal_nof = other.metal_nof;
  metal_non = other.metal_non;
};

MetalDisconnector::~MetalDisconnector(){};

ROMol *MetalDisconnector::getMetalNof() { return metal_nof.get(); }

ROMol *MetalDisconnector::getMetalNon() { return metal_non.get(); }

void MetalDisconnector::setMetalNof(const ROMol &mol) {
  this->metal_nof.reset(new ROMol(mol));
}

void MetalDisconnector::setMetalNon(const ROMol &mol) {
  this->metal_non.reset(new ROMol(mol));
}

ROMol *MetalDisconnector::disconnect(const ROMol &mol) {
  auto *res = new RWMol(mol);
  MetalDisconnector::disconnect(*res);
  return static_cast<ROMol *>(res);
}

void MetalDisconnector::disconnect(RWMol &mol) {
  BOOST_LOG(rdInfoLog) << "Running MetalDisconnector\n";
  std::list<ROMOL_SPTR> metalList = {metal_nof, metal_non};
  std::map<int, NonMetal> nonMetals;
  std::map<int, int> metalChargeExcess;
  for (auto &query : metalList) {
    std::vector<MatchVectType> matches;
    SubstructMatch(mol, *query, matches);

    for (const auto &match : matches) {
      int metal_idx = match[0].second;
      metalChargeExcess[metal_idx] = 0;
      auto metal = mol.getAtomWithIdx(metal_idx);
      int non_idx = match[1].second;
      auto nonMetal = mol.getAtomWithIdx(non_idx);

      Bond *b = mol.getBondBetweenAtoms(metal_idx, non_idx);
      int order = (b->getBondType() >= Bond::DATIVEONE &&
                   b->getBondType() <= Bond::DATIVER)
                      ? 0
                      : static_cast<int>(b->getBondTypeAsDouble());
      // disconnecting metal-R bond
      mol.removeBond(metal_idx, non_idx);
      // increment the cut bond count for this non-metal atom
      // and store the metal it was bonded to. We will need this
      // later to adjust the metal charge
      auto &value = nonMetals[non_idx];
      value.cutBonds += order;
      auto it = std::lower_bound(value.boundMetalIndices.begin(),
                                 value.boundMetalIndices.end(), metal_idx);
      if (it == value.boundMetalIndices.end() || *it != metal_idx) {
        value.boundMetalIndices.insert(it, metal_idx);
      }

      BOOST_LOG(rdInfoLog) << "Removed covalent bond between "
                           << metal->getSymbol() << " and "
                           << nonMetal->getSymbol() << "\n";
    }
    //	std::cout << "After removing bond and charge adjustment: " <<
    // MolToSmiles(mol) << std::endl;
  }

  for (auto it = nonMetals.begin(); it != nonMetals.end(); ++it) {
    auto a = mol.getAtomWithIdx(it->first);
    // do not blindly trust the original formal charge as it is often wrong
    // instead, find out the most appropriate formal charge based
    // on the allowed element valences and on its current valence/lone electrons
    // if there are no standard valences we assume the original
    // valence was correct and the non-metal element was neutral
    int valenceBeforeCut = a->getTotalValence();
    int radBeforeCut = a->getNumRadicalElectrons();
    int fcAfterCut = -it->second.cutBonds;
    int valenceAfterCut = 0;
    int loneElectrons = 0;
    const auto &valens =
        PeriodicTable::getTable()->getValenceList(a->getAtomicNum());
    if (!valens.empty() && valens.front() != -1) {
      for (auto v = valens.begin(); v != valens.end(); ++v) {
        valenceAfterCut = valenceBeforeCut + radBeforeCut - it->second.cutBonds;
        if (valenceAfterCut > *v) {
          if (v + 1 != valens.end()) {
            continue;
          }
          valenceAfterCut = *v;
          break;
        }
        fcAfterCut = valenceAfterCut - *v;
        // if there were radicals before and now we have
        // a negative formal charge, then it's a carbene-like
        // system (e.g., [Me]-[C]=O), or there was something silly
        // such as [Me]-[S]=X or [Me]-[P](X)(X)X
        if ((radBeforeCut % 2) && fcAfterCut < 0) {
          ++fcAfterCut;
          ++loneElectrons;
          // no radical doublets on N and higher
          a->setNumRadicalElectrons(a->getAtomicNum() < 7 ? radBeforeCut + 1
                                                          : 0);
        }
        break;
      }
    }
    // do not put a negative charge on sp2 carbon
    if (fcAfterCut == -1 &&
        valenceAfterCut == static_cast<int>(a->getTotalDegree()) + 1) {
      fcAfterCut = 0;
    }
    a->setFormalCharge(fcAfterCut);
    if (!it->second.cutBonds ||
        (fcAfterCut == -1 && a->getAtomicNum() == 6 &&
         valenceAfterCut == static_cast<int>(a->getTotalDegree()) + 2)) {
      // do not take electrons from the metal if it was a dative bond
      // (e.g., [C-]#[O+] coordinated to metal)
      fcAfterCut = 0;
    }
    std::sort(it->second.boundMetalIndices.begin(),
              it->second.boundMetalIndices.end(),
              [&metalChargeExcess](int a, int b) {
                return (metalChargeExcess.at(a) < metalChargeExcess.at(b));
              });
    fcAfterCut += loneElectrons;
    while (fcAfterCut < 0) {
      for (auto i : it->second.boundMetalIndices) {
        // if the bond was not dative, the non-metal stole electrons
        // from the metal(s), so we need to take electrons from
        // once-bonded metal(s)
        if (fcAfterCut++ >= 0) {
          break;
        }
        ++metalChargeExcess[i];
      }
    }
    a->updatePropertyCache();
  }
  // adjust formal charges of metal atoms
  for (auto it = metalChargeExcess.begin(); it != metalChargeExcess.end();
       ++it) {
    auto a = mol.getAtomWithIdx(it->first);
    auto currentFc = a->getFormalCharge();
    const auto &valens =
        PeriodicTable::getTable()->getValenceList(a->getAtomicNum());
    // valens should have at least -1 in it, as the atom data is currently
    // configured, so max_element should never return valens.end().
    auto max_valence = *std::max_element(valens.begin(), valens.end());
    // Don't go over the maximum real valence.
    if (max_valence != -1 && currentFc >= max_valence) {
      continue;
    }
    int fcAfterCut = it->second;
    if (currentFc > 0) {
      // if the original formal charge on the metal was positive, we trust it
      // and add it to the charge excess
      fcAfterCut += currentFc;
    }
    if (!valens.empty() && valens.front() != -1) {
      for (auto v = valens.begin(); v != valens.end(); ++v) {
        // Some metals (e.g. Mg and Ba) have -1 as a final catchall, which
        // is unhelpful for this.
        if (*v == -1) {
          continue;
        }
        if (fcAfterCut > *v) {
          auto next = v + 1;
          if (next != valens.end() && fcAfterCut >= *v) {
            continue;
          }
          fcAfterCut = *v;
          break;
        }
      }
    }
    if (fcAfterCut > currentFc) {
      a->setFormalCharge(fcAfterCut);
    }
    // make sure that radical electrons on metals are 0
    // and are not added to metals by sanitization
    // by setting NoImplicit to false
    a->setNumRadicalElectrons(0);
    a->setNumExplicitHs(0);
    a->setNoImplicit(false);
    a->updatePropertyCache();
  }
}
}  // namespace MolStandardize
}  // namespace RDKit
