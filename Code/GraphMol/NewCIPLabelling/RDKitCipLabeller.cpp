//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>

#include "Labeller.hpp"
#include "RDKitCipMol.hpp"
#include "configs/Sp2Bond.hpp"
#include "configs/Tetrahedral.hpp"

namespace RDKit {
namespace NewCIPLabelling {

using RDKitCfg = Configuration<RdkA, RdkB>;
using RDKitSeqRule = SequenceRule<RdkA, RdkB>;

namespace {
std::vector<RDKitCfg*> findConfigs(RDKitCipMol& mol) {
  std::vector<RDKitCfg*> configs;

  for (auto& atom : mol.atoms()) {
    auto chiraltag = static_cast<int>(atom->getChiralTag());
    if (chiraltag > 0) {
      auto carriers = std::vector<RdkA>();
      carriers.reserve(4);
      for (RdkB bond : mol.getBonds(atom)) {
        RdkA nbr = mol.getOther(bond, atom);
        carriers.push_back(nbr);
      }
      if (carriers.size() < 4) {
        // Implicit H -- use the central atom instead of a dummy H
        carriers.push_back(atom);
      }

      // RDKit uses SMILES parity, and centres requires SDF parity, which
      // is opposite.
      auto parity = chiraltag ^ 0x3;

      auto cfg = new Tetrahedral<RdkA, RdkB>(atom, std::move(carriers), parity);
      configs.push_back(cfg);
    }
  }

  for (int i = 0; i < mol.getNumBonds(); ++i) {
    RdkB bond = mol.getBond(i);

    if (bond->getBondType() == Bond::DOUBLE) {
      std::vector<RdkA> atoms{{bond->getBeginAtom(), bond->getEndAtom()}};
      std::vector<RdkA> anchors;
      auto previous_bond_dir = Bond::BondDir::NONE;
      auto bond_cfg = 0;
      for (const auto& atom : atoms) {
        for (auto& nbr_bnd : mol.getBonds(atom)) {
          auto nbr_bond_dir = nbr_bnd->getBondDir();
          if (nbr_bond_dir == Bond::BondDir::ENDDOWNRIGHT ||
              nbr_bond_dir == Bond::BondDir::ENDUPRIGHT) {
            anchors.push_back(nbr_bnd->getOtherAtom(atom));
            bond_cfg += nbr_bond_dir != previous_bond_dir;
            std::swap(previous_bond_dir, nbr_bond_dir);
            break;
          }
        }
      }

      if (anchors.size() == 2) {
        auto cfg = new Sp2Bond<RdkA, RdkB>(bond, std::move(atoms),
                                           std::move(anchors), bond_cfg);
        configs.push_back(cfg);
      }
    }
  }

  return configs;
}
}  // namespace
}  // namespace NewCIPLabelling

namespace MolOps {
void assignRealCIPStereo(ROMol& mol) {
  auto cipmol = NewCIPLabelling::RDKitCipMol(&mol);

  auto configs = NewCIPLabelling::findConfigs(cipmol);
  NewCIPLabelling::Labeller::label(&cipmol, configs);

  // Cleanup
  for (auto& c : configs) {
    delete c;
  }
}
}  // namespace MolOps
}  // end of namespace RDKit
