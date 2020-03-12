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

namespace RDKit
{
namespace NewCIPLabelling
{

using RDKitCfg = Configuration<RdkA, RdkB>;
using RDKitSeqRule = SequenceRule<RdkA, RdkB>;

namespace
{
std::vector<RDKitCfg*> findConfigs(RDKitCipMol& mol)
{
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

            auto cfg =
                new Tetrahedral<RdkA, RdkB>(atom, std::move(carriers), parity);
            configs.push_back(cfg);
        }
    }

    for (int i = 0; i < mol.getNumBonds(); ++i) {
        RdkB bond = mol.getBond(i);

        int bond_cfg;
        switch (bond->getStereo()) {
        case Bond::STEREOE:
        case Bond::STEREOTRANS:
            bond_cfg = 1;
            break;
        case Bond::STEREOZ:
        case Bond::STEREOCIS:
            bond_cfg = 2;
            break;
        default:
            continue;
        }
        auto atoms =
            std::vector<RdkA>{{bond->getBeginAtom(), bond->getEndAtom()}};
        auto anchors_idx = bond->getStereoAtoms();
        auto anchors = std::vector<RdkA>{
            {mol.getAtom(anchors_idx[0]), mol.getAtom(anchors_idx[1])}};

        auto cfg = new Sp2Bond<RdkA, RdkB>(bond, std::move(atoms),
                                           std::move(anchors), bond_cfg);
        configs.push_back(cfg);
    }

    return configs;
}
} // namespace
} // namespace NewCIPLabelling

namespace MolOps
{
void assignRealCIPStereo(ROMol& mol)
{
    auto cipmol = NewCIPLabelling::RDKitCipMol(&mol);

    auto configs = NewCIPLabelling::findConfigs(cipmol);
    NewCIPLabelling::Labeller::label(&cipmol, configs);

    // Cleanup
    for (auto& c : configs) {
        delete c;
    }
}
} // namespace MolOps
} // end of namespace RDKit
