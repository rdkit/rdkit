//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MonomerMol.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>

namespace RDKit {


bool isMonomer(const Atom* atom) {
    return atom->hasProp(SMILES_MONOMER);
}

std::string getPolymerId(const Atom* atom) {
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }
    return monomer_info->getChainId();
}

unsigned int getResidueNumber(const Atom* atom) {
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }
    return monomer_info->getResidueNumber();
}

namespace {

void setResidueNumber(Atom* atom, int residue_number) {
    auto* monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error(
            "Atom " + std::to_string(atom->getIdx()) + " does not have monomer info");
    }
    monomer_info->setResidueNumber(residue_number);
}

bool isValidChain(const MonomerMol& monomer_mol, std::string_view polymer_id) {
    // Check that the residue ordering is valid for this polymer. The residues
    // should be in connectivity order
    const auto chain = monomer_mol.getPolymer(polymer_id);
    auto last_residue = chain.atoms[0];
    for (size_t i = 1; i < chain.atoms.size(); ++i) {
        const auto bond =
            monomer_mol.getBondBetweenAtoms(last_residue, chain.atoms[i]);
        if (bond == nullptr) {
            return false;
        }

        // Bond direction should be in the same order as residues
        if (getResidueNumber(bond->getBeginAtom()) >
                getResidueNumber(bond->getEndAtom()) &&
            chain.atoms[i] != bond->getEndAtomIdx()) {
            return false;
        }

        auto linkage = bond->getProp<std::string>(LINKAGE);
        if (linkage == BACKBONE_LINKAGE) {
            last_residue = chain.atoms[i];
        } else if (linkage != BRANCH_LINKAGE) {
            return false;
        }
    }
    return true;
}

Atom* findChainBegin(MonomerMol& monomer_mol) {
    // Find the beginning of the chain by starting at an arbirtary atom
    // and following the backbone backwards until the 'source' (beginning of the
    // chain) is found. If the beginning of the chain is in a cycle, then the
    // last discovered atom will be made the beginning of the chain.
    std::vector<bool> seen(monomer_mol.getNumAtoms(), 0);
    auto chain_begin = monomer_mol.getAtomWithIdx(0);
    bool updated = true;
    while (updated) {
        updated = false;
        seen[chain_begin->getIdx()] = true;
        for (const auto bond : monomer_mol.atomBonds(chain_begin)) {
            auto linkage = bond->getProp<std::string>(LINKAGE);
            if (linkage == "R3-R3") {
                // Don't cross cysteine bridges
                continue;
            } else if (seen[bond->getOtherAtom(chain_begin)->getIdx()]) {
                // this is a loop
                continue;
            }
            // Bonds are always oriented so that the beginAtom->endAtom matches
            // the direction of the chain. If this atom is the 'end' of the
            // bond, it is not the beginning of the chain, so follow the parent
            if (bond->getEndAtom() == chain_begin) {
                auto other = bond->getOtherAtom(chain_begin);
                chain_begin = other;
                updated = true;
                break;
            }
        }
    }
    return chain_begin;
}

void orderResidues(MonomerMol& monomer_mol) {
    // Currently assumes that all monomers are in the same chain. We will
    // eventually want to order residues on a per-chain basis.

    // Find the beginning of the chain (must be a backbone monomer)
    auto chain_begin = findChainBegin(monomer_mol);

    // Now re-order the residues beginning at chain_begin
    std::vector<Atom*> queue;
    std::vector<bool> visited(monomer_mol.getNumAtoms(), 0);
    queue.push_back(chain_begin);
    visited[chain_begin->getIdx()] = true;

    int current_res_idx = 1;
    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        setResidueNumber(current, current_res_idx);
        ++current_res_idx;

        // When ordering residues, sidechain monomers should come before
        // backbone monomers, which is more specific then a regular BFS
        // ordering. For example: A.B(X)C should be ordered as A, B, X, C
        // instead of A, B, C, X
        for (const auto bond : monomer_mol.atomBonds(current)) {
            if (bond->getEndAtom() == current ||
                visited[bond->getOtherAtom(current)->getIdx()]) {
                continue;
            }
            auto linkage = bond->getProp<std::string>(LINKAGE);
            if (linkage == BRANCH_LINKAGE) {
                auto other = bond->getOtherAtom(current);
                setResidueNumber(other, current_res_idx);
                ++current_res_idx;
            } else if (linkage == BACKBONE_LINKAGE) {
                queue.push_back(bond->getOtherAtom(current));
            }
        }
    }
}

} // anonymous namespace

void MonomerMol::assignChains() {
    // Currently, orderResidues only works when there is a single chain
    auto chain_ids = getPolymerIds();
    if (chain_ids.size() == 1 && !isValidChain(*this, chain_ids[0])) {
        orderResidues(*this);
    }

    // Determine and mark the 'connection bonds'
    if (!getRingInfo()->isInitialized()) {
        MolOps::findSSSR(*this, nullptr, /*includeDativeBonds = */ true);
    }
    // get atom rings that belong to a single polymer
    const auto& bnd_rings = getRingInfo()->bondRings();

    for (const auto& ring : bnd_rings) {
        if (std::ranges::all_of(ring, [this](const auto& idx) {
                auto begin_at = getBondWithIdx(idx)->getBeginAtom();
                auto end_at = getBondWithIdx(idx)->getEndAtom();
                return getPolymerId(begin_at) == getPolymerId(end_at);
            })) {
            // break this ring -- find bond between residue #s with largest
            // difference
            unsigned int connection_bond = getNumBonds();
            int max_diff = 0;
            for (const auto& idx : ring) {
                const auto bond = getBondWithIdx(idx);
                auto begin_at = bond->getBeginAtom();
                auto end_at = bond->getEndAtom();
                int diff =
                    getResidueNumber(begin_at) - getResidueNumber(end_at);
                if (std::abs(diff) > max_diff) {
                    connection_bond = bond->getIdx();
                    max_diff = std::abs(diff);
                }
            }
            if (connection_bond != getNumBonds()) {
                std::string linkage = getBondWithIdx(connection_bond)
                                          ->getProp<std::string>(LINKAGE);
                getBondWithIdx(connection_bond)
                    ->setProp(EXTRA_LINKAGE, linkage);
            } else {
                // temporary error handling
                throw std::runtime_error("Could not find connection bond");
            }
        }
    }
    for (auto bond : bonds()) {
        if (getPolymerId(bond->getBeginAtom()) !=
            getPolymerId(bond->getEndAtom())) {
            std::string linkage = bond->getProp<std::string>(LINKAGE);
            bond->setProp(EXTRA_LINKAGE, linkage);
        }
    }
}

} // namespace RDKit
