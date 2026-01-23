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

#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolOps.h>

#include <boost/functional/hash.hpp>

namespace RDKit
{

ChainType toChainType(std::string_view chain_type)
{
    if (chain_type == "PEPTIDE") {
        return ChainType::PEPTIDE;
    } else if (chain_type == "RNA") {
        return ChainType::RNA;
    } else if (chain_type == "DNA") {
        return ChainType::DNA;
    } else if (chain_type == "CHEM") {
        return ChainType::CHEM;
    } else {
        throw std::invalid_argument("Invalid chain type");
    }
}

std::string toString(ChainType chain_type)
{
    switch (chain_type) {
        case ChainType::PEPTIDE:
            return "PEPTIDE";
        case ChainType::RNA:
            return "RNA";
        case ChainType::DNA:
            return "DNA";
        case ChainType::CHEM:
            return "CHEM";
        default:
            throw std::invalid_argument("Invalid chain type");
    }
}

namespace {

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
{
    // Parse linkage in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos || linkage.size() < 5 ||
        linkage[0] != 'R' || linkage[dash + 1] != 'R') {
        throw std::runtime_error("Invalid linkage format: " + linkage);
    }
    return {std::stoi(linkage.substr(1, dash - 1)),
            std::stoi(linkage.substr(dash + 2))};
}

}  // namespace

MonomerMol& MonomerMol::operator=(const MonomerMol& other)
{
    if (this != &other) {
        RWMol::operator=(other);
    }
    return *this;
}

MonomerMol& MonomerMol::operator=(MonomerMol&& other) noexcept
{
    if (this != &other) {
        RWMol::operator=(std::move(other));
    }
    return *this;
}

void MonomerMol::addConnection(size_t monomer1, size_t monomer2,
                               const std::string& linkage,
                               const bool is_custom_bond)
{
    // if bond already exists, extend linkage information
    if (auto bond = getBondBetweenAtoms(monomer1, monomer2);
        bond && is_custom_bond) {
        std::string old_linkage;

        // If the linkage property isn't set, something went wrong
        if (!bond->getPropIfPresent(LINKAGE, old_linkage)) {
            throw std::runtime_error(
                "No linkage property on bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // Make sure we're not recreating this same bond.
        if (old_linkage.find(linkage) != std::string::npos) {
            throw std::runtime_error(
                "Can't duplicate " + linkage + " bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // Since hydrogen bonding and covalent bonds are different types of bond
        // serializing them with the same bond type doesn't make too much sense.
        if (linkage.find("pair") == 0 ||
            old_linkage.find("pair") != std::string::npos) {
            throw std::runtime_error(
                "Multiple bonds can't include hydrogen bond for bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // FIXME: For now, don't allow multiple custom bonds between the same
        // two atoms
        if (bond->hasProp(CUSTOM_BOND)) {
            throw std::runtime_error(
                "Multiple custom bonds not supported for bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // Update the linkage property
        bond->setProp(CUSTOM_BOND, linkage);
    } else {
        auto create_bond = [&](size_t first_monomer, size_t second_monomer,
                               ::RDKit::Bond::BondType bond_type,
                               const std::string& linkage_prop) {
            // Use RWMol::addBond directly to bypass our overridden addBond
            const auto new_total =
                RWMol::addBond(first_monomer, second_monomer, bond_type);
            auto new_bond = getBondWithIdx(new_total - 1);
            new_bond->setProp(LINKAGE, linkage_prop);
            if (is_custom_bond) {
                new_bond->setProp(CUSTOM_BOND, linkage_prop);
            }
        };

        // Connections that use specific and different attachment points (such
        // as R2-R1 or R3-R1) should be set as directional (dative) bonds so
        // that substructure matching and canonicalization can take sequence
        // direction into account. To ensure consistent bond directionality, we
        // always set the direction from the higher attachment point to the
        // lower attachment point.
        bool set_directional_bond = false;
        if (linkage.front() != 'p' && linkage.find('?') == std::string::npos) {
            auto [begin_attchpt, end_attchpt] = getAttchpts(linkage);
            if (begin_attchpt > end_attchpt) {
                create_bond(monomer1, monomer2, ::RDKit::Bond::DATIVE, linkage);
                set_directional_bond = true;
            } else if (begin_attchpt < end_attchpt) {
                auto new_linkage = "R" + std::to_string(end_attchpt) + "-R" +
                                   std::to_string(begin_attchpt);
                create_bond(monomer2, monomer1, ::RDKit::Bond::DATIVE,
                            new_linkage);
                set_directional_bond = true;
            }
        }

        if (!set_directional_bond) {
            auto bond_type = (linkage.front() == 'p' ? ::RDKit::Bond::ZERO
                                                     : ::RDKit::Bond::SINGLE);
            create_bond(monomer1, monomer2, bond_type, linkage);
        }
    }

    if (!is_custom_bond && linkage == BRANCH_LINKAGE) {
        // monomer2 is a branch monomer
        getAtomWithIdx(monomer2)->setProp(BRANCH_MONOMER, true);
    }
}

void MonomerMol::addConnection(size_t monomer1, size_t monomer2,
                               ConnectionType connection_type)
{
    switch (connection_type) {
        case ConnectionType::FORWARD:
            addConnection(monomer1, monomer2, BACKBONE_LINKAGE);
            break;
        case ConnectionType::SIDECHAIN:
            addConnection(monomer1, monomer2, BRANCH_LINKAGE);
            break;
        case ConnectionType::CROSSLINK:
            addConnection(monomer1, monomer2, CROSS_LINKAGE);
            break;
    }
}

bool MonomerMol::isMonomer(unsigned int atomIdx) const
{
    const auto* atom = getAtomWithIdx(atomIdx);
    // An atom is a monomer if it has the ATOM_LABEL property
    return atom->hasProp(ATOM_LABEL);
}

unsigned int MonomerMol::addBond(unsigned int beginAtomIdx,
                                 unsigned int endAtomIdx, Bond::BondType order)
{
    bool beginIsMonomer = isMonomer(beginAtomIdx);
    bool endIsMonomer = isMonomer(endAtomIdx);

    if (beginIsMonomer && endIsMonomer) {
        throw std::invalid_argument(
            "Cannot add bond between two monomers using addBond. "
            "Use addConnection instead.");
    }

    if (beginIsMonomer || endIsMonomer) {
        throw std::invalid_argument(
            "Cannot add bond between monomer and atom without specifying "
            "attachment point. Use addBond(monomerIdx, atomIdx, "
            "attachmentPoint) instead.");
    }

    // Both are regular atoms, delegate to RWMol::addBond
    return RWMol::addBond(beginAtomIdx, endAtomIdx, order);
}

unsigned int MonomerMol::addBond(unsigned int monomerIdx, unsigned int atomIdx,
                                 unsigned int attachmentPoint,
                                 Bond::BondType order)
{
    // TODO: Should this all be under addLinkage?
    bool monomerIsMonomer = isMonomer(monomerIdx);
    bool atomIsMonomer = isMonomer(atomIdx);

    if (monomerIsMonomer && atomIsMonomer) {
        throw std::invalid_argument(
            "Both indices are monomers. Use addConnection for monomer-to-monomer "
            "bonds.");
    }

    if (!monomerIsMonomer && !atomIsMonomer) {
        throw std::invalid_argument(
            "Neither index is a monomer. Use addBond for atom-to-atom "
            "bonds.");
    }

    // Ensure monomerIdx is the monomer
    unsigned int actualMonomerIdx = monomerIsMonomer ? monomerIdx : atomIdx;
    unsigned int actualAtomIdx = monomerIsMonomer ? atomIdx : monomerIdx;

    // Create the bond
    const auto new_total = RWMol::addBond(actualMonomerIdx, actualAtomIdx, order);
    auto* bond = getBondWithIdx(new_total - 1);

    // Set the linkage property - use "RX-" format to indicate monomer-to-atom
    std::string linkage = "R" + std::to_string(attachmentPoint) + "-";
    bond->setProp(LINKAGE, linkage);

    return new_total;
}

[[nodiscard]] std::string getPolymerId(const ::RDKit::Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getChainId();
}


[[nodiscard]] std::vector<std::string> getPolymerIds(const RDKit::ROMol& monomer_mol)
{
    std::vector<std::string> polymer_ids;
    for (auto atom : monomer_mol.atoms()) {
        // Skip non-monomer atoms (regular atoms don't have polymer info)
        if (!atom->hasProp(ATOM_LABEL)) {
            continue;
        }
        auto id = getPolymerId(atom);
        // in vector to preserve order of polymers
        if (std::find(polymer_ids.begin(), polymer_ids.end(), id) ==
            polymer_ids.end()) {
            polymer_ids.push_back(id);
        }
    }
    return polymer_ids;
}

[[nodiscard]] unsigned int getResidueNumber(const ::RDKit::Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getResidueNumber();
}

std::unique_ptr<Atom> makeMonomer(std::string_view name,
                                     std::string_view chain_id,
                                     int residue_number, bool is_smiles)
{
    auto a = std::make_unique<::RDKit::Atom>();
    std::string n{name};
    a->setProp(ATOM_LABEL, n);
    a->setProp("Name", n);
    a->setProp("smilesSymbol", n);
    // Always start with BRANCH_MONOMER as false, will be set to
    // true if branch linkage is made to this monomer
    a->setProp(BRANCH_MONOMER, false);
    a->setProp(SMILES_MONOMER, is_smiles);

    // hack to get some level of canonicalization for monomer mols
    static boost::hash<std::string> hasher;
    a->setIsotope(hasher(n));
    a->setNoImplicit(true);

    auto* residue_info = new ::RDKit::AtomPDBResidueInfo();
    residue_info->setResidueNumber(residue_number);
    residue_info->setResidueName(n);
    residue_info->setChainId(std::string{chain_id});
    a->setMonomerInfo(residue_info);
    return a;
}

size_t MonomerMol::addMonomer(std::string_view name, int residue_number,
                              std::string_view chain_id,
                              MonomerType monomer_type)
{
    auto monomer = makeMonomer(name, chain_id, residue_number,
                               monomer_type == MonomerType::SMILES);
    bool update_label = true;
    bool take_ownership = true;
    auto new_index = addAtom(monomer.release(), update_label, take_ownership);
    return new_index;
}

size_t MonomerMol::addMonomer(std::string_view name, MonomerType monomer_type)
{
    if (getNumAtoms() == 0) {
        throw std::invalid_argument(
            "No atoms in molecule to determine chain ID");
    }
    const auto* last_monomer = getAtomWithIdx(getNumAtoms() - 1);
    const auto chain_id = getPolymerId(last_monomer);
    const auto residue_number = getResidueNumber(last_monomer) + 1;
    return addMonomer(name, residue_number, chain_id, monomer_type);
}

void setResidueNumber(RDKit::Atom* atom, int residue_number)
{
    auto* residue_info =
        static_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    if (residue_info == nullptr) {
        throw std::runtime_error(
            "Atom " + std::to_string(atom->getIdx()) + " does not have residue info");
    }
    residue_info->setResidueNumber(residue_number);
}

Chain getPolymer(const RDKit::ROMol& monomer_mol, std::string_view polymer_id)
{
    std::vector<unsigned int> atoms;
    for (auto atom : monomer_mol.atoms()) {
        // Skip non-monomer atoms
        if (!atom->hasProp(ATOM_LABEL)) {
            continue;
        }
        if (getPolymerId(atom) == polymer_id) {
            atoms.push_back(atom->getIdx());
        }
    }
    // Sort by get_residue_num
    std::sort(atoms.begin(), atoms.end(),
              [&monomer_mol](unsigned int a, unsigned int b) {
                  return getResidueNumber(monomer_mol.getAtomWithIdx(a)) <
                         getResidueNumber(monomer_mol.getAtomWithIdx(b));
              });
    std::vector<unsigned int> bonds;
    for (auto bond : monomer_mol.bonds()) {
        // Skip bonds involving non-monomer atoms
        if (!bond->getBeginAtom()->hasProp(ATOM_LABEL) ||
            !bond->getEndAtom()->hasProp(ATOM_LABEL)) {
            continue;
        }
        if (getPolymerId(bond->getBeginAtom()) == polymer_id &&
            getPolymerId(bond->getEndAtom()) == polymer_id) {
            bonds.push_back(bond->getIdx());
        }
    }

    std::string annotation{};
    for (const auto& sg : ::RDKit::getSubstanceGroups(monomer_mol)) {
        if ((sg.getProp<std::string>("TYPE") != "COP") ||
            !sg.hasProp(ANNOTATION) || !sg.hasProp("ID")) {
            continue;
        }
        if (sg.getProp<std::string>("ID") == polymer_id) {
            annotation = sg.getProp<std::string>(ANNOTATION);
            break;
        }
    }
    return {atoms, bonds, annotation};
}

namespace {

bool isValidChain(const RDKit::MonomerMol& monomer_mol, std::string_view polymer_id)
{
    // Check that the residue ordering is valid for this polymer. The residues
    // should be in connectivity order
    const auto chain = getPolymer(monomer_mol, polymer_id);
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

RDKit::Atom* findChainBegin(RDKit::MonomerMol& monomer_mol)
{
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

void orderResidues(RDKit::MonomerMol& monomer_mol)
{
    // Currently assumes that all monomers are in the same chain. We will
    // eventually want to order residues on a per-chain basis.

    // Find the beginning of the chain (must be a backbone monomer)
    auto chain_begin = findChainBegin(monomer_mol);

    // Now re-order the residues beginning at chain_begin
    std::vector<RDKit::Atom*> queue;
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

}  // namespace

void assignChains(MonomerMol& monomer_mol)
{
    monomer_mol.setProp("HELM_MODEL", true);

    // Currently, orderResidues only works when there is a single chain
    auto chain_ids = getPolymerIds(monomer_mol);
    if (chain_ids.size() == 1 && !isValidChain(monomer_mol, chain_ids[0])) {
        orderResidues(monomer_mol);
    }

    // Determine and mark the 'connection bonds'
    // Include dative bonds in ring finding since we use them for directional
    // monomer-to-monomer connections
    constexpr bool include_dative_bonds = true;
    ::RDKit::MolOps::findSSSR(monomer_mol, /*res=*/nullptr, include_dative_bonds);

    // get atom rings that belong to a single polymer
    const auto& bnd_rings = monomer_mol.getRingInfo()->bondRings();

    for (const auto& ring : bnd_rings) {
        if (std::ranges::all_of(ring, [&](const auto& idx) {
                auto begin_at = monomer_mol.getBondWithIdx(idx)->getBeginAtom();
                auto end_at = monomer_mol.getBondWithIdx(idx)->getEndAtom();
                return getPolymerId(begin_at) == getPolymerId(end_at);
            })) {
            // break this ring -- find bond between residue #s with largest
            // difference
            unsigned int connection_bond = monomer_mol.getNumBonds();
            int max_diff = 0;
            for (const auto& idx : ring) {
                const auto bond = monomer_mol.getBondWithIdx(idx);
                auto begin_at = bond->getBeginAtom();
                auto end_at = bond->getEndAtom();
                int diff =
                    getResidueNumber(begin_at) - getResidueNumber(end_at);
                if (std::abs(diff) > max_diff) {
                    connection_bond = bond->getIdx();
                    max_diff = std::abs(diff);
                }
            }
            if (connection_bond != monomer_mol.getNumBonds()) {
                std::string linkage = monomer_mol.getBondWithIdx(connection_bond)
                                          ->getProp<std::string>(LINKAGE);
                monomer_mol.getBondWithIdx(connection_bond)
                    ->setProp(CUSTOM_BOND, linkage);
            } else {
                // temporary error handling
                throw std::runtime_error("Could not find connection bond");
            }
        }
    }
    for (auto bond : monomer_mol.bonds()) {
        if (getPolymerId(bond->getBeginAtom()) !=
            getPolymerId(bond->getEndAtom())) {
            std::string linkage = bond->getProp<std::string>(LINKAGE);
            bond->setProp(CUSTOM_BOND, linkage);
        }
    }
}

} // namespace RDKit
