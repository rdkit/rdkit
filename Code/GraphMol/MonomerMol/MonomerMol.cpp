/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: tools for Coarse-grain ROMols
 *
 * A "coarse grain" (CG) ROMol uses RDKit atoms to represent monomers. Chains
 * are represented by COP Substance Groups on the ROMol.
 *
 * For use with functionality in schrodinger::rdkit_extensions
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "MonomerMol.h"

#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MonomerInfo.h>

#include <boost/functional/hash.hpp>

template <typename T>
static bool is_dummy_obj(const T& obj)
{
    // placeholder from HELM parser; represents a query/repeated monomer or bond
    return false;
}

namespace RDKit
{

std::string get_polymer_id(const Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getChainId();
}

int get_residue_number(const Atom* atom)
{
    const auto monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom does not have monomer info");
    }

    return static_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info)
        ->getResidueNumber();
}

std::vector<std::string> get_polymer_ids(const RDKit::ROMol& monomer_mol)
{
    std::vector<std::string> polymer_ids;
    for (auto atom : monomer_mol.atoms()) {
        if (is_dummy_obj(*atom)) { // query/repeated monomer
            continue;
        }
        auto id = get_polymer_id(atom);
        // in vector to preseve order of polymers
        if (std::find(polymer_ids.begin(), polymer_ids.end(), id) ==
            polymer_ids.end()) {
            polymer_ids.push_back(id);
        }
    }
    return polymer_ids;
}

Chain get_polymer(const RDKit::ROMol& monomer_mol,
                                       std::string_view polymer_id)
{
    std::vector<unsigned int> atoms;
    for (auto atom : monomer_mol.atoms()) {
        if (is_dummy_obj(*atom)) {
            continue;
        }
        if (get_polymer_id(atom) == polymer_id) {
            atoms.push_back(atom->getIdx());
        }
    }
    // Sort by get_residue_num
    std::sort(atoms.begin(), atoms.end(),
              [&monomer_mol](unsigned int a, unsigned int b) {
                  return get_residue_number(monomer_mol.getAtomWithIdx(a)) <
                         get_residue_number(monomer_mol.getAtomWithIdx(b));
              });
    std::vector<unsigned int> bonds;
    for (auto bond : monomer_mol.bonds()) {
        if (is_dummy_obj(*bond)) {
            continue;
        }
        if (get_polymer_id(bond->getBeginAtom()) == polymer_id &&
            get_polymer_id(bond->getEndAtom()) == polymer_id) {
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

ChainType to_chain_type(std::string_view chain_type)
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

std::string to_string(ChainType chain_type)
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

void add_connection(RDKit::RWMol& monomer_mol, size_t monomer1, size_t monomer2,
                    const std::string& linkage)
{
    const auto new_total =
        monomer_mol.addBond(monomer1, monomer2, ::RDKit::Bond::BondType::SINGLE);
    monomer_mol.getBondWithIdx(new_total - 1)->setProp(LINKAGE, linkage);
    if (linkage == BRANCH_LINKAGE) {
        // monomer2 is a branch monomer
        monomer_mol.getAtomWithIdx(monomer2)->setProp(BRANCH_MONOMER, true);
    }
}

void add_connection(RDKit::RWMol& monomer_mol, size_t monomer1, size_t monomer2,
                    ConnectionType connection_type)
{
    switch (connection_type) {
        case ConnectionType::FORWARD:
            add_connection(monomer_mol, monomer1, monomer2, BACKBONE_LINKAGE);
            break;
        case ConnectionType::SIDECHAIN:
            add_connection(monomer_mol, monomer1, monomer2, BRANCH_LINKAGE);
            break;
    }
}

std::unique_ptr<Monomer> make_monomer(std::string_view name,
                                      std::string_view chain_id,
                                      int residue_number, bool is_smiles)
{
    auto a = std::make_unique<::RDKit::Atom>();
    std::string n{name};
    a->setProp(ATOM_LABEL, n); // temporary, for image generation
    // Always start with BRANCH_MONOMER as false, will be set to
    // true if branch linkage is made to this monomer
    a->setProp(BRANCH_MONOMER, false);
    a->setProp(SMILES_MONOMER, is_smiles);

    // hack to get some level of canonicalization for MonomerMols
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

size_t add_monomer(RDKit::RWMol& monomer_mol, std::string_view name,
                   int residue_number, std::string_view chain_id,
                   MonomerType monomer_type)
{
    auto monomer = make_monomer(name, chain_id, residue_number,
                                monomer_type == MonomerType::SMILES);
    bool update_label = true;
    bool take_ownership = true;
    auto new_index =
        monomer_mol.addAtom(monomer.release(), update_label, take_ownership);
    return new_index;
}

size_t add_monomer(RDKit::RWMol& monomer_mol, std::string_view name,
                   MonomerType monomer_type)
{
    if (monomer_mol.getNumAtoms() == 0) {
        throw std::invalid_argument(
            "No atoms in molecule to determine chain ID");
    }
    const auto* last_monomer = monomer_mol.getAtomWithIdx(monomer_mol.getNumAtoms() - 1);
    const auto chain_id = get_polymer_id(last_monomer);
    const auto residue_number = get_residue_number(last_monomer) + 1;
    return add_monomer(monomer_mol, name, residue_number, chain_id, monomer_type);
}

void set_residue_number(RDKit::Atom* atom, int residue_number)
{
    auto* residue_info =
        static_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    if (residue_info == nullptr) {
        throw std::runtime_error("some Atom does not have residue info");
    }
    residue_info->setResidueNumber(residue_number);
}

bool is_valid_chain(const RDKit::RWMol& monomer_mol, std::string_view polymer_id)
{
    // Check that the residue ordering is valid for this polymer. The residues
    // should be in connectivity order
    const auto chain = get_polymer(monomer_mol, polymer_id);
    auto last_residue = chain.atoms[0];
    for (size_t i = 1; i < chain.atoms.size(); ++i) {
        const auto bond =
            monomer_mol.getBondBetweenAtoms(last_residue, chain.atoms[i]);
        if (bond == nullptr) {
            return false;
        }

        // Bond direction should be in the same order as residues
        if (get_residue_number(bond->getBeginAtom()) >
                get_residue_number(bond->getEndAtom()) &&
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

RDKit::Atom* find_chain_begin(RDKit::RWMol& monomer_mol)
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

void order_residues(RDKit::RWMol& monomer_mol)
{
    // Currently assumes that all monomers are in the same chain. We will
    // eventually want to order residues on a per-chain basis.

    // Find the beginning of the chain (must be a backbone monomer)
    auto chain_begin = find_chain_begin(monomer_mol);

    // Now re-order the residues beginning at chain_begin
    std::vector<RDKit::Atom*> queue;
    std::vector<bool> visited(monomer_mol.getNumAtoms(), 0);
    queue.push_back(chain_begin);
    visited[chain_begin->getIdx()] = true;

    int current_res_idx = 1;
    while (!queue.empty()) {
        auto current = queue.back();
        queue.pop_back();
        set_residue_number(current, current_res_idx);
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
                set_residue_number(other, current_res_idx);
                ++current_res_idx;
            } else if (linkage == BACKBONE_LINKAGE) {
                queue.push_back(bond->getOtherAtom(current));
            }
        }
    }
}

} // namespace RDKit