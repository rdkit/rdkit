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
#include "MonomerLibrary.h"

#include <GraphMol/QueryAtom.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/MolPickler.h>

#include <boost/functional/hash.hpp>

namespace RDKit
{

// Helper functions
namespace {

std::unique_ptr<Atom> makeMonomer(std::string_view name,
                                  std::string_view chain_id,
                                  std::string_view monomer_class,
                                  int residue_number, bool is_smiles)
{
    auto a = std::make_unique<::RDKit::Atom>();
    std::string n{name};
    // Allows monomer names show up in image renderings
    a->setProp(RDKit::common_properties::atomLabel, n);
    // Allows monomer names to be written to SMILES
    a->setProp(RDKit::common_properties::smilesSymbol, n);
    // Always start with BRANCH_MONOMER as false, will be set to
    // true if branch linkage is made to this monomer. Will be used
    // for rendering and monomer ordering purposes.
    a->setProp(BRANCH_MONOMER, false);
    // Provided monomer can be an unnamed monomer represented by a SMILES string,
    // indicated by is_smiles=True
    a->setProp(SMILES_MONOMER, is_smiles);

    // Give canonicalization to monomers based on name
    // TODO: Canonicalize on name and class?
    static boost::hash<std::string> hasher;
    a->setIsotope(hasher(n));
    a->setNoImplicit(true);

    // This is how the information about this monomer is represented
    auto* monomer_info = new ::RDKit::AtomMonomerInfo();
    monomer_info->setResidueNumber(residue_number);
    monomer_info->setResidueName(n);
    monomer_info->setChainId(std::string{chain_id});
    monomer_info->setMonomerClass(std::string{monomer_class});
    a->setMonomerInfo(monomer_info);
    return a;
}

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
{
    // in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos) {
        throw std::runtime_error("Invalid linkage format: " + linkage);
    }
    return {std::stoi(linkage.substr(1, dash - 1)),
            std::stoi(linkage.substr(dash + 2))};
}

} // anonymous namespace

// MonomerMol constructor/assignment implementations

MonomerMol::MonomerMol(const MonomerMol &other)
    : ROMol(other), d_library(other.d_library) {}

MonomerMol &MonomerMol::operator=(const MonomerMol &other)
{
  if (this != &other) {
    this->clear();
    numBonds = 0;
    initFromOther(other, false, -1);
    d_library = other.d_library;
  }
  return *this;
}

MonomerMol::MonomerMol(MonomerMol &&other) noexcept
    : ROMol(std::move(other)), d_library(std::move(other.d_library)) {}

MonomerMol &MonomerMol::operator=(MonomerMol &&other) noexcept
{
    ROMol::operator=(std::move(other));
    d_library = std::move(other.d_library);
    return *this;
}

MonomerMol::MonomerMol(const std::string &binStr)
{
    MolPickler::molFromPickle(binStr, *this);
}

// MonomerLibrary access implementations

MonomerLibrary& MonomerMol::getMonomerLibrary()
{
    if (d_library) {
        return *d_library;
    }
    return MonomerLibrary::getGlobalLibrary();
}

const MonomerLibrary& MonomerMol::getMonomerLibrary() const
{
    if (d_library) {
        return *d_library;
    }
    return MonomerLibrary::getGlobalLibrary();
}

void MonomerMol::setMonomerLibrary(std::shared_ptr<MonomerLibrary> library)
{
    d_library = std::move(library);
}

// MonomerMol member function implementations

void MonomerMol::addConnection(size_t monomer1, size_t monomer2,
                               const std::string& linkage)
{
    if (!isMonomer(getAtomWithIdx(monomer1)) || !isMonomer(getAtomWithIdx(monomer2))) {
        throw std::runtime_error("addConnection must be called with two monomers.");
    }

    // if bond already exists, extend linkage information using 'extra_linkage' property
    if (auto bond = getBondBetweenAtoms(monomer1, monomer2); bond != nullptr) {

        // If the linkage property isn't set, something went wrong
        std::string existing_linkage;
        if (!bond->getPropIfPresent(LINKAGE, existing_linkage)) {
            throw std::runtime_error(
                "No linkage property on bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // Make sure we're not recreating this same bond.
        if (existing_linkage.find(linkage) != std::string::npos) {
            throw std::runtime_error(
                "Can't duplicate " + linkage + " bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // For now, only set 'extra linkage' on linkages that already exist as a part
        // of the backbone
        if (existing_linkage != BACKBONE_LINKAGE) {
            throw std::runtime_error(
                "Only dative bonds can have multiple linkages between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // Since hydrogen bonding and covalent bonds are different types of bond
        // serializing them with the same bond type doesn't make too much sense.
        if (linkage == HYDROGEN_LINKAGE) {
            throw std::runtime_error(
                "Multiple bonds can't include hydrogen bond for bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }

        // For now, don't allow >2 additional linkages between the same
        // two monomers
        if (bond->hasProp(EXTRA_LINKAGE)) {
            throw std::runtime_error(
                "Multiple custom bonds not supported for bond between atom=" +
                std::to_string(monomer1) + " and atom=" + std::to_string(monomer2));
        }


        bond->setProp(EXTRA_LINKAGE, linkage);
    } else {
        auto create_bond = [&, this](unsigned int first_monomer,
                               unsigned int second_monomer,
                               ::RDKit::Bond::BondType bond_type,
                               const std::string& linkage_prop) {
            const auto new_total = this->addBond(first_monomer, second_monomer, bond_type);
            auto bond = this->getBondWithIdx(new_total - 1);
            bond->setProp(LINKAGE, linkage_prop);
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
                auto new_linkage = "R" + std::to_string(end_attchpt) + "-R" + std::to_string(begin_attchpt);
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

        if (linkage == BRANCH_LINKAGE) {
            // monomer2 is a branch monomer
            this->getAtomWithIdx(monomer2)->setProp(BRANCH_MONOMER, true);
        }
    }
}

void MonomerMol::addAtomMonomerConnection(size_t atom_idx, size_t monomer_idx,
                               const std::string& linkage, Bond::BondType bond_type)
{
    if (!isMonomer(this->getAtomWithIdx(monomer_idx)) || isMonomer(this->getAtomWithIdx(atom_idx))) {
        throw std::runtime_error(
            "addAtomMonomerConnection must be called with an atom and a monomer.");
    }

    const auto new_total = this->addBond(atom_idx, monomer_idx, bond_type);
    auto bond = this->getBondWithIdx(new_total - 1);
    bond->setProp(LINKAGE, linkage);
}

unsigned int MonomerMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
                            Bond::BondType bondType) {
  // if the atom indices are bad, the next two calls will catch that.
  auto beginAtom = getAtomWithIdx(atomIdx1);
  auto endAtom = getAtomWithIdx(atomIdx2);
  PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
  PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
               "bond already exists");

  auto *b = new Bond(bondType);
  b->setOwningMol(this);
  if (bondType == Bond::AROMATIC) {
    b->setIsAromatic(1);
    //
    // assume that aromatic bonds connect aromatic atoms
    //   This is relevant for file formats like MOL, where there
    //   is no such thing as an aromatic atom, but bonds can be
    //   marked aromatic.
    //
    beginAtom->setIsAromatic(1);
    endAtom->setIsAromatic(1);
  }
  auto [which, ok] = boost::add_edge(atomIdx1, atomIdx2, d_graph);
  d_graph[which] = b;
  ++numBonds;
  b->setIdx(numBonds - 1);
  b->setBeginAtomIdx(atomIdx1);
  b->setEndAtomIdx(atomIdx2);

  // the valence values on the begin and end atoms need to be updated:
  beginAtom->clearPropertyCache();
  endAtom->clearPropertyCache();

  // we're in a batch edit, and at least one of the bond ends is scheduled
  // for deletion, so mark the new bond for deletion too:
  if (dp_delAtoms &&
      ((atomIdx1 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx1)) ||
       (atomIdx2 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx2)))) {
    if (dp_delBonds->size() < numBonds) {
      dp_delBonds->resize(numBonds);
    }
    dp_delBonds->set(numBonds - 1);
  }

  return numBonds;
}

std::vector<std::string> MonomerMol::getPolymerIds() const
{
    std::vector<std::string> polymer_ids;
    for (auto atom : atoms()) {
        auto id = getPolymerId(atom);
        // in vector to preseve order of polymers
        if (std::find(polymer_ids.begin(), polymer_ids.end(), id) ==
            polymer_ids.end()) {
            polymer_ids.push_back(id);
        }
    }
    return polymer_ids;
}

size_t MonomerMol::addMonomer(std::string_view name, int residue_number, std::string_view monomer_class,
                              std::string_view chain_id, MonomerType monomer_type)
{
    auto monomer = makeMonomer(name, chain_id, monomer_class, residue_number,
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
    const auto monomer_class = last_monomer->getMonomerInfo()->getMonomerClass();
    return addMonomer(name, residue_number, monomer_class, chain_id, monomer_type);
}

Chain MonomerMol::getPolymer(std::string_view polymer_id) const
{
    std::vector<unsigned int> chain_atoms;
    for (auto atom : atoms()) {
        if (getPolymerId(atom) == polymer_id) {
            chain_atoms.push_back(atom->getIdx());
        }
    }
    std::sort(chain_atoms.begin(), chain_atoms.end(),
              [this](unsigned int a, unsigned int b) {
                  return getResidueNumber(getAtomWithIdx(a)) <
                         getResidueNumber(getAtomWithIdx(b));
              });
    std::vector<unsigned int> chain_bonds;
    for (auto bond : bonds()) {
        if (getPolymerId(bond->getBeginAtom()) == polymer_id &&
            getPolymerId(bond->getEndAtom()) == polymer_id) {
            chain_bonds.push_back(bond->getIdx());
        }
    }

    // Annotations stored as a COP substance group ("heavy chain", "light chain")
    std::string annotation{};
    for (const auto& sg : ::RDKit::getSubstanceGroups(*this)) {
        if ((sg.getProp<std::string>("TYPE") != "COP") ||
            !sg.hasProp(ANNOTATION) || !sg.hasProp("ID")) {
            continue;
        }
        if (sg.getProp<std::string>("ID") == polymer_id) {
            annotation = sg.getProp<std::string>(ANNOTATION);
            break;
        }
    }
    return {chain_atoms, chain_bonds, annotation};
}

} // namespace RDKit
