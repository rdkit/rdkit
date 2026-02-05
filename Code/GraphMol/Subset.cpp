//
//  Copyright (C) 2025 Hussein Faara, Brian Kelley and other RDKit Contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include "Subset.h"

namespace RDKit {
namespace {

inline void copyComputedProps(const ROMol &src, ROMol &dst) {
  dst.updateProps(src);
  for (auto &v : dst.getPropList(true, false)) {
    if (v != RDKit::detail::computedPropName) dst.clearProp(v);
  }
}

static void copySelectedAtomsAndBonds(RWMol &extracted_mol,
                                      const RDKit::ROMol &reference_mol,
                                      SubsetInfo &selection_info,
                                      const SubsetOptions &options) {
  auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
      selection_info;
  for (const auto &ref_atom : reference_mol.atoms()) {
    if (!selectedAtoms[ref_atom->getIdx()]) {
      continue;
    }

    std::unique_ptr<Atom> extracted_atom{
        options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
    extracted_atom->clearComputedProps();

    constexpr bool updateLabel = false;
    constexpr bool takeOwnership = true;
    atomMapping[ref_atom->getIdx()] = extracted_mol.addAtom(
        extracted_atom.release(), updateLabel, takeOwnership);
  }

  for (const auto &ref_bond : reference_mol.bonds()) {
    if (!selectedBonds[ref_bond->getIdx()]) {
      continue;
    }
    if (atomMapping.find(ref_bond->getBeginAtomIdx()) == atomMapping.end() ||
	atomMapping.find(ref_bond->getEndAtomIdx()) == atomMapping.end()) {
      throw ValueErrorException("copyMolSubset: subset bonds contain atoms not contained in subset atoms");
    }
	
    std::unique_ptr<Bond> extracted_bond{
        options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};

    // Check the stereo atoms
    auto &atoms = extracted_bond->getStereoAtoms();
    if (atoms.size() == 2) {
      auto map1 = atomMapping.find(atoms[0]);
      auto map2 = atomMapping.find(atoms[1]);
      if (map1 != atomMapping.end() && map2 != atomMapping.end()) {
        atoms[0] = map1->second;
        atoms[1] = map2->second;
      } else {
        atoms.clear();  // We couldn't map the stereo atoms
      }
    }

    for (auto &atomidx : atoms) {
      auto map = atomMapping.find(atomidx);
      if (map != atomMapping.end()) {
        atomidx = map->second;
      }
    }

    extracted_bond->setBeginAtomIdx(atomMapping[ref_bond->getBeginAtomIdx()]);
    extracted_bond->setEndAtomIdx(atomMapping[ref_bond->getEndAtomIdx()]);

    constexpr bool takeOwnership = true;
    auto num_bonds =
        extracted_mol.addBond(extracted_bond.release(), takeOwnership);
    bondMapping[ref_bond->getIdx()] = num_bonds - 1;
  }
  // we need to update rings now

  if (selectedBonds.any() && reference_mol.getRingInfo()->isInitialized()) {
    extracted_mol.getRingInfo()->reset();
  }
}

static bool isSelectedSGroup(const SubstanceGroup &sgroup,
                             const SubsetInfo &selection_info) {
  auto is_selected_component = [](auto &indices, auto &selection_test) {
    return indices.empty() ||
           std::all_of(indices.begin(), indices.end(), selection_test);
  };

  auto atom_test = [&](int idx) { return selection_info.selectedAtoms[idx]; };
  auto bond_test = [&](int idx) { return selection_info.selectedBonds[idx]; };
  return is_selected_component(sgroup.getAtoms(), atom_test) &&
         is_selected_component(sgroup.getBonds(), bond_test) &&
         is_selected_component(sgroup.getParentAtoms(), atom_test);
}

static void copySelectedSubstanceGroups(RWMol &extracted_mol,
                                        const RDKit::ROMol &reference_mol,
                                        const SubsetInfo &selection_info,
                                        const SubsetOptions &) {
  auto update_indices = [](auto &sgroup, auto getter, auto setter,
                           auto &mapping) {
    auto indices = getter(sgroup);
    std::for_each(indices.begin(), indices.end(),
                  [&](auto &idx) { idx = mapping.at(idx); });
    setter(sgroup, std::move(indices));
  };

  const auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
      selection_info;
  for (const auto &sgroup : getSubstanceGroups(reference_mol)) {
    if (!isSelectedSGroup(sgroup, selection_info)) {
      continue;
    }

    SubstanceGroup extracted_sgroup(sgroup);
    extracted_sgroup.setOwningMol(&extracted_mol);

    update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getAtoms),
                   std::mem_fn(&SubstanceGroup::setAtoms), atomMapping);
    update_indices(extracted_sgroup,
                   std::mem_fn(&SubstanceGroup::getParentAtoms),
                   std::mem_fn(&SubstanceGroup::setParentAtoms), atomMapping);
    update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getBonds),
                   std::mem_fn(&SubstanceGroup::setBonds), bondMapping);

    addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
  }
}

static void copySelectedStereoGroups(RWMol &extracted_mol,
                                     const RDKit::ROMol &reference_mol,
                                     const SubsetInfo &selection_info) {
  auto is_selected_component = [](auto &objects, auto &selected_indices) {
    return objects.empty() ||
           std::any_of(objects.begin(), objects.end(), [&](auto &object) {
             return selected_indices[object->getIdx()];
           });
  };

  auto is_selected_stereo_group = [&](const auto &stereo_group) {
    return is_selected_component(stereo_group.getAtoms(),
                                 selection_info.selectedAtoms) &&
           is_selected_component(stereo_group.getBonds(),
                                 selection_info.selectedBonds);
  };

  std::vector<Atom *> extracted_atoms(extracted_mol.getNumAtoms());
  for (const auto &atom : extracted_mol.atoms()) {
    extracted_atoms[atom->getIdx()] = atom;
  }

  std::vector<Bond *> extracted_bonds(extracted_mol.getNumBonds());
  for (const auto &bond : extracted_mol.bonds()) {
    extracted_bonds[bond->getIdx()] = bond;
  }

  const auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
      selection_info;
  std::vector<StereoGroup> extracted_stereo_groups;
  for (const auto &stereo_group : reference_mol.getStereoGroups()) {
    if (!is_selected_stereo_group(stereo_group)) {
      continue;
    }

    std::vector<Atom *> atoms;
    for (const auto &atom : stereo_group.getAtoms()) {
      auto mapping = atomMapping.find(atom->getIdx());
      if (mapping != atomMapping.end()) {
        atoms.push_back(extracted_atoms[mapping->second]);
      }
    }

    std::vector<Bond *> bonds;
    for (const auto &bond : stereo_group.getBonds()) {
      auto mapping = bondMapping.find(bond->getIdx());
      if (mapping != bondMapping.end()) {
        bonds.push_back(extracted_bonds[mapping->second]);
      }
    }

    extracted_stereo_groups.push_back({stereo_group.getGroupType(),
                                       std::move(atoms), std::move(bonds),
                                       stereo_group.getReadId()});
    extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
  }

  extracted_mol.setStereoGroups(std::move(extracted_stereo_groups));
}

static void getSubsetInfo(SubsetInfo &selection_info, const RDKit::ROMol &mol,
                          const std::vector<unsigned int> &path,
                          const SubsetOptions &options) {
  const auto num_atoms = mol.getNumAtoms();
  const auto num_bonds = mol.getNumBonds();
  selection_info.selectedAtoms.clear();
  selection_info.selectedAtoms.resize(num_atoms);
  selection_info.selectedBonds.clear();
  selection_info.selectedBonds.resize(num_bonds);

  selection_info.atomMapping.clear();
  selection_info.bondMapping.clear();

  auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
      selection_info;

  if (options.method == SubsetMethod::BONDS_BETWEEN_ATOMS) {
    for (const auto &atom_idx : path) {
      if (atom_idx < num_atoms) {
        selectedAtoms.set(atom_idx);
      }
    }
    for (const auto &bond : mol.bonds()) {
      if (selectedAtoms[bond->getBeginAtomIdx()] &&
          selectedAtoms[bond->getEndAtomIdx()]) {
        selectedBonds.set(bond->getIdx());
      }
    }
  } else if (options.method == SubsetMethod::BONDS) {
    for (const auto &bond_idx : path) {
      if (bond_idx < num_bonds) {
        selectedBonds.set(bond_idx);
        const auto &bnd = mol.getBondWithIdx(bond_idx);
        selectedAtoms.set(bnd->getBeginAtomIdx());
        selectedAtoms.set(bnd->getEndAtomIdx());
      }
    }
  }
}

void copyCoords(RDKit::RWMol &copy, const RDKit::RWMol &mol,
                SubsetInfo &subset_info) {
  if (mol.getNumConformers()) {
    // copy coordinates over:
    for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
         ++confIt) {
      auto *conf = new Conformer(copy.getNumAtoms());
      conf->set3D((*confIt)->is3D());
      for (auto &mapping : subset_info.atomMapping) {
        conf->setAtomPos(mapping.second, (*confIt)->getAtomPos(mapping.first));
      }
      conf->setId((*confIt)->getId());
      copy.addConformer(conf, false);
    }
  }
}

std::unique_ptr<RDKit::RWMol> copyMolSubset(const RDKit::ROMol &mol,
                                            SubsetInfo &selection_info,
                                            const SubsetOptions &options) {
  auto extracted_mol = std::make_unique<RWMol>();
  copySelectedAtomsAndBonds(*extracted_mol, mol, selection_info, options);
  copySelectedSubstanceGroups(*extracted_mol, mol, selection_info, options);
  copySelectedStereoGroups(*extracted_mol, mol, selection_info);
  if (options.copyCoordinates) {
    copyCoords(*extracted_mol, mol, selection_info);
  }

  if (options.sanitize) {
    MolOps::sanitizeMol(*extracted_mol);
  }

  if (options.clearComputedProps) {
    // this clears atom/bond and molecule computed props
    extracted_mol->clearComputedProps();
  } else {
    // this copies the mol computed props over, the atoms/bonds already have
    // their's
    copyComputedProps(mol, *extracted_mol);
  }

  return extracted_mol;
}

}  // namespace

std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &atoms,
    const std::vector<unsigned int> &bonds, const SubsetOptions &options) {
  SubsetInfo mappings;
  return copyMolSubset(mol, atoms, bonds, mappings, options);
}

std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &atoms,
    const std::vector<unsigned int> &bonds, SubsetInfo &selection_info,
    const SubsetOptions &options) {
  const auto natoms = mol.getNumAtoms();
  const auto nbonds = mol.getNumBonds();
  if ((atoms.size() == natoms && bonds.size() == nbonds)) {
    // optimization to copy the entire thing
    return std::make_unique<RDKit::RWMol>(mol);
  }

  selection_info.selectedAtoms.reset();
  selection_info.selectedAtoms.resize(natoms);
  selection_info.selectedBonds.reset();
  selection_info.selectedBonds.resize(nbonds);
  selection_info.atomMapping.clear();
  selection_info.bondMapping.clear();

  for (auto v : atoms) {
    if (v < natoms) {
      selection_info.selectedAtoms.set(v);
    } else {
      throw IndexErrorException(static_cast<int>(v));
    }
  }
  for (auto v : bonds) {
    if (v < nbonds) {
      selection_info.selectedBonds.set(v);
    } else {
      throw IndexErrorException(static_cast<int>(v));
    }
  }

  auto res = copyMolSubset(mol, selection_info, options);
  return res;
}

std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    const SubsetOptions &options) {
  SubsetInfo selection_info;
  return copyMolSubset(mol, path, selection_info, options);
}

std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    SubsetInfo &selection_info, const SubsetOptions &options) {
  getSubsetInfo(selection_info, mol, path, options);
  auto res = copyMolSubset(mol, selection_info, options);
  return res;
}

}  // namespace RDKit
