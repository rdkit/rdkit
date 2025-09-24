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
  for(auto &v : dst.getPropList(true, false)) {
    if(v != RDKit::detail::computedPropName)
      dst.clearProp(v);
  }
}

static void copySelectedAtomsAndBonds(::RDKit::RWMol &extracted_mol,
                                      const RDKit::ROMol &reference_mol,
                                      SubsetInfo &selection_info,
				      const SubsetOptions &options) {
  auto &[selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
      selection_info;
  for (const auto &ref_atom : reference_mol.atoms()) {
    if (!selected_atoms.at(ref_atom->getIdx())) {
      continue;
    }

    std::unique_ptr<::RDKit::Atom> extracted_atom{
      options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy() };

    static constinit bool updateLabel = false;
    static constinit bool takeOwnership = true;
    atom_mapping[ref_atom->getIdx()] = extracted_mol.addAtom(
        extracted_atom.release(), updateLabel, takeOwnership);
  }

  for (const auto &ref_bond : reference_mol.bonds()) {
    if (!selected_bonds[ref_bond->getIdx()]) {
      continue;
    }

    std::unique_ptr<::RDKit::Bond> extracted_bond{
      options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy() };

    // Check the stereo atoms
    auto &atoms = extracted_bond->getStereoAtoms();
    if(atoms.size() == 2) {
      auto map1 = atom_mapping.find(atoms[0]);
      auto map2 = atom_mapping.find(atoms[1]);
      if(map1 != atom_mapping.end() && map2 != atom_mapping.end()) {
	atoms[0] = map1->second;
	atoms[1] = map2->second;
      } else {
	atoms.clear(); // We couldn't map the stereo atoms
      }
    }
    
    for(auto &atomidx : atoms) {
      auto map = atom_mapping.find(atomidx);
      if(map != atom_mapping.end())
	atomidx = map->second;
    }
    
    extracted_bond->setBeginAtomIdx(atom_mapping[ref_bond->getBeginAtomIdx()]);
    extracted_bond->setEndAtomIdx(atom_mapping[ref_bond->getEndAtomIdx()]);

    static constinit bool takeOwnership = true;
    auto num_bonds =
        extracted_mol.addBond(extracted_bond.release(), takeOwnership);
    bond_mapping[ref_bond->getIdx()] = num_bonds - 1;

    // we need to update rings now
    if(reference_mol.getRingInfo()->isInitialized()) {
      extracted_mol.getRingInfo()->reset();
    }
  }
}

[[nodiscard]] static bool is_selected_sgroup(
    const ::RDKit::SubstanceGroup &sgroup,
    const SubsetInfo &selection_info) {
  auto is_selected_component = [](auto &indices, auto &selection_test) {
    return indices.empty() ||
           std::all_of(indices.begin(), indices.end(), selection_test);
  };

  // clang-format off
    auto atom_test = [&](int idx) { return selection_info.selected_atoms[idx]; };
    auto bond_test = [&](int idx) { return selection_info.selected_bonds[idx]; };
    return is_selected_component(sgroup.getAtoms(), atom_test) &&
           is_selected_component(sgroup.getBonds(), bond_test) &&
           is_selected_component(sgroup.getParentAtoms(), atom_test);
  // clang-format on
}

static void copySelectedSubstanceGroups(
    ::RDKit::RWMol &extracted_mol, const RDKit::ROMol &reference_mol,
    const SubsetInfo &selection_info, const SubsetOptions &options) {
  auto update_indices = [](auto &sgroup, auto getter, auto setter,
                           auto &mapping) {
    auto indices = getter(sgroup);
    std::for_each(indices.begin(), indices.end(),
                  [&](auto &idx) { idx = mapping.at(idx); });
    setter(sgroup, std::move(indices));
  };

  const auto &[selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
      selection_info;
  for (const auto &sgroup : ::RDKit::getSubstanceGroups(reference_mol)) {
    if (!is_selected_sgroup(sgroup, selection_info)) {
      continue;
    }

    ::RDKit::SubstanceGroup extracted_sgroup(sgroup);
    extracted_sgroup.setOwningMol(&extracted_mol);

    update_indices(
        extracted_sgroup, std::mem_fn(&::RDKit::SubstanceGroup::getAtoms),
        std::mem_fn(&::RDKit::SubstanceGroup::setAtoms), atom_mapping);
    update_indices(
        extracted_sgroup, std::mem_fn(&::RDKit::SubstanceGroup::getParentAtoms),
        std::mem_fn(&::RDKit::SubstanceGroup::setParentAtoms), atom_mapping);
    update_indices(
        extracted_sgroup, std::mem_fn(&::RDKit::SubstanceGroup::getBonds),
        std::mem_fn(&::RDKit::SubstanceGroup::setBonds), bond_mapping);

    ::RDKit::addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
  }
}

static void copySelectedStereoGroups(::RDKit::RWMol &extracted_mol,
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
                                 selection_info.selected_atoms) &&
           is_selected_component(stereo_group.getBonds(),
                                 selection_info.selected_bonds);
  };

  std::vector<::RDKit::Atom *> extracted_atoms(extracted_mol.getNumAtoms());
  for (const auto &atom : extracted_mol.atoms()) {
    extracted_atoms[atom->getIdx()] = atom;
  }

  std::vector<::RDKit::Bond *> extracted_bonds(extracted_mol.getNumBonds());
  for (const auto &bond : extracted_mol.bonds()) {
    extracted_bonds[bond->getIdx()] = bond;
  }

  const auto &[selected_atoms, selected_bonds, atom_mapping, bond_mapping] =
      selection_info;
  std::vector<::RDKit::StereoGroup> extracted_stereo_groups;
  for (const auto &stereo_group : reference_mol.getStereoGroups()) {
    if (!is_selected_stereo_group(stereo_group)) {
      continue;
    }

    std::vector<::RDKit::Atom *> atoms;
    for (const auto &atom : stereo_group.getAtoms()) {
      auto mapping = atom_mapping.find(atom->getIdx());
      if(mapping != atom_mapping.end())
	atoms.push_back(extracted_atoms[mapping->second]);
    }

    std::vector<::RDKit::Bond *> bonds;
    for (const auto &bond : stereo_group.getBonds()) {
      auto mapping = bond_mapping.find(bond->getIdx());
      if(mapping != bond_mapping.end()) {
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

static void getSubsetInfo(SubsetInfo &selection_info,
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    const SubsetOptions &options) {
  const auto num_atoms = mol.getNumAtoms();
  const auto num_bonds = mol.getNumBonds();
  selection_info.selected_atoms = std::vector<bool>(num_atoms);
  selection_info.selected_bonds = std::vector<bool>(num_bonds);
  selection_info.atom_mapping .clear();
  selection_info.bond_mapping.clear();
  
  auto &[selected_atoms, selected_bonds, atom_mapping, bond_mapping] = selection_info;
  
  if (options.method == SubsetMethod::BONDS_BETWEEN_ATOMS) {
    for (const auto &atom_idx : path) {
      if(atom_idx < num_atoms) {
	selected_atoms[atom_idx] = true;
      }
    }
    for (const auto &bond : mol.bonds()) {
      if (selected_atoms.at(bond->getBeginAtomIdx()) &&
          selected_atoms.at(bond->getEndAtomIdx()) ){
	selected_bonds.at(bond->getIdx()) = true;
      }
    }
  } else if (options.method == SubsetMethod::BONDS) {
    for(const auto &bond_idx : path) {
      if(bond_idx < num_bonds) {
	selected_bonds.at(bond_idx) = true;
	const auto &bnd = mol.getBondWithIdx(bond_idx);
	selected_atoms.at(bnd->getBeginAtomIdx()) = true;
	selected_atoms.at(bnd->getEndAtomIdx()) = true;
      }
    }
  }
}

void copyCoords(RDKit::RWMol &copy, const RDKit::RWMol &mol, SubsetInfo &subset_info) {
    if (mol.getNumConformers()) {
    // copy coordinates over:
    for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
         ++confIt) {
      auto *conf = new Conformer(copy.getNumAtoms());
      conf->set3D((*confIt)->is3D());
      for (auto &mapping : subset_info.atom_mapping) { 
        conf->setAtomPos(mapping.second, (*confIt)->getAtomPos(mapping.first));
      }
      conf->setId((*confIt)->getId());
      copy.addConformer(conf, false);
    }
  }
}
  
std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, SubsetInfo & selection_info,
    const SubsetOptions &options) {
  auto extracted_mol = std::make_unique<::RDKit::RWMol>();
  copySelectedAtomsAndBonds(*extracted_mol, mol, selection_info, options);
  copySelectedSubstanceGroups(*extracted_mol, mol, selection_info, options);
  copySelectedStereoGroups(*extracted_mol, mol, selection_info);
  if(options.copyCoordinates)
    copyCoords(*extracted_mol, mol, selection_info);
  
  if (options.sanitize) {
    ::RDKit::MolOps::sanitizeMol(*extracted_mol);
  }
  
  if (options.clearComputedProps) {
    // this clears atom/bond and molecule computed props
    extracted_mol->clearComputedProps();
  } else {
    // this copies the mol computed props over, the atoms/bonds already have their's
    copyComputedProps(mol, *extracted_mol);
  }

  return extracted_mol;
}
  
}

std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol,
    const std::vector<unsigned int> &atoms,
    const std::vector<unsigned int> &bonds,
    const SubsetOptions &options, SubsetInfo *mappings) {
  SubsetInfo *selection_info;
  if(mappings)
    selection_info = mappings;
  else
    selection_info = new SubsetInfo;

  const auto natoms = mol.getNumAtoms();
  const auto nbonds = mol.getNumBonds();
  
  if(
     (atoms.size() == natoms && (
				 options.method == SubsetMethod::BONDS_BETWEEN_ATOMS ||
				 bonds.size() == nbonds) ||
      (bonds.size() == nbonds && options.method == SubsetMethod::BONDS ))
     ) {
    // optimization to copy the entire thing
    std::cerr << " ----- subset full copy ---- " << std::endl;
    return std::make_unique<RDKit::RWMol>(mol);
  }
  

  selection_info->selected_atoms = std::vector<bool>(natoms);
  selection_info->selected_bonds = std::vector<bool>(nbonds);
  selection_info->atom_mapping.clear();
  selection_info->bond_mapping.clear();
  for(auto v: atoms) {
    selection_info->selected_atoms.at(v) = true;
  }
  for(auto v: bonds) {
    selection_info->selected_bonds.at(v) = true;
  }

  
  auto res = copyMolSubset(mol, *selection_info, options);
  if(!mappings)
    delete selection_info;

  return res;
}
  
std::unique_ptr<RDKit::RWMol> copyMolSubset(
    const RDKit::ROMol &mol, const std::vector<unsigned int> &path,
    const SubsetOptions &options, SubsetInfo *mappings) {
  SubsetInfo *selection_info = mappings==nullptr ? new SubsetInfo : mappings;
  
  getSubsetInfo(*selection_info, mol, path, options);
  auto res = copyMolSubset(mol, *selection_info, options);
  if (!mappings)
    delete selection_info;
  return res;
}
  
}
