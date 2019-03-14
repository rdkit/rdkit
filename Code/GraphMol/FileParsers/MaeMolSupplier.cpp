//
//  Copyright (C) 2018 Pat Lorton
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <fstream>
#include <map>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolInterchange/details.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <maeparser/Reader.hpp>

namespace RDKit {

using RDKit::MolInterchange::bolookup;

MaeMolSupplier::MaeMolSupplier(std::shared_ptr<std::istream> inStream,
                               bool sanitize, bool removeHs) {
  PRECONDITION(inStream, "bad stream");
  dp_sInStream = inStream;
  dp_inStream = inStream.get();
  df_owner = true;
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
  d_next_struct = d_reader->next("f_m_ct");
}

MaeMolSupplier::MaeMolSupplier(std::istream *inStream, bool takeOwnership,
                               bool sanitize, bool removeHs) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(takeOwnership, "takeOwnership is required for MaeMolSupplier");
  dp_inStream = inStream;
  dp_sInStream.reset(dp_inStream);
  df_owner = takeOwnership;  // always true
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
  d_next_struct = d_reader->next("f_m_ct");
}

MaeMolSupplier::MaeMolSupplier(const std::string &fileName, bool sanitize,
                               bool removeHs) {
  df_owner = true;
  auto *ifs = new std::ifstream(fileName.c_str(), std::ios_base::binary);
  if (!ifs || !(*ifs) || ifs->bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  dp_inStream = static_cast<std::istream *>(ifs);
  dp_sInStream.reset(dp_inStream);
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new schrodinger::mae::Reader(dp_sInStream));
  d_next_struct = d_reader->next("f_m_ct");
}

void MaeMolSupplier::init() {}
void MaeMolSupplier::reset() {}

ROMol *MaeMolSupplier::next() {
  if (d_next_struct == nullptr) {
    throw FileParseException("All structures read from Maestro file");
  }
  // Make sure even if later calls except, we're ready to read the next struct
  auto current_struct = d_next_struct;
  d_next_struct = d_reader->next("f_m_ct");

  auto mol = new RWMol();
  auto mol_title = current_struct->getStringProperty("s_m_title");
  mol->setProp(common_properties::_Name, mol_title);
  // Atom data is in the m_atom indexed block
  {
    const auto atom_data = current_struct->getIndexedBlock("m_atom");
    // All atoms are guaranteed to have these three field names:
    const auto atomic_numbers = atom_data->getIntProperty("i_m_atomic_number");
    const auto xs = atom_data->getRealProperty("r_m_x_coord");
    const auto ys = atom_data->getRealProperty("r_m_y_coord");
    const auto zs = atom_data->getRealProperty("r_m_z_coord");
    const auto size = atomic_numbers->size();
    std::shared_ptr<schrodinger::mae::IndexedIntProperty> atomic_charges;
    try {
      atomic_charges = atom_data->getIntProperty("i_m_formal_charge");
    } catch (std::out_of_range &) {
    }

    // atomic numbers, and x, y, and z coordinates
    auto conf = new RDKit::Conformer(size);
    conf->set3D(true);
    conf->setId(0);
    for (size_t i = 0; i < size; ++i) {
      Atom *atom = new Atom(atomic_numbers->at(i));
      mol->addAtom(atom, true, true);
      if (atomic_charges) {
        atom->setFormalCharge(atomic_charges->at(i));
      }

      RDGeom::Point3D pos;
      pos.x = xs->at(i);
      pos.y = ys->at(i);
      pos.z = zs->at(i);
      conf->setAtomPos(i, pos);
    }
    mol->addConformer(conf, false);
  }

  // Bond data is in the m_bond indexed block
  {
    const auto bond_data = current_struct->getIndexedBlock("m_bond");
    // All bonds are guaranteed to have these three field names:
    auto from_atoms = bond_data->getIntProperty("i_m_from");
    auto to_atoms = bond_data->getIntProperty("i_m_to");
    auto orders = bond_data->getIntProperty("i_m_order");
    const auto size = from_atoms->size();

    for (size_t i = 0; i < size; ++i) {
      // Maestro atoms are 1 indexed!
      const auto from_atom = from_atoms->at(i) - 1;
      const auto to_atom = to_atoms->at(i) - 1;
      const auto order = bolookup.find(orders->at(i))->second;
      if (from_atom > to_atom) continue;  // Maestro files double-list bonds

      auto bond = new Bond(order);
      bond->setOwningMol(mol);
      bond->setBeginAtomIdx(from_atom);
      bond->setEndAtomIdx(to_atom);
      mol->addBond(bond, true);
    }
  }

  if (df_sanitize) {
    if (df_removeHs) {
      MolOps::removeHs(*mol, false, false);
    } else {
      MolOps::sanitizeMol(*mol);
    }
  } else {
    // we need some properties for the chiral setup
    mol->updatePropertyCache(false);
  }

  /* Set tetrahedral chirality from 3D co-ordinates */
  MolOps::assignChiralTypesFrom3D(*mol);
  MolOps::detectBondStereochemistry(*mol);

  return (ROMol *)mol;
}

bool MaeMolSupplier::atEnd() { return d_next_struct == nullptr; }
}  // namespace RDKit
