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
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <maeparser/Reader.hpp>

using std::shared_ptr;
using namespace schrodinger::mae;

namespace RDKit {

using RDKit::MolInterchange::bolookup;

MaeMolSupplier::MaeMolSupplier(shared_ptr<std::istream> inStream,
                               bool sanitize, bool removeHs) {
  PRECONDITION(inStream, "bad stream");
  dp_sInStream = inStream;
  dp_inStream = inStream.get();
  df_owner = true;
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new Reader(dp_sInStream));
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

  d_reader.reset(new Reader(dp_sInStream));
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

  d_reader.reset(new Reader(dp_sInStream));
  d_next_struct = d_reader->next("f_m_ct");
}

void MaeMolSupplier::init() {}
void MaeMolSupplier::reset() {}

class PDBInfo
{
public:
    PDBInfo(const shared_ptr<const IndexedBlock>& atom_data);

    shared_ptr<IndexedStringProperty> m_atom_name;
    shared_ptr<IndexedStringProperty> m_residue_name;
    shared_ptr<IndexedStringProperty> m_chain_id;
    shared_ptr<IndexedStringProperty> m_insertion_code;

    shared_ptr<IndexedIntProperty> m_resnum;

    shared_ptr<IndexedRealProperty> m_occupancy;
    shared_ptr<IndexedRealProperty> m_tempfac;

};

PDBInfo::PDBInfo(const shared_ptr<const IndexedBlock>& atom_data)
{
    try {
    m_atom_name = atom_data->getStringProperty("s_m_pdb_atom_name");
    } catch (std::out_of_range &) { }

    try {
    m_residue_name = atom_data->getStringProperty("s_m_pdb_residue_name");
    } catch (std::out_of_range &) { }

    try {
    m_chain_id = atom_data->getStringProperty("s_m_chain_name");
    } catch (std::out_of_range &) { }

    try {
    m_insertion_code = atom_data->getStringProperty("s_m_insertion_code");
    } catch (std::out_of_range &) { }

    try {
    m_resnum = atom_data->getIntProperty("i_m_residue_number");
    } catch (std::out_of_range &) { }

    try {
    m_occupancy = atom_data->getRealProperty("r_m_pdb_occupancy");
    } catch (std::out_of_range &) { }

    try {
    m_tempfac = atom_data->getRealProperty("r_m_pdb_tfactor");
    } catch (std::out_of_range &) { }
}

void addAtomPDBData(const PDBInfo& pdb_info, Atom* atom, size_t atom_num)
{
    if(!pdb_info.m_atom_name->isDefined(atom_num)) {
        return;  // Need a PDB atom name to populate info
    }
    AtomPDBResidueInfo *rd_info = new AtomPDBResidueInfo(
            pdb_info.m_atom_name->at(atom_num));

    atom->setMonomerInfo(rd_info);

    if(pdb_info.m_residue_name && pdb_info.m_residue_name->isDefined(atom_num))
    {
        rd_info->setResidueName(pdb_info.m_residue_name->at(atom_num));
    }

    if(pdb_info.m_chain_id && pdb_info.m_chain_id->isDefined(atom_num))
    {
        rd_info->setChainId(pdb_info.m_chain_id->at(atom_num));
    }

    if(pdb_info.m_insertion_code &&
            pdb_info.m_insertion_code->isDefined(atom_num))
    {
        rd_info->setInsertionCode(pdb_info.m_insertion_code->at(atom_num));
    }

    if(pdb_info.m_resnum && pdb_info.m_resnum->isDefined(atom_num))
    {
        rd_info->setResidueNumber(pdb_info.m_resnum->at(atom_num));
    }

    if(pdb_info.m_occupancy && pdb_info.m_occupancy->isDefined(atom_num))
    {
        rd_info->setOccupancy(pdb_info.m_occupancy->at(atom_num));
    }

    if(pdb_info.m_tempfac && pdb_info.m_tempfac->isDefined(atom_num))
    {
        rd_info->setTempFactor(pdb_info.m_tempfac->at(atom_num));
    }
}

void addAtoms(shared_ptr<Block>& current_struct, RWMol& mol)
{
  // Atom data is in the m_atom indexed block
  {
    const auto atom_data = current_struct->getIndexedBlock("m_atom");
    // All atoms are guaranteed to have these three field names:
    const auto atomic_numbers = atom_data->getIntProperty("i_m_atomic_number");
    const auto xs = atom_data->getRealProperty("r_m_x_coord");
    const auto ys = atom_data->getRealProperty("r_m_y_coord");
    const auto zs = atom_data->getRealProperty("r_m_z_coord");
    const auto atom_count = atomic_numbers->size();
    shared_ptr<IndexedIntProperty> atomic_charges;
    try {
      atomic_charges = atom_data->getIntProperty("i_m_formal_charge");
    } catch (std::out_of_range &) {
    }

    // atomic numbers, and x, y, and z coordinates
    auto conf = new RDKit::Conformer(atom_count);
    conf->set3D(true);
    conf->setId(0);
    PDBInfo pdb_info(atom_data);
    for (size_t i = 0; i < atom_count; ++i) {
      Atom *atom = new Atom(atomic_numbers->at(i));
      mol.addAtom(atom, true, true);
      if (atomic_charges) {
        atom->setFormalCharge(atomic_charges->at(i));
      }

      RDGeom::Point3D pos;
      pos.x = xs->at(i);
      pos.y = ys->at(i);
      pos.z = zs->at(i);
      conf->setAtomPos(i, pos);

      if(pdb_info.m_atom_name) {
          addAtomPDBData(pdb_info, atom, i);
      }
    }
    mol.addConformer(conf, false);
  }

}

void addBonds(shared_ptr<Block>& current_struct, RWMol& mol)
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
        mol.addBond(bond, true);
    }
}


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

  addAtoms(current_struct, *mol);
  addBonds(current_struct, *mol);

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
