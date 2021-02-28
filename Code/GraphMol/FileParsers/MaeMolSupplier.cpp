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
#include <GraphMol/FileParsers/FileParserUtils.h>

#include <boost/tokenizer.hpp>

#include <maeparser/MaeConstants.hpp>
#include <maeparser/Reader.hpp>

using namespace schrodinger;
using RDKit::MolInterchange::bolookup;

namespace RDKit {

namespace {

const std::string PDB_ATOM_NAME = "s_m_pdb_atom_name";
const std::string PDB_RESIDUE_NAME = "s_m_pdb_residue_name";
const std::string PDB_CHAIN_NAME = "s_m_chain_name";
const std::string PDB_INSERTION_CODE = "s_m_insertion_code";
const std::string PDB_RESIDUE_NUMBER = "i_m_residue_number";
const std::string PDB_OCCUPANCY = "r_m_pdb_occupancy";
const std::string PDB_TFACTOR = "r_m_pdb_tfactor";

class PDBInfo {
 public:
  PDBInfo(const mae::IndexedBlock &atom_block) {
    try {
      m_atom_name = atom_block.getStringProperty(PDB_ATOM_NAME);
    } catch (std::out_of_range &) {
    }

    try {
      m_residue_name = atom_block.getStringProperty(PDB_RESIDUE_NAME);
    } catch (std::out_of_range &) {
    }

    try {
      m_chain_id = atom_block.getStringProperty(PDB_CHAIN_NAME);
    } catch (std::out_of_range &) {
    }

    try {
      m_insertion_code = atom_block.getStringProperty(PDB_INSERTION_CODE);
    } catch (std::out_of_range &) {
    }

    try {
      m_resnum = atom_block.getIntProperty(PDB_RESIDUE_NUMBER);
    } catch (std::out_of_range &) {
    }

    try {
      m_occupancy = atom_block.getRealProperty(PDB_OCCUPANCY);
    } catch (std::out_of_range &) {
    }

    try {
      m_tempfac = atom_block.getRealProperty(PDB_TFACTOR);
    } catch (std::out_of_range &) {
    }
  }

  void addPDBData(Atom *atom, size_t atom_num) {
    if (!m_atom_name || !m_atom_name->isDefined(atom_num)) {
      return;  // Need a PDB atom name to populate info
    }
    AtomPDBResidueInfo *rd_info =
        new AtomPDBResidueInfo(m_atom_name->at(atom_num));

    atom->setMonomerInfo(rd_info);

    if (m_residue_name && m_residue_name->isDefined(atom_num)) {
      rd_info->setResidueName(m_residue_name->at(atom_num));
    }

    if (m_chain_id && m_chain_id->isDefined(atom_num)) {
      rd_info->setChainId(m_chain_id->at(atom_num));
    }

    if (m_insertion_code && m_insertion_code->isDefined(atom_num)) {
      rd_info->setInsertionCode(m_insertion_code->at(atom_num));
    }

    if (m_resnum && m_resnum->isDefined(atom_num)) {
      rd_info->setResidueNumber(m_resnum->at(atom_num));
    }

    if (m_occupancy && m_occupancy->isDefined(atom_num)) {
      rd_info->setOccupancy(m_occupancy->at(atom_num));
    }

    if (m_tempfac && m_tempfac->isDefined(atom_num)) {
      rd_info->setTempFactor(m_tempfac->at(atom_num));
    }
  }

 private:
  std::shared_ptr<mae::IndexedStringProperty> m_atom_name;
  std::shared_ptr<mae::IndexedStringProperty> m_residue_name;
  std::shared_ptr<mae::IndexedStringProperty> m_chain_id;
  std::shared_ptr<mae::IndexedStringProperty> m_insertion_code;

  std::shared_ptr<mae::IndexedIntProperty> m_resnum;

  std::shared_ptr<mae::IndexedRealProperty> m_occupancy;
  std::shared_ptr<mae::IndexedRealProperty> m_tempfac;
};

bool streamIsGoodOrExhausted(std::istream *stream) {
  PRECONDITION(stream, "bad stream");
  return stream->good() || (stream->eof() && stream->fail() && !stream->bad());
}

void parseChiralityLabel(RWMol &mol, const std::string &stereo_prop) {
  boost::char_separator<char> sep{"_"};
  boost::tokenizer<boost::char_separator<char>> tokenizer{stereo_prop, sep};

  auto tItr = tokenizer.begin();

  const int chiral_idx = FileParserUtils::toInt(*tItr) - 1;
  Atom *chiral_atom = mol.getAtomWithIdx(chiral_idx);
  CHECK_INVARIANT(chiral_atom != nullptr, "bad prop value");

  unsigned nSwaps = 2;
  const char rotation_direction = (++tItr)->back();
  switch (rotation_direction) {
    case 'R':  // R, ANR
      nSwaps = 0;
      break;
    case 'S':  // S, ANS
      nSwaps = 1;
      break;
    case '?':  // Undefined
      return;
    default:
      break;
  }
  CHECK_INVARIANT(nSwaps < 2, "bad prop value");

  INT_LIST bond_indexes;
  for (++tItr; tItr != tokenizer.end(); ++tItr) {
    const int nbr_idx = FileParserUtils::toInt(*tItr) - 1;
    const Bond *bnd = mol.getBondBetweenAtoms(chiral_idx, nbr_idx);
    CHECK_INVARIANT(bnd, "bad chiral bond");
    bond_indexes.push_back(bnd->getIdx());
  }
  CHECK_INVARIANT(bond_indexes.size() == chiral_atom->getDegree(),
                  "bad prop value");

  nSwaps += chiral_atom->getPerturbationOrder(bond_indexes);
  switch (nSwaps % 2) {
    case 0:
      chiral_atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      break;
    case 1:
      chiral_atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
      break;
  }
}

void parseStereoBondLabel(RWMol &mol, const std::string &stereo_prop) {
  boost::char_separator<char> sep{"_"};
  boost::tokenizer<boost::char_separator<char>> tokenizer{stereo_prop, sep};

  Bond::BondStereo type = Bond::STEREONONE;
  std::vector<int> atom_indexes;
  for (const auto &t : tokenizer) {
    if (t == "E") {
      type = Bond::STEREOTRANS;
    } else if (t == "Z") {
      type = Bond::STEREOCIS;
    } else {
      // Atom indexes are 0-based in RDKit, and 1-based in Mae.
      atom_indexes.push_back(FileParserUtils::toInt(t) - 1);
    }
  }
  CHECK_INVARIANT(type != Bond::STEREONONE, "bad prop value");

  // We currently don't support allenes or allene-likes
  if (atom_indexes.size() != 4) {
    return;
  }

  auto *bond = mol.getBondBetweenAtoms(atom_indexes[1], atom_indexes[2]);
  CHECK_INVARIANT(bond, "bad stereo bond");
  CHECK_INVARIANT(bond->getBondType() == Bond::DOUBLE, "bad stereo bond");

  bond->setStereoAtoms(atom_indexes[0], atom_indexes[3]);
  bond->setStereo(type);
}

//! Copy over the structure properties, including stereochemistry.
void set_mol_properties(RWMol &mol, const mae::Block &ct_block) {
  for (const auto &prop : ct_block.getProperties<std::string>()) {
    if (prop.first == mae::CT_TITLE) {
      mol.setProp(common_properties::_Name, prop.second);
    } else if (prop.first.find(mae::CT_CHIRALITY_PROP_PREFIX) == 0 ||
               prop.first.find(mae::CT_PSEUDOCHIRALITY_PROP_PREFIX) == 0) {
      parseChiralityLabel(mol, prop.second);
    } else if (prop.first.find(mae::CT_EZ_PROP_PREFIX) == 0) {
      parseStereoBondLabel(mol, prop.second);
    } else {
      mol.setProp(prop.first, prop.second);
    }
  }
  for (const auto &prop : ct_block.getProperties<double>()) {
    mol.setProp(prop.first, prop.second);
  }
  for (const auto &prop : ct_block.getProperties<int>()) {
    mol.setProp(prop.first, prop.second);
  }
  for (const auto &prop : ct_block.getProperties<mae::BoolProperty>()) {
    mol.setProp(prop.first, static_cast<bool>(prop.second));
  }
}

//! Set atom properties. Some of these have already been parsed to construct the
//! Atom object, and should be skipped. Also, atom properties may be undefined
//! for an atom, and should also be skipped for that atom.
void set_atom_properties(Atom &atom, const mae::IndexedBlock &atom_block,
                         size_t i) {
  for (const auto &prop : atom_block.getProperties<std::string>()) {
    if (prop.first == PDB_ATOM_NAME || prop.first == PDB_RESIDUE_NAME ||
        prop.first == PDB_CHAIN_NAME || prop.first == PDB_INSERTION_CODE) {
      // PDB information is parsed separately.
      continue;
    } else if (!prop.second->isDefined(i)) {
      continue;
    }

    atom.setProp(prop.first, prop.second->at(i));
  }

  for (const auto &prop : atom_block.getProperties<double>()) {
    if (prop.first == mae::ATOM_X_COORD || prop.first == mae::ATOM_Y_COORD ||
        prop.first == mae::ATOM_Z_COORD) {
      // Coordinates are used in defining a conformation, and should not be
      // set on the atom.
      continue;
    } else if (prop.first == PDB_OCCUPANCY || prop.first == PDB_TFACTOR) {
      // PDB information is parsed separately.
      continue;
    } else if (!prop.second->isDefined(i)) {
      continue;
    }

    atom.setProp(prop.first, prop.second->at(i));
  }
  for (const auto &prop : atom_block.getProperties<int>()) {
    if (prop.first == mae::ATOM_ATOMIC_NUM) {
      // Atomic number was already used in the creation of the atom
      continue;
    } else if (prop.first == PDB_RESIDUE_NUMBER) {
      // PDB information is parsed separately.
      continue;
    } else if (!prop.second->isDefined(i)) {
      continue;
    }

    if (prop.first == mae::ATOM_FORMAL_CHARGE) {
      // Formal charge has a specific setter
      atom.setFormalCharge(prop.second->at(i));
    } else {
      atom.setProp(prop.first, prop.second->at(i));
    }
  }
  for (const auto &prop : atom_block.getProperties<mae::BoolProperty>()) {
    if (!prop.second->isDefined(i)) {
      continue;
    }

    atom.setProp(prop.first, static_cast<bool>(prop.second->at(i)));
  }
}

void addAtoms(const mae::IndexedBlock &atom_block, RWMol &mol) {
  // All atoms are guaranteed to have these three field names:
  const auto atomic_numbers = atom_block.getIntProperty(mae::ATOM_ATOMIC_NUM);
  const auto xs = atom_block.getRealProperty(mae::ATOM_X_COORD);
  const auto ys = atom_block.getRealProperty(mae::ATOM_Y_COORD);
  const auto zs = atom_block.getRealProperty(mae::ATOM_Z_COORD);

  // atomic numbers, and x, y, and z coordinates
  const auto size = atomic_numbers->size();
  auto conf = new RDKit::Conformer(size);
  conf->set3D(true);
  conf->setId(0);

  PDBInfo pdb_info(atom_block);

  for (size_t i = 0; i < size; ++i) {
    Atom *atom = new Atom(atomic_numbers->at(i));
    mol.addAtom(atom, true, true);

    pdb_info.addPDBData(atom, i);
    set_atom_properties(*atom, atom_block, i);

    RDGeom::Point3D pos;
    pos.x = xs->at(i);
    pos.y = ys->at(i);
    pos.z = zs->at(i);
    conf->setAtomPos(i, pos);
  }
  mol.addConformer(conf, false);
}

void addBonds(const mae::IndexedBlock &bond_block, RWMol &mol) {
  // All bonds are guaranteed to have these three field names:
  const auto from_atoms = bond_block.getIntProperty(mae::BOND_ATOM_1);
  const auto to_atoms = bond_block.getIntProperty(mae::BOND_ATOM_2);
  const auto orders = bond_block.getIntProperty(mae::BOND_ORDER);
  const auto size = from_atoms->size();

  for (size_t i = 0; i < size; ++i) {
    // Maestro atoms are 1 indexed!
    const auto from_atom = from_atoms->at(i) - 1;
    const auto to_atom = to_atoms->at(i) - 1;
    const auto order = bolookup.find(orders->at(i))->second;
    if (from_atom > to_atom) {
      continue;  // Maestro files may double-list some bonds
    }

    auto bond = new Bond(order);
    bond->setOwningMol(mol);
    bond->setBeginAtomIdx(from_atom);
    bond->setEndAtomIdx(to_atom);
    mol.addBond(bond, true);
  }
}

void build_mol(RWMol &mol, mae::Block &structure_block, bool sanitize,
               bool removeHs) {
  const auto &atom_block = structure_block.getIndexedBlock(mae::ATOM_BLOCK);
  addAtoms(*atom_block, mol);

  const auto &bond_block = structure_block.getIndexedBlock(mae::BOND_BLOCK);
  addBonds(*bond_block, mol);

  // These properties need to be set last, as stereochemistry is defined here,
  // and it requires atoms and bonds to be available.
  set_mol_properties(mol, structure_block);

  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(mol, false, false);
    } else {
      MolOps::sanitizeMol(mol);
    }
  } else {
    // we need some properties for the chiral setup
    mol.updatePropertyCache(false);
  }

  // If there are 3D coordinates, try to read more chiralities from them, but do
  // not override the ones that were read from properties
  bool replaceExistingTags = false;
  if (mol.getNumConformers() && mol.getConformer().is3D()) {
    MolOps::assignChiralTypesFrom3D(mol, -1, replaceExistingTags);
  }

  // Find more stereo bonds, assign labels, but don't replace the existing ones
  MolOps::detectBondStereochemistry(mol, replaceExistingTags);
  MolOps::assignStereochemistry(mol, replaceExistingTags);
}

}  // namespace

MaeMolSupplier::MaeMolSupplier(std::shared_ptr<std::istream> inStream,
                               bool sanitize, bool removeHs) {
  PRECONDITION(inStream, "bad stream");
  dp_sInStream = inStream;
  dp_inStream = inStream.get();
  df_owner = true;
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new mae::Reader(dp_sInStream));
  CHECK_INVARIANT(streamIsGoodOrExhausted(dp_inStream), "bad instream");

  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    throw FileParseException(e.what());
  }
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

  d_reader.reset(new mae::Reader(dp_sInStream));
  CHECK_INVARIANT(streamIsGoodOrExhausted(dp_inStream), "bad instream");

  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    throw FileParseException(e.what());
  }
}

MaeMolSupplier::MaeMolSupplier(const std::string &fileName, bool sanitize,
                               bool removeHs) {
  df_owner = true;
  dp_inStream = openAndCheckStream(fileName);
  dp_sInStream.reset(dp_inStream);
  df_sanitize = sanitize;
  df_removeHs = removeHs;

  d_reader.reset(new mae::Reader(dp_sInStream));
  CHECK_INVARIANT(streamIsGoodOrExhausted(dp_inStream), "bad instream");

  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    throw FileParseException(e.what());
  }
}

void MaeMolSupplier::init() {}
void MaeMolSupplier::reset() {}

ROMol *MaeMolSupplier::next() {
  PRECONDITION(dp_sInStream != nullptr, "no stream");
  if (!d_stored_exc.empty()) {
    throw FileParseException(d_stored_exc);
  } else if (atEnd()) {
    throw FileParseException("All structures read from Maestro file");
  }

  auto mol = new RWMol;

  try {
    build_mol(*mol, *d_next_struct, df_sanitize, df_removeHs);
  } catch (...) {
    delete mol;
    moveToNextBlock();
    throw;
  }

  moveToNextBlock();

  return static_cast<ROMol *>(mol);
}

void MaeMolSupplier::moveToNextBlock() {
  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    d_stored_exc = e.what();
  }
}

bool MaeMolSupplier::atEnd() { return d_next_struct == nullptr; }
}  // namespace RDKit