//
//  Copyright (C) 2018 Pat Lorton
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>
#include <iostream>
#include <fstream>

#include <boost/unordered_set.hpp>

#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/MolInterchange/details.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/FileParsers/MaestroProperties.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParserUtils.h>

#include <boost/tokenizer.hpp>

#include <maeparser/MaeConstants.hpp>
#include <maeparser/Reader.hpp>

using namespace schrodinger;
using namespace RDKit::FileParsers::schrodinger;
using RDKit::MolInterchange::bolookup;

namespace RDKit {

namespace v2 {
namespace FileParsers {
namespace {

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

std::string strip_prefix_from_mae_property(const std::string &propName) {
  const char *propNamePtr = propName.c_str();
  if (*propNamePtr == 'b' || *propNamePtr == 'i' || *propNamePtr == 'r' ||
      *propNamePtr == 's') {
    ++propNamePtr;
    if (strncmp(propNamePtr, "_rdk_", 5) == 0) {
      return propName.substr(6);
    } else if (strncmp(propNamePtr, "_rdkit_", 7) == 0) {
      return propName.substr(8);
    }
  }
  return propName;
}

bool is_ignored_property(const std::string &prop) {
  static const boost::unordered_set<std::string> ignored_properties = {
      MAE_ENHANCED_STEREO_STATUS,
      MAE_STEREO_STATUS,
  };

  return ignored_properties.find(prop) != ignored_properties.end();
}

//! Copy over the structure properties, including stereochemistry.
void set_mol_properties(RWMol &mol, const mae::Block &ct_block) {
  for (const auto &prop : ct_block.getProperties<std::string>()) {
    if (is_ignored_property(prop.first)) {
      continue;
    }

    if (prop.first == mae::CT_TITLE) {
      mol.setProp(common_properties::_Name, prop.second);
    } else if (prop.first.find(mae::CT_CHIRALITY_PROP_PREFIX) == 0 ||
               prop.first.find(mae::CT_PSEUDOCHIRALITY_PROP_PREFIX) == 0) {
      parseChiralityLabel(mol, prop.second);
    } else if (prop.first.find(mae::CT_EZ_PROP_PREFIX) == 0) {
      parseStereoBondLabel(mol, prop.second);
    } else {
      auto propName = strip_prefix_from_mae_property(prop.first);
      mol.setProp(propName, prop.second);
    }
  }
  for (const auto &prop : ct_block.getProperties<double>()) {
    if (is_ignored_property(prop.first)) {
      continue;
    }

    auto propName = strip_prefix_from_mae_property(prop.first);
    mol.setProp(propName, prop.second);
  }
  for (const auto &prop : ct_block.getProperties<int>()) {
    if (is_ignored_property(prop.first)) {
      continue;
    }

    auto propName = strip_prefix_from_mae_property(prop.first);
    mol.setProp(propName, prop.second);
  }
  for (const auto &prop : ct_block.getProperties<mae::BoolProperty>()) {
    if (is_ignored_property(prop.first)) {
      continue;
    }

    auto propName = strip_prefix_from_mae_property(prop.first);
    mol.setProp(propName, static_cast<bool>(prop.second));
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

    auto propName = strip_prefix_from_mae_property(prop.first);
    atom.setProp(propName, prop.second->at(i));
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

    auto propName = strip_prefix_from_mae_property(prop.first);
    atom.setProp(propName, prop.second->at(i));
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
    } else if (prop.first == MAE_RGROUP_LABEL) {
      // Schrodinger adopted RDKit's Group label property,
      // but with a "i_sd_" prefix instead of the usual "i_rdkit_"
      atom.setProp(common_properties::_MolFileRLabel, prop.second->at(i));
    } else {
      auto propName = strip_prefix_from_mae_property(prop.first);
      atom.setProp(propName, prop.second->at(i));
    }
  }
  for (const auto &prop : atom_block.getProperties<mae::BoolProperty>()) {
    if (!prop.second->isDefined(i)) {
      continue;
    }

    auto propName = strip_prefix_from_mae_property(prop.first);
    atom.setProp(propName, static_cast<bool>(prop.second->at(i)));
  }
}

void addAtoms(const mae::IndexedBlock &atom_block, RWMol &mol,
              std::vector<unsigned int> &atomsToRemove) {
  // All atoms are guaranteed to have these three field names:
  const auto atomicNumbers = atom_block.getIntProperty(mae::ATOM_ATOMIC_NUM);
  const auto xs = atom_block.getRealProperty(mae::ATOM_X_COORD);
  const auto ys = atom_block.getRealProperty(mae::ATOM_Y_COORD);
  const auto zs = atom_block.getRealProperty(mae::ATOM_Z_COORD);

  // atomic numbers, and x, y, and z coordinates
  const auto size = atomicNumbers->size();
  auto conf = new RDKit::Conformer(size);
  conf->setId(0);

  PDBInfo pdb_info(atom_block);

  bool nonzeroZ = false;
  for (size_t i = 0; i < size; ++i) {
    bool removeAtom = false;
    auto atomicNumber = atomicNumbers->at(i);
    if (atomicNumber == 0 || atomicNumber == -1 || atomicNumber == -3) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: atom " << (i + 1)
          << " in input Maestro file has atomic number '" << atomicNumber
          << "', which is reserved for internal use, and not allowed in inputs."
          << " The atom will be ignored.";

      // removing the atom now would be problematic for parsing the bonds
      // (especially if the atom is bonded!), so we'll just add a dummy
      // atom instead, and remove it again once we have finished parsing
      // the file.
      atomsToRemove.push_back(i);
      removeAtom = true;
      atomicNumber = 0;
    } else if (atomicNumber == -2) {
      // Maestro files use atomic number -2 to indicate a dummy atom.
      atomicNumber = 0;
    }

    Atom *atom = new Atom(atomicNumber);
    mol.addAtom(atom, true, true);

    RDGeom::Point3D pos;
    pos.x = xs->at(i);
    pos.y = ys->at(i);
    pos.z = zs->at(i);
    conf->setAtomPos(i, pos);

    // If the atom is going to be removed, don't bother with pdb info or
    // properties, and also don't consider it for planarity
    if (!removeAtom) {
      pdb_info.addPDBData(atom, i);
      set_atom_properties(*atom, atom_block, i);

      nonzeroZ |= (std::abs(pos.z) > 1.e-4);
    }
  }

  conf->set3D(nonzeroZ);
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
    if (auto bond = mol.getBondBetweenAtoms(from_atom, to_atom);
        bond != nullptr) {
      if (order != bond->getBondType()) {
        BOOST_LOG(rdWarningLog)
            << "WARNING: bond between atoms " << from_atom << " and " << to_atom
            << " is defined more than once with different bond orders. "
            << "The first definition will be honored, and the rest will be ignored.";
      }
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
  std::vector<unsigned int> atomsToRemove;
  const auto &atom_block = structure_block.getIndexedBlock(mae::ATOM_BLOCK);
  addAtoms(*atom_block, mol, atomsToRemove);

  std::shared_ptr<const mae::IndexedBlock> bond_block{nullptr};
  try {
    bond_block = structure_block.getIndexedBlock(mae::BOND_BLOCK);
  } catch (const std::out_of_range &) {
    // In Maestro files, the atom block is mandatory, but the bond block is not.
  }
  if (bond_block != nullptr) {
    addBonds(*bond_block, mol);
  }

  // These properties need to be set last, as stereochemistry is defined here,
  // and it requires atoms and bonds to be available.
  set_mol_properties(mol, structure_block);

  bool replaceExistingTags = false;
  if (sanitize) {
    if (removeHs) {
      // Bond stereo detection must happen before H removal, or
      // else we might be removing stereogenic H atoms in double
      // bonds (e.g. imines). But before we run stereo detection,
      // we need to run mol cleanup so don't have trouble with
      // e.g. nitro groups. Sadly, this a;; means we will find
      // run both cleanup and ring finding twice (a fast find
      // rings in bond stereo detection, and another in
      // sanitization's SSSR symmetrization).
      unsigned int failedOp = 0;
      MolOps::sanitizeMol(mol, failedOp, MolOps::SANITIZE_CLEANUP);
      MolOps::detectBondStereochemistry(mol);
      MolOps::removeHs(mol, false, false);
    } else {
      MolOps::sanitizeMol(mol);
      MolOps::detectBondStereochemistry(mol, replaceExistingTags);
    }
  } else {
    // we need some properties for the chiral setup
    mol.updatePropertyCache(false);
    MolOps::detectBondStereochemistry(mol, replaceExistingTags);
  }

  // If there are 3D coordinates, try to read more chiralities from them, but do
  // not override the ones that were read from properties

  if (mol.getNumConformers() && mol.getConformer().is3D()) {
    MolOps::assignChiralTypesFrom3D(mol, -1, replaceExistingTags);
  }

  // Assign labels, but don't replace the existing ones
  MolOps::assignStereochemistry(mol, replaceExistingTags);

  // If we saw any invalid atoms, remove them now
  if (!atomsToRemove.empty()) {
    mol.beginBatchEdit();
    for (auto aidx : atomsToRemove) {
      mol.removeAtom(aidx);
    }
    mol.commitBatchEdit();
  }
}

void throw_idx_error(unsigned idx) {
  std::ostringstream errout;
  errout << "ERROR: Index error (idx = " << idx << ") : "
         << " we do no have enough ct blocks";
  throw FileParseException(errout.str());
}

}  // namespace

MaeMolSupplier::MaeMolSupplier(std::shared_ptr<std::istream> inStream,
                               const MaeMolSupplierParams &params) {
  PRECONDITION(inStream, "bad stream");
  dp_sInStream = inStream;
  dp_inStream = inStream.get();
  df_owner = true;
  d_params = params;

  init();
}

MaeMolSupplier::MaeMolSupplier(std::istream *inStream, bool takeOwnership,
                               const MaeMolSupplierParams &params) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(takeOwnership, "takeOwnership is required for MaeMolSupplier");
  dp_inStream = inStream;
  dp_sInStream.reset(dp_inStream);
  df_owner = takeOwnership;  // always true
  d_params = params;

  init();
}

MaeMolSupplier::MaeMolSupplier(const std::string &fileName,
                               const MaeMolSupplierParams &params) {
  df_owner = true;
  dp_inStream = openAndCheckStream(fileName);
  dp_sInStream.reset(dp_inStream);
  d_params = params;

  init();
}

void MaeMolSupplier::init() {
  PRECONDITION(dp_sInStream, "no input stream")
  d_reader.reset(new mae::Reader(dp_sInStream));
  CHECK_INVARIANT(streamIsGoodOrExhausted(dp_inStream), "bad instream");

  d_position = 0;
  d_length = 0;

  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    throw FileParseException(e.what());
  }
}
void MaeMolSupplier::reset() {
  dp_inStream->clear();
  dp_inStream->seekg(0, std::ios::beg);

  auto length = d_length;
  init();
  d_length = length;
}

void MaeMolSupplier::setData(const std::string &text,
                             const MaeMolSupplierParams &params) {
  dp_inStream = static_cast<std::istream *>(
      new std::istringstream(text, std::ios_base::binary));
  dp_sInStream.reset(dp_inStream);
  df_owner = true;  // maeparser requires ownership
  d_params = params;
  init();
}

std::unique_ptr<RWMol> MaeMolSupplier::next() {
  PRECONDITION(dp_sInStream != nullptr, "no stream");
  if (!d_stored_exc.empty()) {
    throw FileParseException(d_stored_exc);
  } else if (atEnd()) {
    throw FileParseException("All structures read from Maestro file");
  }

  auto mol = std::make_unique<RWMol>();

  try {
    build_mol(*mol, *d_next_struct, d_params.sanitize, d_params.removeHs);
  } catch (const std::exception &e) {
    moveToNextBlock();
    throw FileParseException(e.what());
  }

  moveToNextBlock();

  return mol;
}

void MaeMolSupplier::moveToNextBlock() {
  try {
    d_next_struct = d_reader->next(mae::CT_BLOCK);
  } catch (const mae::read_exception &e) {
    d_stored_exc = e.what();
  }
  ++d_position;
}

bool MaeMolSupplier::atEnd() {
  if (d_next_struct == nullptr) {
    d_length = d_position;
    return true;
  }
  return false;
}

unsigned int MaeMolSupplier::length() {
  PRECONDITION(dp_inStream, "no stream");

  if (d_length == 0 && !atEnd()) {
    // maeparser has an internal buffer, so we can't just iterate over
    // block till we reach the end of the file. So we have to rewind
    // the input stream, use it to create a separate parser, fast
    // forward this one to the end of the data, and then get the length
    // from that parser. Then we can restore the input stream to
    // the position where it was before, so that it is still in
    // sync with maeparser's internal buffer.

    dp_sInStream->clear();
    auto current_position = dp_sInStream->tellg();
    dp_sInStream->seekg(0, std::ios::beg);

    MaeMolSupplier tmp_supplier(dp_sInStream);
    while (!tmp_supplier.atEnd()) {
      tmp_supplier.moveToNextBlock();
    }

    d_length = tmp_supplier.length();
    dp_sInStream->seekg(current_position, std::ios::beg);
  }

  return d_length;
}

void MaeMolSupplier::moveTo(unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");

  if (d_length > 0 && idx > d_length) {
    throw_idx_error(idx);
  }

  if (idx < d_position) {
    reset();
  }

  while (idx > d_position) {
    moveToNextBlock();

    if (atEnd()) {
      throw_idx_error(idx);
    }
  }
}

std::unique_ptr<RWMol> MaeMolSupplier::operator[](unsigned int idx) {
  PRECONDITION(dp_inStream, "no stream");
  moveTo(idx);
  return next();
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
