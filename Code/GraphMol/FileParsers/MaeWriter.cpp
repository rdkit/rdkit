//
//
//  Copyright (C) 2023 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolWriters.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <maeparser/MaeBlock.hpp>
#include <maeparser/MaeConstants.hpp>
#include <maeparser/Writer.hpp>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>

using namespace schrodinger;

namespace RDKit {

namespace {

template <typename T>
std::shared_ptr<mae::IndexedProperty<T>> getIndexedProperty(
    mae::IndexedBlock& indexedBlock, const std::string& propName,
    size_t numAtoms,
    std::shared_ptr<mae::IndexedProperty<T>> (mae::IndexedBlock::*getterFunc)(
        const std::string&) const,
    void (mae::IndexedBlock::*setterFunc)(
        const std::string&, std::shared_ptr<mae::IndexedProperty<T>>)) {
  auto prop = (indexedBlock.*getterFunc)(propName);
  if (prop != nullptr) {
    return prop;
  }

  std::vector<T> boolValues(numAtoms);
  auto nullsBitset = new boost::dynamic_bitset<>(numAtoms);
  nullsBitset->set();
  auto newProp =
      std::make_shared<mae::IndexedProperty<T>>(boolValues, nullsBitset);

  (indexedBlock.*setterFunc)(propName, newProp);
  return newProp;
}

template <typename T>
std::shared_ptr<mae::IndexedProperty<T>> getIndexedProperty(
    mae::IndexedBlock& indexedBlock, const std::string& propName,
    size_t numAtoms);

template <>
std::shared_ptr<mae::IndexedProperty<mae::BoolProperty>>
getIndexedProperty<mae::BoolProperty>(mae::IndexedBlock& indexedBlock,
                                      const std::string& propName,
                                      size_t numAtoms) {
  if (propName[0] != 'b') {
    auto msg =
        std::string("Property '") + propName + "' is not a boolean value";
    throw std::runtime_error(msg);
  }

  return getIndexedProperty<mae::BoolProperty>(
      indexedBlock, propName, numAtoms, &mae::IndexedBlock::getBoolProperty,
      &mae::IndexedBlock::setBoolProperty);
}

template <>
std::shared_ptr<mae::IndexedProperty<int>> getIndexedProperty<int>(
    mae::IndexedBlock& indexedBlock, const std::string& propName,
    size_t numAtoms) {
  if (propName[0] != 'i') {
    auto msg =
        std::string("Property '") + propName + "' is not an integer value";
    throw std::runtime_error(msg);
  }

  return getIndexedProperty<int>(indexedBlock, propName, numAtoms,
                                 &mae::IndexedBlock::getIntProperty,
                                 &mae::IndexedBlock::setIntProperty);
}

template <>
std::shared_ptr<mae::IndexedProperty<double>> getIndexedProperty<double>(
    mae::IndexedBlock& indexedBlock, const std::string& propName,
    size_t numAtoms) {
  if (propName[0] != 'r') {
    auto msg = std::string("Property '") + propName + "' is not a real value";
    throw std::runtime_error(msg);
  }

  return getIndexedProperty<double>(indexedBlock, propName, numAtoms,
                                    &mae::IndexedBlock::getRealProperty,
                                    &mae::IndexedBlock::setRealProperty);
}

template <>
std::shared_ptr<mae::IndexedProperty<std::string>>
getIndexedProperty<std::string>(mae::IndexedBlock& indexedBlock,
                                const std::string& propName, size_t numAtoms) {
  if (propName[0] != 's') {
    auto msg = std::string("Property '") + propName + "' is not a string value";
    throw std::runtime_error(msg);
  }

  return getIndexedProperty<std::string>(indexedBlock, propName, numAtoms,
                                         &mae::IndexedBlock::getStringProperty,
                                         &mae::IndexedBlock::setStringProperty);
}

class UnsupportedBondException : public std::runtime_error {
 public:
  UnsupportedBondException(const char* msg) : std::runtime_error(msg) {}
  UnsupportedBondException(const std::string& msg) : std::runtime_error(msg) {}
};

unsigned bondTypeToOrder(const Bond& bond) {
  switch (bond.getBondType()) {
    case Bond::BondType::SINGLE:
      return 1;
    case Bond::BondType::DOUBLE:
      return 2;
    case Bond::BondType::TRIPLE:
      return 3;
    case Bond::BondType::ZERO:
    case Bond::BondType::DATIVE:
      return 0;
    default:
      throw UnsupportedBondException(
          "Bond " + std::to_string(bond.getIdx()) +
          " has a type that is not supported by maeparser.");
  }
}

void mapMolProperties(const ROMol& mol, const STR_VECT& propNames,
                      mae::Block& stBlock) {
  (void)mol;
  (void)propNames;
  (void)stBlock;

#if 0
  for (auto& prop : propNames) {
    if (!mol.hasProp(prop)) {
      // We don't raise any warnings if we don't find a property because it
      // might be an atom property
      continue;
    }

    auto prop = mol.getProp(prop)
  }

  for (const auto& propName : st.getPropertyNames(PropertyType::Bool)) {
    stBlock.setBoolProperty(propName, st.getProperty<bool>(propName));
  }
  for (const auto& propName : st.getPropertyNames(PropertyType::Int)) {
    stBlock.setIntProperty(propName, st.getProperty<int>(propName));
  }
  for (const auto& propName : st.getPropertyNames(PropertyType::Real)) {
    stBlock.setRealProperty(propName, st.getProperty<double>(propName));
  }
  for (const auto& propName : st.getPropertyNames(PropertyType::String)) {
    stBlock.setStringProperty(propName,
                               st.getProperty<std::string>(propName));
  }
#endif
}

void mapAtom(const Conformer& conformer, const Atom& atom,
             const STR_VECT& propNames, mae::IndexedBlock& atomBlock,
             size_t numAtoms) {
  auto idx = atom.getIdx();
  auto coordinates = conformer.getAtomPos(idx);

  auto xCoord =
      getIndexedProperty<double>(atomBlock, mae::ATOM_X_COORD, numAtoms);
  xCoord->set(idx, coordinates.x);

  auto yCoord =
      getIndexedProperty<double>(atomBlock, mae::ATOM_Y_COORD, numAtoms);
  yCoord->set(idx, coordinates.y);

  auto zCoord =
      getIndexedProperty<double>(atomBlock, mae::ATOM_Z_COORD, numAtoms);
  zCoord->set(idx, coordinates.z);

  auto atomicNum =
      getIndexedProperty<int>(atomBlock, mae::ATOM_ATOMIC_NUM, numAtoms);
  atomicNum->set(idx, atom.getAtomicNum());

  // Default heavy atoms to be drawn in grey in Maestro, and H atoms in white.
  std::string color = (atom.getAtomicNum() == 1 ? "FFFFFF" : "A0A0A0");
  auto atomRgbColor =
      getIndexedProperty<std::string>(atomBlock, "s_m_color_rgb", numAtoms);
  atomRgbColor->set(idx, color);

  auto formalCharge =
      getIndexedProperty<int>(atomBlock, mae::ATOM_FORMAL_CHARGE, numAtoms);
  formalCharge->set(idx, atom.getFormalCharge());

#if 0
  auto residue_id = atom.getResidueId();
  auto residue_num =
      getIndexedProperty<int>(atomBlock, M2IO_DATA_RES_NUM, numAtoms);
  residue_num->set(idx, residue_id.residue_number);
  auto insertion_code = getIndexedProperty<std::string>(
      atomBlock, M2IO_DATA_INSERTION_CODE, numAtoms);
  insertion_code->set(idx, std::string{residue_id.insertion_code});

  auto mmod_residue = getIndexedProperty<std::string>(
      atomBlock, M2IO_DATA_MMOD_RES, numAtoms);
  mmod_residue->set(idx, std::string{atom.getMacromodelResidue()});

  auto chain_name =
      getIndexedProperty<std::string>(atomBlock, M2IO_DATA_CHAIN, numAtoms);
  chain_name->set(idx, atom.getChain());


  auto pdb_residue = getIndexedProperty<std::string>(
      atomBlock, M2IO_DATA_PDB_RES, numAtoms);
  pdb_residue->set(idx, atom.getPDBResidue());

  auto pdb_atom = getIndexedProperty<std::string>(
      atomBlock, M2IO_DATA_PDB_ATOM, numAtoms);
  pdb_atom->set(idx, atom.getPDBAtomName());



  // Custom properties
  for (const auto& propName : atom.getPropertyNames(PropertyType::Bool)) {
    auto idxd_prop = getIndexedProperty<mae::BoolProperty>(
        atomBlock, propName, numAtoms);
    idxd_prop->set(idx, atom.getProperty<bool>(propName));
  }
  for (const auto& propName : atom.getPropertyNames(PropertyType::Int)) {
    // M2IO_DATA_ATOM_CHIRALITY should not be written to .mae files
    // (see mmct_ct_m2io_close_atomBlock in mmct.cpp)
    if (propName == M2IO_DATA_ATOM_CHIRALITY) {
      continue;
    }

    auto idxd_prop =
        getIndexedProperty<int>(atomBlock, propName, numAtoms);
    idxd_prop->set(idx, atom.getProperty<int>(propName));
  }
  for (const auto& propName : atom.getPropertyNames(PropertyType::Real)) {
    auto idxd_prop =
        getIndexedProperty<double>(atomBlock, propName, numAtoms);
    idxd_prop->set(idx, atom.getProperty<double>(propName));
  }
  for (const auto& propName : atom.getPropertyNames(PropertyType::String)) {
    auto idxd_prop =
        getIndexedProperty<std::string>(atomBlock, propName, numAtoms);
    idxd_prop->set(idx, atom.getProperty<std::string>(propName));
  }

#endif
}

void mapAtoms(const ROMol& mol, const STR_VECT& propNames, int confId,
              mae::IndexedBlockMap& indexedBlockMap) {
  auto atomBlock = std::make_shared<mae::IndexedBlock>(mae::ATOM_BLOCK);
  auto conformer = mol.getConformer(confId);

  auto numAtoms = mol.getNumAtoms();
  for (auto& atom : mol.atoms()) {
    mapAtom(conformer, *atom, propNames, *atomBlock, numAtoms);
  }

  indexedBlockMap.addIndexedBlock(mae::ATOM_BLOCK, atomBlock);
}

void mapBonds(const ROMol& mol, const STR_VECT& propNames,
              mae::IndexedBlockMap& indexedBlockMap) {
  auto bondBlock = std::make_shared<mae::IndexedBlock>(mae::BOND_BLOCK);

  auto numBonds = mol.getNumBonds();
  auto bondAtomFrom =
      getIndexedProperty<int>(*bondBlock, mae::BOND_ATOM_1, numBonds);
  auto bondAtomTo =
      getIndexedProperty<int>(*bondBlock, mae::BOND_ATOM_2, numBonds);
  auto bondOrder =
      getIndexedProperty<int>(*bondBlock, mae::BOND_ORDER, numBonds);

  std::shared_ptr<mae::IndexedProperty<mae::BoolProperty>> dativeBondMark =
      nullptr;
  for (auto& bond : mol.bonds()) {
    if (bond->getBondType() == Bond::BondType::DATIVE) {
      dativeBondMark = getIndexedProperty<mae::BoolProperty>(
          *bondBlock, "b_sPrivate_dative_bond", numBonds);
      break;
    }
  }

  for (auto& bond : mol.bonds()) {
    auto idx = bond->getIdx();
    bondAtomFrom->set(idx, bond->getBeginAtomIdx());
    bondAtomTo->set(idx, bond->getEndAtomIdx());
    bondOrder->set(idx, bondTypeToOrder(*bond));

    if (dativeBondMark != nullptr) {
      dativeBondMark->set(idx, (bond->getBondType() == Bond::BondType::DATIVE));
    }
  }

  indexedBlockMap.addIndexedBlock(mae::BOND_BLOCK, bondBlock);
}
}  // namespace

MaeWriter::MaeWriter(const std::string& fileName) {
  auto* tmpStream = new std::ofstream(fileName.c_str());
  if (!(*tmpStream) || (tmpStream->bad())) {
    delete tmpStream;
    std::ostringstream errout;
    errout << "Bad output file " << fileName;
    throw BadFileException(errout.str());
  }
  dp_ostream.reset(static_cast<std::ostream*>(tmpStream));
}

MaeWriter::MaeWriter(std::ostream* outStream) : dp_ostream{outStream} {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
}

MaeWriter::MaeWriter(std::shared_ptr<std::ostream> outStream)
    : dp_ostream{std::move(outStream)} {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
}

MaeWriter::~MaeWriter() { close(); };

void MaeWriter::open() { dp_writer.reset(new mae::Writer(dp_ostream)); }

void MaeWriter::setProps(const STR_VECT& propNames) {
  if (d_molid > 0) {
    BOOST_LOG(rdWarningLog) << "WARNING: Setting property list after a few "
                               "molecules have been written\n";
  }
  d_props = propNames;
}

void MaeWriter::flush() {
  PRECONDITION(dp_ostream, "no output stream");
  try {
    dp_ostream->flush();
  } catch (...) {
    try {
      if (dp_ostream->good()) {
        dp_ostream->setstate(std::ios::badbit);
      }
    } catch (const std::runtime_error&) {
    }
  }
}

void MaeWriter::close() {
  if (dp_writer) {
    dp_writer.reset();
  }
  if (dp_ostream) {
    flush();
  }
  dp_ostream.reset();
}

void MaeWriter::write(const ROMol& mol, int confId) {
  PRECONDITION(dp_ostream, "no output stream");

  RWMol tmpMol(mol);
  if (!MolOps::KekulizeIfPossible(tmpMol)) {
    BOOST_LOG(rdErrorLog)
        << "ERROR: the mol cannot be kekulized, and will not be written to the output file.\n";
    return;
  }
  if (mol.getNumConformers() == 0) {
    // make sure there's at least one conformer we can write
    RDDepict::compute2DCoords(tmpMol);
  }

  auto stBlock = std::make_shared<mae::Block>(mae::CT_BLOCK);

  mapMolProperties(tmpMol, d_props, *stBlock);

  auto indexedBlockMap = std::make_shared<mae::IndexedBlockMap>();

  mapAtoms(tmpMol, d_props, confId, *indexedBlockMap);

  try {
    mapBonds(tmpMol, d_props, *indexedBlockMap);
  } catch (const UnsupportedBondException& exc) {
    BOOST_LOG(rdErrorLog)
        << "ERROR: " << exc.what()
        << " The mol will not be written to the output file.\n";
    return;
  }

  stBlock->setIndexedBlockMap(indexedBlockMap);

  if (!dp_writer) {
    open();
  }

  dp_writer->write(stBlock);
  ++d_molid;
}

}  // namespace RDKit