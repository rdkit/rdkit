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
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <maeparser/MaeBlock.hpp>
#include <maeparser/MaeConstants.hpp>
#include <maeparser/Writer.hpp>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>

using namespace schrodinger;

namespace RDKit {

namespace {
const std::string MAE_ATOM_RGB_COLOR = "s_m_color_rgb";
const std::string MAE_BOND_DATIVE_MARK = "b_sPrivate_dative_bond";
const std::string PDB_ATOM_NAME = "s_m_pdb_atom_name";
const std::string PDB_RESIDUE_NAME = "s_m_pdb_residue_name";
const std::string PDB_CHAIN_NAME = "s_m_chain_name";
const std::string PDB_INSERTION_CODE = "s_m_insertion_code";
const std::string PDB_RESIDUE_NUMBER = "i_m_residue_number";
const std::string PDB_OCCUPANCY = "r_m_pdb_occupancy";
const std::string PDB_TFACTOR = "r_m_pdb_tfactor";

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

  std::vector<T> values(numAtoms);
  auto nullsBitset = new boost::dynamic_bitset<>(numAtoms);
  nullsBitset->set();
  auto newProp = std::make_shared<mae::IndexedProperty<T>>(values, nullsBitset);

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

void copyProperties(
    const RDProps& origin, const STR_VECT& propNames, unsigned idx,
    std::function<void(const std::string&, unsigned, bool)> boolSetter,
    std::function<void(const std::string&, unsigned, int)> intSetter,
    std::function<void(const std::string&, unsigned, double)> realSetter,
    std::function<void(const std::string&, unsigned, const std::string&)>
        stringSetter) {
  // Map other properties, but first clear out the computed ones,
  // since we don't want to export these.
  origin.clearComputedProps();

  for (const auto& prop : origin.getDict().getData()) {
    // Skip the property holding the names of the computed properties
    if (prop.key == detail::computedPropName) {
      continue;
    }

    // Also skip the property if we have a list of properties we want to export
    // and this one is not one of them.
    if (!propNames.empty() && (std::find(propNames.begin(), propNames.end(),
                                         prop.key) == propNames.end())) {
      continue;
    }

    switch (prop.val.getTag()) {
      case RDTypeTag::BoolTag: {
        auto propName = std::string("b_rdk_") + prop.key;
        boolSetter(propName, idx, rdvalue_cast<bool>(prop.val));
        break;
      }

      case RDTypeTag::IntTag:
      case RDTypeTag::UnsignedIntTag: {
        auto propName = std::string("i_rdk_") + prop.key;
        intSetter(propName, idx, rdvalue_cast<int>(prop.val));
        break;
      }

      case RDTypeTag::DoubleTag:
      case RDTypeTag::FloatTag: {
        auto propName = std::string("r_rdk_") + prop.key;
        realSetter(propName, idx, rdvalue_cast<double>(prop.val));
        break;
      }

      case RDTypeTag::StringTag: {
        auto propName = std::string("s_rdk_") + prop.key;
        stringSetter(propName, idx, rdvalue_cast<std::string>(prop.val));
        break;
      }
      default:
        // AnyTag, EmptyTag, any kind of Vector.
        // Should we support vectors?
        BOOST_LOG(rdWarningLog)
            << "WARNING: the property " << prop.key
            << " has an unsupported type, and will not be exported.";
        continue;
    }
  }
}

template <typename T>
void setPropertyValue(mae::IndexedBlock& block, const std::string& propName,
                      size_t numValues, unsigned idx, const T& value) {
  auto property = getIndexedProperty<T>(block, propName, numValues);
  property->set(idx, value);
}

class UnsupportedBondException : public std::runtime_error {
 public:
  UnsupportedBondException(const char* msg) : std::runtime_error(msg) {}
  UnsupportedBondException(const std::string& msg) : std::runtime_error(msg) {}
};

int bondTypeToOrder(const Bond& bond) {
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
  // We always write a title, even if the mol doesn't have one
  // (in such case, we add an empty string).
  std::string molName;
  mol.getPropIfPresent(common_properties::_Name, molName);
  stBlock.setStringProperty(mae::CT_TITLE, molName);
  mol.clearProp(common_properties::_Name);

  // TO DO: Map stereo properties

  auto boolSetter = [&stBlock](const std::string& prop, unsigned, bool value) {
    stBlock.setBoolProperty(prop, value);
  };
  auto intSetter = [&stBlock](const std::string& prop, unsigned, int value) {
    stBlock.setIntProperty(prop, value);
  };
  auto realSetter = [&stBlock](const std::string& prop, unsigned,
                               double value) {
    stBlock.setRealProperty(prop, value);
  };
  auto stringSetter = [&stBlock](const std::string& prop, unsigned,
                                 const std::string& value) {
    stBlock.setStringProperty(prop, value);
  };

  int fake_idx = 0;
  copyProperties(mol, propNames, fake_idx, boolSetter, intSetter, realSetter,
                 stringSetter);
}
void mapAtom(
    const Conformer& conformer, const Atom& atom, const STR_VECT& propNames,
    const std::string& heavyAtomColor, mae::IndexedBlock& atomBlock,
    size_t numAtoms,
    std::function<void(const std::string&, unsigned, bool)> boolSetter,
    std::function<void(const std::string&, unsigned, int)> intSetter,
    std::function<void(const std::string&, unsigned, double)> realSetter,
    std::function<void(const std::string&, unsigned, const std::string&)>
        stringSetter) {
  auto idx = atom.getIdx();
  auto coordinates = conformer.getAtomPos(idx);

  // Required properties
  setPropertyValue(atomBlock, mae::ATOM_X_COORD, numAtoms, idx, coordinates.x);
  setPropertyValue(atomBlock, mae::ATOM_Y_COORD, numAtoms, idx, coordinates.y);
  setPropertyValue(atomBlock, mae::ATOM_Z_COORD, numAtoms, idx, coordinates.z);

  auto atomic_num = static_cast<int>(atom.getAtomicNum());
  setPropertyValue(atomBlock, mae::ATOM_ATOMIC_NUM, numAtoms, idx, atomic_num);

  // Default heavy atoms to be drawn in grey in Maestro, and H atoms in white.
  std::string color = (atom.getAtomicNum() == 1 ? "FFFFFF" : heavyAtomColor);
  setPropertyValue(atomBlock, MAE_ATOM_RGB_COLOR, numAtoms, idx, color);

  setPropertyValue(atomBlock, mae::ATOM_FORMAL_CHARGE, numAtoms, idx,
                   atom.getFormalCharge());

  // Residue information
  auto monomerInfo =
      static_cast<const AtomPDBResidueInfo*>(atom.getMonomerInfo());
  if (monomerInfo != nullptr) {
    setPropertyValue(atomBlock, PDB_ATOM_NAME, numAtoms, idx,
                     monomerInfo->getName());

    setPropertyValue(atomBlock, PDB_RESIDUE_NAME, numAtoms, idx,
                     monomerInfo->getResidueName());

    setPropertyValue(atomBlock, PDB_CHAIN_NAME, numAtoms, idx,
                     monomerInfo->getChainId());

    setPropertyValue(atomBlock, PDB_INSERTION_CODE, numAtoms, idx,
                     monomerInfo->getInsertionCode());

    setPropertyValue(atomBlock, PDB_RESIDUE_NUMBER, numAtoms, idx,
                     monomerInfo->getResidueNumber());

    setPropertyValue(atomBlock, PDB_OCCUPANCY, numAtoms, idx,
                     monomerInfo->getOccupancy());

    setPropertyValue(atomBlock, PDB_TFACTOR, numAtoms, idx,
                     monomerInfo->getTempFactor());
  }

  // Custom properties
  copyProperties(atom, propNames, idx, boolSetter, intSetter, realSetter,
                 stringSetter);
}

void mapAtoms(const ROMol& mol, const STR_VECT& propNames,
              const std::string& heavyAtomColor, int confId,
              mae::IndexedBlockMap& indexedBlockMap) {
  auto atomBlock = std::make_shared<mae::IndexedBlock>(mae::ATOM_BLOCK);
  auto conformer = mol.getConformer(confId);
  auto numAtoms = mol.getNumAtoms();

  auto boolSetter = [&atomBlock, &numAtoms](const std::string& prop,
                                            unsigned idx,
                                            mae::BoolProperty value) {
    setPropertyValue(*atomBlock, prop, numAtoms, idx, value);
  };
  auto intSetter = [&atomBlock, &numAtoms](const std::string& prop,
                                           unsigned idx, int value) {
    setPropertyValue(*atomBlock, prop, numAtoms, idx, value);
  };
  auto realSetter = [&atomBlock, &numAtoms](const std::string& prop,
                                            unsigned idx, double value) {
    setPropertyValue(*atomBlock, prop, numAtoms, idx, value);
  };
  auto stringSetter = [&atomBlock, &numAtoms](const std::string& prop,
                                              unsigned idx,
                                              const std::string& value) {
    setPropertyValue(*atomBlock, prop, numAtoms, idx, value);
  };

  for (auto atom : mol.atoms()) {
    mapAtom(conformer, *atom, propNames, heavyAtomColor, *atomBlock, numAtoms,
            boolSetter, intSetter, realSetter, stringSetter);
  }

  indexedBlockMap.addIndexedBlock(mae::ATOM_BLOCK, atomBlock);
}

void mapBond(
    const Bond& bond,
    std::shared_ptr<mae::IndexedProperty<mae::BoolProperty>>& dativeBondMark,
    const STR_VECT& propNames, mae::IndexedBlock& bondBlock, size_t numBonds,
    std::function<void(const std::string&, unsigned, bool)> boolSetter,
    std::function<void(const std::string&, unsigned, int)> intSetter,
    std::function<void(const std::string&, unsigned, double)> realSetter,
    std::function<void(const std::string&, unsigned, const std::string&)>
        stringSetter) {
  auto idx = bond.getIdx();

  // Indexes in the atom block are 1-based
  auto bondTo = static_cast<int>(bond.getBeginAtomIdx()) + 1;
  auto bondFrom = static_cast<int>(bond.getEndAtomIdx()) + 1;

  // There is no bond directionality in Maestro, and atom indexes
  // in bonds are usually written in ascending order
  if (bondFrom > bondTo) {
    std::swap(bondFrom, bondTo);
  }

  setPropertyValue(bondBlock, mae::BOND_ATOM_1, numBonds, idx, bondTo);
  setPropertyValue(bondBlock, mae::BOND_ATOM_2, numBonds, idx, bondFrom);
  setPropertyValue(bondBlock, mae::BOND_ORDER, numBonds, idx,
                   bondTypeToOrder(bond));

  if (dativeBondMark != nullptr) {
    dativeBondMark->set(idx, (bond.getBondType() == Bond::BondType::DATIVE));
  }

  // Custom properties
  copyProperties(bond, propNames, idx, boolSetter, intSetter, realSetter,
                 stringSetter);
}

void mapBonds(const ROMol& mol, const STR_VECT& propNames,
              mae::IndexedBlockMap& indexedBlockMap) {
  auto bondBlock = std::make_shared<mae::IndexedBlock>(mae::BOND_BLOCK);

  auto numBonds = mol.getNumBonds();

  std::shared_ptr<mae::IndexedProperty<mae::BoolProperty>> dativeBondMark =
      nullptr;
  for (auto& bond : mol.bonds()) {
    if (bond->getBondType() == Bond::BondType::DATIVE) {
      dativeBondMark = getIndexedProperty<mae::BoolProperty>(
          *bondBlock, MAE_BOND_DATIVE_MARK, numBonds);
      break;
    }
  }

  auto boolSetter = [&bondBlock, &numBonds](const std::string& prop,
                                            unsigned idx,
                                            mae::BoolProperty value) {
    setPropertyValue(*bondBlock, prop, numBonds, idx, value);
  };
  auto intSetter = [&bondBlock, &numBonds](const std::string& prop,
                                           unsigned idx, int value) {
    setPropertyValue(*bondBlock, prop, numBonds, idx, value);
  };
  auto realSetter = [&bondBlock, &numBonds](const std::string& prop,
                                            unsigned idx, double value) {
    setPropertyValue(*bondBlock, prop, numBonds, idx, value);
  };
  auto stringSetter = [&bondBlock, &numBonds](const std::string& prop,
                                              unsigned idx,
                                              const std::string& value) {
    setPropertyValue(*bondBlock, prop, numBonds, idx, value);
  };

  for (auto bond : mol.bonds()) {
    mapBond(*bond, dativeBondMark, propNames, *bondBlock, numBonds, boolSetter,
            intSetter, realSetter, stringSetter);
  }

  indexedBlockMap.addIndexedBlock(mae::BOND_BLOCK, bondBlock);
}

std::shared_ptr<mae::Block> _MolToMaeCtBlock(const ROMol& mol,
                                             const std::string& heavyAtomColor,
                                             int confId,
                                             const STR_VECT& propNames) {
  if (mol.getNumAtoms() == 0) {
    BOOST_LOG(rdErrorLog)
        << "ERROR: molecules without atoms cannot be exported to Maestro files.\n";
    return nullptr;
  }

  RWMol tmpMol(mol);
  if (!MolOps::KekulizeIfPossible(tmpMol)) {
    BOOST_LOG(rdErrorLog)
        << "ERROR: the mol cannot be kekulized, and will not be written to the output file.\n";
    return nullptr;
  }
  if (mol.getNumConformers() == 0) {
    // make sure there's at least one conformer we can write
    RDDepict::compute2DCoords(tmpMol);
  }

  auto stBlock = std::make_shared<mae::Block>(mae::CT_BLOCK);

  mapMolProperties(tmpMol, propNames, *stBlock);

  auto indexedBlockMap = std::make_shared<mae::IndexedBlockMap>();

  mapAtoms(tmpMol, propNames, heavyAtomColor, confId, *indexedBlockMap);

  if (mol.getNumBonds() > 0) {
    try {
      mapBonds(tmpMol, propNames, *indexedBlockMap);

    } catch (const UnsupportedBondException& exc) {
      BOOST_LOG(rdErrorLog)
          << "ERROR: " << exc.what()
          << " The mol will not be written to the output file.\n";
      return nullptr;
    }
  }

  stBlock->setIndexedBlockMap(indexedBlockMap);

  return stBlock;
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
  write(mol, defaultMaeHeavyAtomColor, confId);
}

void MaeWriter::write(const ROMol& mol, const std::string& heavyAtomColor,
                      int confId) {
  PRECONDITION(dp_ostream, "no output stream");

  auto block = _MolToMaeCtBlock(mol, heavyAtomColor, confId, d_props);

  if (block != nullptr) {
    if (!dp_writer) {
      open();
    }

    block->write(*dp_ostream);
    ++d_molid;
  }
}

std::string MaeWriter::getText(const ROMol& mol,
                               const std::string& heavyAtomColor, int confId,
                               const STR_VECT& propNames) {
  std::stringstream sstr;
  auto block = _MolToMaeCtBlock(mol, heavyAtomColor, confId, propNames);
  if (block != nullptr) {
    block->write(sstr);
  }

  return sstr.str();
}

}  // namespace RDKit
