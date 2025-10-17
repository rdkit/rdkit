//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include "GraphMolEnums.h"

#include "details.h"
#include "RingInfo.h"
#include "StereoGroup.h"
#include "SubstanceGroup.h"

#include <RDGeneral/export.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/rdmol_throw.h>
#include <GraphMol/MonomerInfo.h>
#include <Query/Query.h>

#include <boost/dynamic_bitset.hpp>

#include <atomic>
#include <limits>
#include <mutex>
#include <stdint.h>
#include <vector>

namespace RDKit {

class RDMol;
class ConstRDMolAtom;
class RDMolAtom;
class ConstRDMolBond;
class RDMolBond;

class Atom;
class Bond;
class Conformer;
class ROMol;
class RWMol;

class AtomData {
  static constexpr uint8_t unsetValenceVal =
      std::numeric_limits<uint8_t>::max();

  uint8_t atomicNum = 0;
  uint8_t numExplicitHs = 0;
  uint8_t explicitValence = unsetValenceVal;
  uint8_t implicitValence = unsetValenceVal;
  bool isAromatic = false;
  bool noImplicit = false;
  int8_t formalCharge = 0;
  uint8_t chiralTag = AtomEnums::ChiralType::CHI_UNSPECIFIED;
  uint64_t d_flags = 0;
  uint8_t numRadicalElectrons = 0;
  uint8_t hybridization = AtomEnums::HybridizationType::UNSPECIFIED;
  uint16_t isotope = 0;

  friend class RDMol;
  friend class MolPickler;
  friend class Atom;

public:
  uint32_t getAtomicNum() const { return atomicNum; }
  void setAtomicNum(uint32_t num) { atomicNum = num; }
  uint32_t getNumExplicitHs() const {
    return numExplicitHs;
  }
  void setNumExplicitHs(uint32_t num) {
    PRECONDITION(num <= unsetValenceVal,
                 "setNumExplicitHs called with too large number");
    numExplicitHs = num;
  }
  uint32_t getExplicitValence() const {
    PRECONDITION(explicitValence != unsetValenceVal,
                 "getExplicitValence() called without preceding call to "
                 "calcExplicitValence()");
    return explicitValence;
  }
  uint32_t getImplicitValence() const {
    // TODO: This PRECONDITION check was added to Atom::getValence, in the
    // commit adding Atom::getValence, but there's code in Chirality.cpp and
    // Atropisomers.cpp that checks for Atom::getImplicitValence returning -1,
    // which can only happen if the PRECONDITION check fails.
    // Is this correct?  Should the code checking for -1 be updated?
    PRECONDITION(noImplicit || implicitValence != unsetValenceVal,
                 "getImplicitValence() called without preceding call to "
                 "calcImplicitValence()");
    return noImplicit ? 0 : implicitValence;
  }
  uint32_t getTotalValence() const {
    return getExplicitValence() + getImplicitValence();
  }

  enum class ValenceType : std::uint8_t {
    IMPLICIT = 0,
    EXPLICIT
  };

  uint32_t getValence(ValenceType which) {
    if (which == ValenceType::EXPLICIT) {
      return getExplicitValence();
    }
    return getImplicitValence();
  }

  bool getNoImplicit() const { return noImplicit; }
  void setNoImplicit(bool value) { noImplicit = value; }

  // This is just a synonym of getImplicitValence
  uint32_t getNumImplicitHs() const {
    PRECONDITION(noImplicit || (implicitValence != unsetValenceVal),
                 "getNumImplicitHs() called without preceding call to "
                 "calcImplicitValence()");
    return getImplicitValence();
  }
  int8_t getFormalCharge() const {
    return formalCharge;
  }
  void setFormalCharge(int8_t charge) { formalCharge = charge; }

  uint32_t getIsotope() const { return isotope; }
  void setIsotope(uint32_t value) { isotope = value; }

  AtomEnums::ChiralType getChiralTag() const {
    return AtomEnums::ChiralType(chiralTag);
  }
  void setChiralTag(AtomEnums::ChiralType tag) {
    chiralTag = uint8_t(tag);
  }

  AtomEnums::HybridizationType getHybridization() const {
    return AtomEnums::HybridizationType(hybridization);
  }
  void setHybridization(AtomEnums::HybridizationType value) {
    hybridization = uint8_t(value);
  }

  bool getIsAromatic() const { return isAromatic; }
  void setIsAromatic(bool value) { isAromatic = value; }

  uint32_t getNumRadicalElectrons() const { return numRadicalElectrons; }
  void setNumRadicalElectrons(uint32_t num) { numRadicalElectrons = num; }

  bool needsUpdatePropertyCache() const {
    return explicitValence == unsetValenceVal ||
           (!noImplicit && implicitValence == unsetValenceVal);
  }
  void clearPropertyCache() {
    explicitValence = unsetValenceVal;
    implicitValence = unsetValenceVal;
  }

  double getMass() const;

  //! Flags that can be used by to store information on atoms.
  //!   These are not serialized and should be treated as temporary values.
  //!   No guarantees are made about preserving these flags across library
  //!   calls.
  void setFlags(std::uint64_t flags) { d_flags = flags; }
  std::uint64_t getFlags() const { return d_flags; }
  std::uint64_t &getFlags() { return d_flags; }
};

inline uint32_t getTwiceBondType(BondEnums::BondType type) {
  using BondEnums::BondType;
  switch (type) {
    case BondType::UNSPECIFIED:
    case BondType::IONIC:
    case BondType::ZERO:
      return 0;
    case BondType::SINGLE:
      return 2;
    case BondType::DOUBLE:
      return 4;
    case BondType::TRIPLE:
      return 6;
    case BondType::QUADRUPLE:
      return 8;
    case BondType::QUINTUPLE:
      return 10;
    case BondType::HEXTUPLE:
      return 12;
    case BondType::ONEANDAHALF:
      return 3;
    case BondType::TWOANDAHALF:
      return 5;
    case BondType::THREEANDAHALF:
      return 7;
    case BondType::FOURANDAHALF:
      return 9;
    case BondType::FIVEANDAHALF:
      return 11;
    case BondType::AROMATIC:
      return 3;
    case BondType::DATIVEONE:
      return 2; // FIX: this should probably be different
    case BondType::DATIVE:
      return 2; // FIX: again probably wrong
    case BondType::HYDROGEN:
      return 0;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
      return 0;
  }
}

class BondData {
  uint8_t bondType = BondEnums::BondType::UNSPECIFIED;
  uint8_t dirTag = BondEnums::BondDir::NONE;
  uint8_t stereo = BondEnums::BondStereo::STEREONONE;
  bool isAromatic = false;
  bool isConjugated = false;
  atomindex_t beginAtomIdx = atomindex_t(-1);
  atomindex_t endAtomIdx = atomindex_t(-1);
  atomindex_t stereoAtoms[2] = {atomindex_t(-1), atomindex_t(-1)};
  uint64_t d_flags = 0;

  friend class RDMol;
  friend class Bond;
 public:

  BondEnums::BondType getBondType() const {
    return BondEnums::BondType(bondType);
  }
  void setBondType(BondEnums::BondType type) { bondType = uint8_t(type); }

  BondEnums::BondStereo getStereo() const {
    return BondEnums::BondStereo(stereo);
  }
  void setStereo(BondEnums::BondStereo value) {
    //PRECONDITION(value <= BondEnums::BondStereo::STEREOE || hasStereoAtoms(),
    //             "Stereo atoms should be specified before specifying CIS/TRANS "
    //             "bond stereochemistry")
    stereo = uint8_t(value);
  }

  atomindex_t getBeginAtomIdx() const { return beginAtomIdx; }
  void setBeginAtomIdx(atomindex_t atomIdx) { beginAtomIdx = atomIdx; }
  atomindex_t getEndAtomIdx() const { return endAtomIdx; }
  void setEndAtomIdx(atomindex_t atomIdx) { endAtomIdx = atomIdx; }

  atomindex_t getOtherAtomIdx(const atomindex_t thisIdx) const {
    PRECONDITION(beginAtomIdx == thisIdx || endAtomIdx == thisIdx,
                 "bad index");
    return (beginAtomIdx == thisIdx) ? endAtomIdx : beginAtomIdx;
  }

  BondEnums::BondDir getBondDir() const { return BondEnums::BondDir(dirTag); }
  void setBondDir(BondEnums::BondDir dir) { dirTag = BondEnums::BondDir(dir); }

  bool getIsAromatic() const { return isAromatic; }
  void setIsAromatic(bool value) { isAromatic = value; }

  bool getIsConjugated() const { return isConjugated; }
  void setIsConjugated(bool value) { isConjugated = value; }

  uint32_t getTwiceValenceContrib(const atomindex_t atomIndex) const {
    if (atomIndex != beginAtomIdx && atomIndex != endAtomIdx) {
      return 0;
    }
    using BondEnums::BondType;
    if ((bondType == BondType::DATIVE || bondType == BondType::DATIVEONE) &&
        atomIndex != endAtomIdx) {
      return 0;
    } else if (bondType == BondType::HYDROGEN) {
      return 0;
    }
    return getTwiceBondType(BondType(bondType));
  }

  //! Flags that can be used by to store information on bonds.
  //!   These are not serialized and should be treated as temporary values.
  //!   No guarantees are made about preserving these flags across library
  //!   calls.
  void setFlags(std::uint64_t flags) { d_flags = flags; }
  std::uint64_t getFlags() const { return d_flags; }
  std::uint64_t &getFlags() { return d_flags; }
};

class RingInfoCache {
 public:
  // Indices into atomsInRings and bondsInRings where each ring begins
  // Size is numRings+1
  std::vector<uint32_t> ringBegins;
  // Size is the combined size of all rings (ringBegins.back())
  // Each ring starts at ringBegins[ringIndex]
  std::vector<uint32_t> atomsInRings;
  // Size is the combined size of all rings (ringBegins.back())
  // Each ring starts at ringBegins[ringIndex]
  // TODO: Consider combining bondsInRings with atomsInRings, since they're always kept the same length.
  std::vector<uint32_t> bondsInRings;
  // Size is numAtoms+1
  std::vector<uint32_t> atomMembershipBegins;
  // Size is the combined size of all rings (atomMembershipBegins.back())
  // Each atom's list of rings starts at atomMembershipBegins[atomIndex]
  std::vector<uint32_t> atomMemberships;
  // Size is numBonds+1
  std::vector<uint32_t> bondMembershipBegins;
  // Size is the combined size of all rings (bondMembershipBegins.back())
  // Each bond's list of rings starts at bondMembershipBegins[bondIndex]
  std::vector<uint32_t> bondMemberships;

  // 2D bit vector representing whether each pair of rings shares at least one bond.
  // Size is numRings x numRings, in one contiguous bit vector
  std::vector<bool> areRingsFused;
  // The number of bonds in each ring that are shared with any other rings.
  // Size is numRings
  std::vector<uint32_t> numFusedBonds;

  // TODO: Add support for ring families

  INT_VECT tempPerAtomInts;
  std::vector<std::pair<uint32_t, uint32_t>> tempPerAtomPairs;

  RDKit::FIND_RING_TYPE findRingType = FIND_RING_TYPE_OTHER_OR_UNKNOWN;

  bool isInit = false;

  bool isInitialized() const { return isInit; }

  void reset() {
    ringBegins.resize(0);
    atomsInRings.resize(0);
    bondsInRings.resize(0);
    atomMembershipBegins.resize(0);
    atomMemberships.resize(0);
    bondMembershipBegins.resize(0);
    bondMemberships.resize(0);
    areRingsFused.resize(0);
    numFusedBonds.resize(0);
    tempPerAtomInts.resize(0);
    tempPerAtomPairs.resize(0);
    findRingType = FIND_RING_TYPE_OTHER_OR_UNKNOWN;
    isInit = false;
  }

  uint32_t numRings() const {
    return (ringBegins.size() == 0) ? 0 : (ringBegins.size() - 1);
  }

  uint32_t numAtomRings(atomindex_t atomIndex) const {
    return atomMembershipBegins[atomIndex + 1] -
           atomMembershipBegins[atomIndex];
  }

  uint32_t numBondRings(uint32_t bondIndex) const {
    return bondMembershipBegins[bondIndex + 1] -
           bondMembershipBegins[bondIndex];
  }

  bool isAtomInRingOfSize(atomindex_t atomIndex, uint32_t size) const {
    const uint32_t* atomMembershipBegin =
        atomMemberships.data() + atomMembershipBegins[atomIndex];
    const uint32_t* atomMembershipEnd =
        atomMemberships.data() + atomMembershipBegins[atomIndex + 1];
    for (; atomMembershipBegin != atomMembershipEnd; ++atomMembershipBegin) {
      const uint32_t ring = *atomMembershipBegin;
      if (ringBegins[ring + 1] - ringBegins[ring] == size) {
        return true;
      }
    }
    return false;
  }

  bool isBondInRingOfSize(uint32_t bondIndex, uint32_t size) const {
    const uint32_t* bondMembershipBegin =
        bondMemberships.data() + bondMembershipBegins[bondIndex];
    const uint32_t* bondMembershipEnd =
        bondMemberships.data() + bondMembershipBegins[bondIndex + 1];
    for (; bondMembershipBegin != bondMembershipEnd; ++bondMembershipBegin) {
      const uint32_t ring = *bondMembershipBegin;
      if (ringBegins[ring + 1] - ringBegins[ring] == size) {
        return true;
      }
    }
    return false;
  }

  uint32_t minAtomRingSize(atomindex_t atomIndex) const {
    uint32_t minSize = std::numeric_limits<uint32_t>::max();
    const uint32_t* atomMembershipBegin =
        atomMemberships.data() + atomMembershipBegins[atomIndex];
    const uint32_t* atomMembershipEnd =
        atomMemberships.data() + atomMembershipBegins[atomIndex + 1];
    for (; atomMembershipBegin != atomMembershipEnd; ++atomMembershipBegin) {
      const uint32_t ring = *atomMembershipBegin;
      const uint32_t size = ringBegins[ring + 1] - ringBegins[ring];
      if (size < minSize) {
        minSize = size;
      }
    }
    return minSize;
  }

  uint32_t minBondRingSize(uint32_t bondIndex) const {
    uint32_t minSize = std::numeric_limits<uint32_t>::max();
    const uint32_t* bondMembershipBegin =
        bondMemberships.data() + bondMembershipBegins[bondIndex];
    const uint32_t* bondMembershipEnd =
        bondMemberships.data() + bondMembershipBegins[bondIndex + 1];
    for (; bondMembershipBegin != bondMembershipEnd; ++bondMembershipBegin) {
      const uint32_t ring = *bondMembershipBegin;
      const uint32_t size = ringBegins[ring + 1] - ringBegins[ring];
      if (size < minSize) {
        minSize = size;
      }
    }
    return minSize;
  }
};

struct RDKIT_GRAPHMOL_EXPORT StereoGroups {
  static constexpr uint32_t undefinedGroupId = 0;

  // 1 per group
  std::vector<StereoGroupType> stereoTypes;
  std::vector<uint32_t> readIds;
  std::vector<uint32_t> writeIds;

  // Size n_groups + 1
  std::vector<uint32_t> atomBegins = {0};
  std::vector<uint32_t> bondBegins = {0};

  // Size atomBegins[n_atoms]
  std::vector<uint32_t> atoms;
  // Size bondBegins[n_bonds]
  std::vector<uint32_t> bonds;


  //! Returns the number of stereo groups.
  uint32_t getNumGroups() const {return stereoTypes.size();}
  //! Removes a stereo group by index.
  void removeGroup(uint32_t groupIndex);

  void addGroup(StereoGroupType type, const std::vector<uint32_t>& atomIndices,
                const std::vector<uint32_t>& bondIndices, uint32_t readId = undefinedGroupId);

  StereoGroupType getGroupType(uint32_t groupIndex) const {
    URANGE_CHECK(groupIndex, stereoTypes.size());
    return stereoTypes[groupIndex];
  }

  uint32_t getReadID(uint32_t groupIndex) const {
    URANGE_CHECK(groupIndex, readIds.size());
    return readIds[groupIndex];
  }

  uint32_t getWriteID(uint32_t groupIndex) const {
    URANGE_CHECK(groupIndex, writeIds.size());
    return writeIds[groupIndex];
  }
  void setWriteID(uint32_t groupIndex, uint32_t value) {
    URANGE_CHECK(groupIndex, writeIds.size());
    writeIds[groupIndex] = value;
  }

  std::pair<uint32_t*, uint32_t*> getAtoms(uint32_t groupIndex) {
    URANGE_CHECK(groupIndex + 1, atomBegins.size());
    return {atoms.data() + atomBegins[groupIndex], atoms.data() + atomBegins[groupIndex + 1]};
  }
  std::pair<const uint32_t*, const uint32_t*> getAtoms(uint32_t groupIndex) const {
    URANGE_CHECK(groupIndex + 1, atomBegins.size());
    return {atoms.data() + atomBegins[groupIndex], atoms.data() + atomBegins[groupIndex + 1]};
  }
  std::pair<uint32_t*, uint32_t*> getBonds(uint32_t groupIndex) {
    URANGE_CHECK(groupIndex + 1, bondBegins.size());
    return {bonds.data() + bondBegins[groupIndex], bonds.data() + bondBegins[groupIndex + 1]};
  }
  std::pair<const uint32_t*, const uint32_t*> getBonds(uint32_t groupIndex) const {
    URANGE_CHECK(groupIndex + 1, bondBegins.size());
    return {bonds.data() + bondBegins[groupIndex], bonds.data() + bondBegins[groupIndex + 1]};
  }

  //! Compare two groups for equality. Does not compare read and write IDs.
  bool GroupsEq(uint32_t index1, uint32_t index2);
  //! Removes atom index from all groups where present. If decrementIndices is true,
  //! all atom indices greater than the removed index are decremented.
  void removeAtomFromGroups(uint32_t atomIndex, bool decrementIndices = false);
  //! Removes bond index from all groups where present. If decrementIndices is
  //! true, all bond indices greater than the removed index are decremented.
  void removeBondFromGroups(uint32_t bondIndex, bool decrementIndices = false);
  //! Remove any group containing the specified atom.
  void removeGroupsWithAtom(uint32_t atomIndex);
  //! Remove any group containing any of the specified atoms.
  void removeGroupsWithAtoms(const std::vector<uint32_t>& atomIndices);
  //! Remove any group containing the specified bond.
  void removeGroupsWithBond(uint32_t bondIndex);
  //! Remove any group containing any of the specified bonds.
  void removeGroupsWithBonds(const std::vector<uint32_t>& bondIndices);
  //! Populates write IDs for stereo groups that don't have them.
  void assignStereoGroupIds();
  //! Populates all write IDs from read IDs.
  void forwardStereoGroupIds();
};

struct ConformerIdAndFlags {
  uint32_t id;
  bool is3D;
};

enum class PropertyType : uint8_t {
  CHAR,     // 8 bit
  BOOL,     // 8 bit,
  INT32,    // 32 bit
  UINT32,   // 32 bit
  INT64,    // 64 bit
  UINT64,   // 64 bit
  FLOAT,    // 32 bit
  DOUBLE,   // 64 bit
  ANY       // RDValue
};

template <typename T>
struct TypeToPropertyType {
  static constexpr PropertyType family = PropertyType::ANY;
};
template<> struct TypeToPropertyType<bool> { static constexpr PropertyType family = PropertyType::BOOL;};
template<> struct TypeToPropertyType<char> { static constexpr PropertyType family = PropertyType::CHAR;};
template<> struct TypeToPropertyType<int> { static constexpr PropertyType family = PropertyType::INT32;};
template<> struct TypeToPropertyType<unsigned int> { static constexpr PropertyType family = PropertyType::UINT32;};
template<> struct TypeToPropertyType<int64_t> { static constexpr PropertyType family = PropertyType::INT64;};
template<> struct TypeToPropertyType<uint64_t> { static constexpr PropertyType family = PropertyType::UINT64;};
template<> struct TypeToPropertyType<float> { static constexpr PropertyType family = PropertyType::FLOAT;};
template<> struct TypeToPropertyType<double> { static constexpr PropertyType family = PropertyType::DOUBLE;};

inline PropertyType RDValueTagToPropertyType(
    decltype(std::declval<RDValue>().getTag()) tag) {
  switch (tag) {
    case RDTypeTag::BoolTag:
      return PropertyType::BOOL;
    case RDTypeTag::IntTag:
      return PropertyType::INT32;
    case RDTypeTag::UnsignedIntTag:
      return PropertyType::UINT32;
    case RDTypeTag::FloatTag:
      return PropertyType::FLOAT;
    case RDTypeTag::DoubleTag:
      return PropertyType::DOUBLE;
    default:
      return PropertyType::ANY;
  }
}

//! Variant for holding scalar array properties
struct RDKIT_GRAPHMOL_EXPORT PropArray {
  //! Number of elements
  uint32_t size = 0;
  //! Scalar family
  PropertyType family;
  //! Data pointer
  void* data = nullptr;
  //! Mask for set values (bool);
  std::unique_ptr<bool[]> isSetMask = nullptr;
  uint32_t numSet = 0;

  PropArray();
  PropArray(uint32_t size, PropertyType family, bool isSet);
  PropArray(const PropArray& other);
  PropArray& operator=(const PropArray& other);
  PropArray(PropArray&& other);
  PropArray& operator=(PropArray&& other);
  ~PropArray() noexcept;

  //! Set up array data.
  void construct(bool isSet);

  //! Tear down array data
  void destroy();

  //! Cast an array data point to an RDValue containing that data point.
  RDValue toRDValue(uint32_t idx) const;

  template<typename T>
  T getValueAs(uint32_t index) const {
    using RawType = std::remove_cv_t<std::remove_reference_t<T>>;
    // There are 4 combinations of scalar/RDvalue request/content we need to
    // account for.
    // Note that any invalid conversions, such as int to vector of int, will be
    // caught via RDValue cast.
    if constexpr (TypeToPropertyType<RawType>::family == PropertyType::ANY) {
      // Case where we need an RDValue and have one already.
      if (family == PropertyType::ANY) {
        RDValue &caster = static_cast<RDValue *>(data)[index];
        if constexpr (std::is_same_v<RawType, std::string>) {
          std::string res;
          rdvalue_tostring(caster, res);
          return res;
        } else {
          return from_rdvalue<T>(caster);
        }
      }
      // Case where we need an RDValue but have a scalar.
      else {
        RDValue caster = toRDValue(index);
        if constexpr (std::is_same_v<RawType, std::string>) {
          std::string res;
          try {
            rdvalue_tostring(caster, res);
            RDValue::cleanup_rdvalue(caster);
            return res;
          } catch (...) {
            RDValue::cleanup_rdvalue(caster);
            throw;
          }
        } else {
          try {
            T res = rdvalue_cast<T>(caster);
            RDValue::cleanup_rdvalue(caster);
            return res;
          } catch (...) {
            RDValue::cleanup_rdvalue(caster);
            throw;
          }
        }
      }
    } else {
      // Case where we need a scalar but have an RDValue.
      if (family == PropertyType::ANY) {
        RDValue &caster = static_cast<RDValue *>(data)[index];
        return from_rdvalue<T>(caster);
      }
      // Case where we have a scalar and need a scalar.
      else {
        if (family == TypeToPropertyType<RawType>::family) {
          return static_cast<T *>(data)[index];
        }
        RDValue tempVal = toRDValue(index);
        try {
          T res = rdvalue_cast<T>(tempVal);
          RDValue::cleanup_rdvalue(tempVal);
          return res;
        } catch (...) {
            RDValue::cleanup_rdvalue(tempVal);
            throw;
        };

      }
    }
  }

  PropertyType getType(uint32_t index) const {
    if (family != PropertyType::ANY) {
      return family;
    }
    return RDValueTagToPropertyType(
        static_cast<RDValue *>(data)[index].getTag());
  }

  // For code requesting values as if from an RDValue to know
  // what type to request.
  decltype(std::declval<RDValue>().getTag()) getRDValueTag(
      uint32_t index) const {
    switch (family) {
      case PropertyType::CHAR:
        // TODO: Should char be treated as an integer or a string?
        // String would require changes in PropArray::toRDValue,
        // for getValueAs to work correctly.
        return RDTypeTag::AnyTag;
      case PropertyType::BOOL:
        return RDTypeTag::BoolTag;
      case PropertyType::INT32:
        return RDTypeTag::IntTag;
      case PropertyType::UINT32:
        return RDTypeTag::UnsignedIntTag;
      case PropertyType::INT64:
        // Be cautious of possible integer overflow
        return RDTypeTag::IntTag;
      case PropertyType::UINT64:
        // Be cautious of possible integer overflow
        return RDTypeTag::UnsignedIntTag;
      case PropertyType::FLOAT:
        return RDTypeTag::FloatTag;
      case PropertyType::DOUBLE:
        return RDTypeTag::DoubleTag;
      case PropertyType::ANY:
        return static_cast<RDValue *>(data)[index].getTag();
    }
  }

  //! Adds a single new property to the end of the array.
  void appendElement();
  //! Removes an element from the array.
  void removeElement(uint32_t index);

  //! Converts storage type to RDValue, family PropertyType::ANY,
  //! or does nothing if already RDValue.
  void convertToRDValue();
};

namespace Ranges {
template <bool ISBOND, bool ISCONST>
class IndexRange;
template <bool ISBOND, bool ISCONST>
class IndirectRange;
}  // namespace Ranges

class RDKIT_GRAPHMOL_EXPORT RDMol {
  std::vector<AtomData> atomData;
  std::vector<uint32_t> atomBondStarts = {0u};
  std::vector<uint32_t> otherAtomIndices;
  std::vector<uint32_t> bondDataIndices;
  // Mutable only for compatibility synchronization in getBondStereoAtoms below
  mutable std::vector<BondData> bondData;

  std::unique_ptr<StereoGroups> stereoGroups;
  // Mutable only for compatibility synchronization while locked.
  mutable RingInfoCache ringInfo;

  // FIXME: Sort out a more efficient representation for substance groups.
  std::vector<SubstanceGroup> substanceGroups;

  size_t numConformers = 0;
  size_t conformerAtomCapacity = 0;

  // All conformers are stored in conformerPositionData, which consists of
  // numConformers x conformerAtomCapacity x 3 doubles.
  // It's resized when adding conformers or when the number of atoms exceeds
  // conformerAtomCapacity.
  // Mutable only for compatibility synchronization while locked.
  mutable std::vector<double> conformerPositionData;

  // This is unused if all conformer IDs are equal to their index and
  // all conformers are 3D.
  // Mutable only for compatibility synchronization while locked.
  mutable std::vector<ConformerIdAndFlags> conformerIdsAndFlags;

public:
  using QUERYATOM_QUERY = Queries::Query<int, ConstRDMolAtom, true>;
  using QUERYBOND_QUERY = Queries::Query<int, ConstRDMolBond, true>;

private:
  std::vector<std::unique_ptr<QUERYATOM_QUERY>> atomQueries;
  std::vector<std::unique_ptr<QUERYBOND_QUERY>> bondQueries;

  bool isStereochemDone = false;

  //! Atom tracker for batch edits. A batch edit is in progress if non-null.
  std::unique_ptr<boost::dynamic_bitset<>> dp_delAtoms;
  //! Bond tracker for batch edits. A batch edit is in progress if non-null.
  std::unique_ptr<boost::dynamic_bitset<>> dp_delBonds;

  struct CompatibilityData;
  // CompatibilityData needs private access to Conformer for setOwningMol,
  // but CompatibilityData itself is private, so Conformer needs private access here.
  friend class Conformer;
  mutable std::mutex compatibilityMutex;
  mutable std::atomic<CompatibilityData*> compatibilityData = nullptr;

 public:
  enum class Scope : uint8_t {
    MOL,   // Just 1 value (can be an array value)
    ATOM,  // 1 value per atom
    BOND   // 1 value per bond
  };

  class RDKIT_GRAPHMOL_EXPORT PropIterator;

  class RDKIT_GRAPHMOL_EXPORT Property {
    PropToken d_name;
    Scope d_scope;
    bool d_isComputed;
    RDValue d_inPlaceData;
    PropArray d_arrayData;

    friend class RDKit::RDMol::PropIterator;
    friend class RDKit::RDMol;

    // Don't access the constructors or assignment operators outside of RDMol.
    // They're public so that the standard library functions and classes have
    // access, for example std::vector and its use of std::move.
   public:
    // Default initialization just needs to indicate a type with no allocation.
    Property() : d_name(), d_scope(Scope::MOL), d_isComputed(false) {}
    Property(Property &&other) noexcept
        : d_name(std::move(other.d_name)),
          d_scope(other.d_scope),
          d_isComputed(other.d_isComputed),
          d_inPlaceData(std::move(other.d_inPlaceData)),
          d_arrayData(std::move(other.d_arrayData)) {
      // If a move constructor is ever added to RDValue, setting the type here
      // may not be required and may be invalid, but the static_assert will
      // force someone to check what the correct behaviour should be.
      static_assert(std::is_trivially_move_constructible_v<RDValue>);
      if (std::is_trivially_move_constructible_v<RDValue>) {
        other.d_inPlaceData.type = RDTypeTag::EmptyTag;
      }
    }
    ~Property() { RDValue::cleanup_rdvalue(d_inPlaceData); }
    Property& operator=(Property&& other) noexcept {
      if (this == &other) {
        return *this;
      }
      d_name = std::move(other.d_name);
      d_scope = other.d_scope;
      d_isComputed = other.d_isComputed;
      d_inPlaceData = std::move(other.d_inPlaceData);
      // If a move constructor is ever added to RDValue, setting the type here
      // may not be required and may be invalid, but the static_assert will
      // force someone to check what the correct behaviour should be.
      static_assert(std::is_trivially_move_constructible_v<RDValue>);
      if (std::is_trivially_move_constructible_v<RDValue>) {
        other.d_inPlaceData.type = RDTypeTag::EmptyTag;
      }
      d_arrayData = std::move(other.d_arrayData);
      return *this;
    }
    Property(const Property &other)
        : d_name(other.d_name),
          d_scope(other.d_scope),
          d_isComputed(other.d_isComputed),
          d_arrayData(other.d_arrayData) {
      copy_rdvalue(d_inPlaceData, other.d_inPlaceData);
    }
    Property &operator=(const Property &other) {
      if (this == &other) {
        return *this;
      }
      d_name = other.d_name;
      d_scope = other.d_scope;
      d_isComputed = other.d_isComputed;
      copy_rdvalue(d_inPlaceData, other.d_inPlaceData);
      d_arrayData = other.d_arrayData;
      return *this;
    }

    const PropToken &name() const { return d_name; }
    Scope scope() const { return d_scope; }
    bool isComputed() const { return d_isComputed; }

    template<typename T>
    void setVal(const T& val) {
      PRECONDITION(d_scope == Scope::MOL, "In-place data only for MOL scope");
      RDValue::cleanup_rdvalue(d_inPlaceData);
      d_inPlaceData = val;
    }
    void setVal(const RDValue &val) {
      PRECONDITION(d_scope == Scope::MOL, "In-place data only for MOL scope");
      copy_rdvalue(d_inPlaceData, val);
    }
    void setVal(const char *val) {
      std::string s(val);
      setVal(s);
    }

    template<typename T>
    T getVal() {
      PRECONDITION(d_scope == Scope::MOL, "In-place data only for MOL scope");
      return rdvalue_cast<T>(d_inPlaceData);
    }

    PropertyType getRawType() const {
      if (d_scope == Scope::MOL) {
        return PropertyType::ANY;
      }
      return d_arrayData.family;
    }

    PropertyType getType() const {
      PRECONDITION(d_scope == Scope::MOL, "In-place data only for MOL scope");
      return RDValueTagToPropertyType(d_inPlaceData.getTag());
    }
    PropertyType getType(uint32_t index) const {
      PRECONDITION(d_scope != Scope::MOL, "Array data only for non-MOL scope");
      return d_arrayData.getType(index);
    }

    auto getRDValueTag() const {
      PRECONDITION(d_scope == Scope::MOL, "In-place data only for MOL scope");
      return d_inPlaceData.getTag();
    }
    auto getRDValueTag(uint32_t index) const {
      PRECONDITION(d_scope != Scope::MOL, "Array data only for non-MOL scope");
      return d_arrayData.getRDValueTag(index);
    }
  };

  class RDKIT_GRAPHMOL_EXPORT PropIterator {
    std::vector<Property>::const_iterator d_current;
    std::vector<Property>::const_iterator d_end;
    bool d_includeComputed;
    Scope d_scope;
    uint32_t d_index;

    bool currentScopeIsCorrect() const {
      if (d_current->scope() != d_scope) {
        return false;
      }
      if (!d_includeComputed && d_current->isComputed()) {
        return false;
      }
      if (d_scope != Scope::MOL && d_index != anyIndexMarker &&
          !d_current->d_arrayData.isSetMask[d_index]) {
        return false;
      }
      return true;
    }

   public:
    static constexpr uint32_t anyIndexMarker = uint32_t(-1);

    using difference_type = size_t;
    using value_type = PropToken;
    using pointer = PropToken*;
    using reference = PropToken&;
    using iterator_category = std::forward_iterator_tag;

    PropIterator() = default;
    PropIterator(PropIterator&&) = default;
    PropIterator(const PropIterator&) = default;
    PropIterator &operator=(PropIterator&&) = default;
    PropIterator &operator=(const PropIterator&) = default;

    PropIterator(std::vector<Property>::const_iterator current,
                 std::vector<Property>::const_iterator end,
                 bool includeComputed, Scope scope, uint32_t index)
        : d_current(current),
          d_end(end),
          d_includeComputed(includeComputed),
          d_scope(scope),
          d_index(index) {
      // Iterate to the first matching property.
      while (d_current != d_end && !currentScopeIsCorrect()) {
        ++d_current;
      }
    }

    const Property &operator*() const {
      PRECONDITION(d_current != d_end, "PropIterator::operator* end");
      PRECONDITION(d_current->scope() == d_scope,
                   "PropIterator::operator* scope");
      PRECONDITION(d_includeComputed || !d_current->isComputed(),
                   "PropIterator::operator* isComputed");
      PRECONDITION(d_scope == Scope::MOL || d_index == anyIndexMarker ||
                       (d_index < d_current->d_arrayData.size &&
                        d_current->d_arrayData.isSetMask[d_index]),
                   "PropIterator::operator* index")
      return *d_current;
    }
    const Property *operator->() const { return &**this; }

    PropIterator &operator++() {
      PRECONDITION(d_current != d_end, "PropIterator::operator++ end");
      PRECONDITION(d_scope == Scope::MOL || d_index == anyIndexMarker ||
                       d_index < d_current->d_arrayData.size,
                   "PropIterator::operator++ index");
      while (true) {
        ++d_current;
        if (d_current == d_end) {
          break;
        }
        if (currentScopeIsCorrect()) {
          break;
        }
      }
      return *this;
    }

    PropIterator operator++(int) {
      auto copy = *this;
      ++(*this);
      return copy;
    }

    bool operator==(const PropIterator &that) const {
      return d_current == that.d_current;
    }
    bool operator!=(const PropIterator &that) const { return !(*this == that); }
  };

 private:

  std::vector<Property> properties;

  //! Map of bookmarks to atom indices associated with that bookmark
  std::unordered_map<int, std::vector<int>> atomBookmarks;
  //! Map of bookmarks to bond indices associated with that bookmark
  std::unordered_map<int, std::vector<int>> bondBookmarks;

  std::unordered_map<int, std::unique_ptr<AtomMonomerInfo>> monomerInfo;
  friend class ROMol;
  friend class RWMol;
  friend class Atom;
  friend class Bond;
  friend class QueryAtom;
  friend class QueryBond;

 public:
  using ADJ_ITER_PAIR = std::pair<const uint32_t*,const uint32_t*>;

  RDMol() = default;
  explicit RDMol(const RDMol &other);
  RDMol(const RDMol &other, bool quickCopy, int confId = -1);
  RDMol &operator=(const RDMol &other);
  RDMol(RDMol &&other) noexcept;
  RDMol &operator=(RDMol &&other) noexcept;

  ~RDMol();
  void clear();

  // Returns the number of explicit atoms
  uint32_t getNumAtoms() const {
    return uint32_t(atomData.size());
  }
  uint32_t getNumAtoms(bool onlyExplicit) const {
    const uint32_t numAtoms = uint32_t(atomData.size());
    uint32_t res = numAtoms;
    if (!onlyExplicit) {
      for (uint32_t atom = 0; atom < numAtoms; ++atom) {
        res += getTotalNumHs(atom, false);
      }
    }
    return res;
  }
  uint32_t getNumHeavyAtoms() const {
    const uint32_t numAtoms = uint32_t(atomData.size());
    uint32_t res = 0;
    for (uint32_t i = 0; i != numAtoms; ++i) {
      if (atomData[i].getAtomicNum() > 1) {
        ++res;
      }
    }
    return res;
  };

  const AtomData& getAtom(uint32_t atomIndex) const {
    URANGE_CHECK(atomIndex, getNumAtoms());
    return atomData[atomIndex];
  }
  AtomData& getAtom(uint32_t atomIndex) {
    URANGE_CHECK(atomIndex, getNumAtoms());
    return atomData[atomIndex];
  }
  const BondData& getBond(uint32_t bondIndex) const {
    URANGE_CHECK(bondIndex, getNumBonds());
    return bondData[bondIndex];
  }
  BondData& getBond(uint32_t bondIndex) {
    URANGE_CHECK(bondIndex, getNumBonds());
    return bondData[bondIndex];
  }

  uint32_t getTotalNumHs(uint32_t atomIndex, bool includeNeighbors = false) const {
    PRECONDITION(atomIndex < getNumAtoms(), "RDMol::getTotalNumHs called with out of bounds atomIndex");
    const AtomData& atom = atomData[atomIndex];
    uint32_t res = atom.getNumExplicitHs() + atom.getNumImplicitHs();
    if (includeNeighbors) {
      uint32_t atomBegin = atomBondStarts[atomIndex];
      const uint32_t atomEnd = atomBondStarts[atomIndex+1];
      for (; atomBegin < atomEnd; ++atomBegin) {
        const AtomData& nbr = atomData[otherAtomIndices[atomBegin]];
        if (nbr.getAtomicNum() == 1) {
          ++res;
        }
      }
    }
    return res;
  }
  uint32_t getAtomDegree(uint32_t atomIndex) const {
    URANGE_CHECK(atomIndex, getNumAtoms());
    uint32_t atomBegin = atomBondStarts[atomIndex];
    uint32_t atomEnd = atomBondStarts[atomIndex + 1];
    return atomEnd - atomBegin;
  }
  uint32_t getAtomTotalDegree(uint32_t atomIndex) const {
    return getTotalNumHs(atomIndex, false) + getAtomDegree(atomIndex);
  }
  uint32_t getNumBonds(bool onlyHeavy = true) const {
    // By default return the bonds that connect only the heavy atoms
    // hydrogen connecting bonds are ignores
    uint32_t numBonds = uint32_t(bondData.size());
    if (!onlyHeavy) {
      // If we need hydrogen connecting bonds add them up
      for (uint32_t atom = 0, numAtoms = getNumAtoms(); atom < numAtoms; ++atom) {
        numBonds += getTotalNumHs(atom, false);
      }
    }
    return numBonds;
  }

  bool isAromaticAtom(atomindex_t atomIndex) const {
    PRECONDITION(atomIndex < getNumAtoms(),
                 "RDMol::isAromaticAtom called with out of bounds atomIndex");
    const AtomData& atom = atomData[atomIndex];
    if (atom.getIsAromatic()) {
      return true;
    }
    uint32_t atomBegin = atomBondStarts[atomIndex];
    const uint32_t atomEnd = atomBondStarts[atomIndex + 1];
    for (; atomBegin < atomEnd; ++atomBegin) {
      const BondData& bond = bondData[bondDataIndices[atomBegin]];
      if (bond.getIsAromatic() ||
          bond.getBondType() == BondEnums::BondType::AROMATIC) {
        return true;
      }
    }

    return false;
  }

  // These bond stereo atoms functions are on RDMol instead of BondData so that
  // it has access to the compatibility fallback.  This fallback is necessary
  // because the original interface returns references to INT_VECT, and the
  // compatibility interface needs to still provide that without breaking
  // the access to values via this new interface.

  void clearBondStereoAtoms(uint32_t bondIndex) {
    URANGE_CHECK(bondIndex, getNumBonds());
    BondData& bond = bondData[bondIndex];
    bond.stereoAtoms[0] = atomindex_t(-1);
    bond.stereoAtoms[1] = atomindex_t(-1);
    if (hasCompatibilityData()) {
      clearBondStereoAtomsCompat(bondIndex);
    }
  }
  bool hasBondStereoAtoms(uint32_t bondIndex) const {
    URANGE_CHECK(bondIndex, getNumBonds());
    const BondData& bond = bondData[bondIndex];
    if (!hasCompatibilityData()) {
      PRECONDITION(
          (bond.stereoAtoms[0] == atomindex_t(-1)) ==
              (bond.stereoAtoms[1] == atomindex_t(-1)),
          "BondData::hasStereoAtoms called and only one valid stereo atom index");
      return (bond.stereoAtoms[0] != atomindex_t(-1));
    }
    return hasBondStereoAtomsCompat(bondIndex);
  }
  void setBondStereoAtoms(uint32_t bondIndex, atomindex_t a, atomindex_t b, bool checkAtomIndices = true) {
    URANGE_CHECK(bondIndex, getNumBonds());
    if (checkAtomIndices) {
      PRECONDITION(a != atomindex_t(-1) && b != atomindex_t(-1),
                   "BondData::setStereoAtoms called with invalid indices");
      URANGE_CHECK(a, getNumAtoms());
      URANGE_CHECK(b, getNumAtoms());
    }
    BondData& bond = bondData[bondIndex];
    bond.stereoAtoms[0] = a;
    bond.stereoAtoms[1] = b;
    if (hasCompatibilityData()) {
      auto *compatStereo = getBondStereoAtomsCompat(bondIndex);
      CHECK_INVARIANT(compatStereo, "No stereo atoms in compatibility data");
      compatStereo->clear();
      compatStereo->push_back(a);
      compatStereo->push_back(b);
    }
  }
  const atomindex_t* getBondStereoAtoms(uint32_t bondIndex) const {
    URANGE_CHECK(bondIndex, getNumBonds());
    if (hasCompatibilityData()) {
      const INT_VECT* data = getBondStereoAtomsCompat(bondIndex);
      PRECONDITION(data, "bond stereo compat data not initialized");
      if (data->size() > 0) {
        PRECONDITION(data->size() == 2, "bond stereo compat data not size 2");
        static_assert(sizeof(atomindex_t) == sizeof((*data)[0]));
        return reinterpret_cast<const atomindex_t*>(data->data());
      }

      // If the non-const overload of Bond::getStereoAtoms is called and the
      // caller clears the vector, this will fall back to BondData::stereoAtoms,
      // but it must be set to -1 to indicate being empty.  This should be
      // threadsdafe for reading cases, because all threads will write the same
      // values here and won't return until they've written the values.
      BondData &bond = bondData[bondIndex];
      bond.stereoAtoms[0] = atomindex_t(-1);
      bond.stereoAtoms[1] = atomindex_t(-1);
      return bond.stereoAtoms;
    }
    const BondData& bond = bondData[bondIndex];
    return bond.stereoAtoms;
  }

  uint32_t calcExplicitValence(atomindex_t atomIndex, bool strict);
  uint32_t calcImplicitValence(atomindex_t atomIndex, bool strict);

  bool hasValenceViolation(atomindex_t atomIndex) const;

  //! inverts the atom's \c chiralTag, returns whether or not a change was made
  bool invertAtomChirality(atomindex_t atomIndex);

  ADJ_ITER_PAIR getAtomNeighbors(atomindex_t atomIndex) const {
    URANGE_CHECK(atomIndex, getNumAtoms());
    uint32_t atomBegin = atomBondStarts[atomIndex];
    uint32_t atomEnd = atomBondStarts[atomIndex + 1];
    return ADJ_ITER_PAIR(otherAtomIndices.data() + atomBegin,
                         otherAtomIndices.data() + atomEnd);
  }
  ADJ_ITER_PAIR getAtomBonds(atomindex_t atomIndex) const {
    URANGE_CHECK(atomIndex, getNumAtoms());
    uint32_t atomBegin = atomBondStarts[atomIndex];
    uint32_t atomEnd = atomBondStarts[atomIndex + 1];
    return ADJ_ITER_PAIR(bondDataIndices.data() + atomBegin,
                         bondDataIndices.data() + atomEnd);
  }
  uint32_t getBondIndexBetweenAtoms(atomindex_t atom0, atomindex_t atom1) const {
    URANGE_CHECK(atom0, getNumAtoms());
    URANGE_CHECK(atom1, getNumAtoms());
    uint32_t atomBegin = atomBondStarts[atom0];
    const uint32_t atomEnd = atomBondStarts[atom0 + 1];
    for (; atomBegin < atomEnd; ++atomBegin) {
      uint32_t otherAtomIndex = otherAtomIndices[atomBegin];
      if (otherAtomIndex == atom1) {
        return bondDataIndices[atomBegin];
      }
    }
    return std::numeric_limits<std::uint32_t>::max();
  }
  const uint32_t* getAtomBondStarts() const {
    return atomBondStarts.data();
  }
  const uint32_t* getOtherAtomIndices() const {
    return otherAtomIndices.data();
  }
  const uint32_t* getBondDataIndices() const {
    return bondDataIndices.data();
  }

  // For RDKit use only
  std::vector<AtomData>& getAtomDataVector() { return atomData; }
  // For RDKit use only
  std::vector<BondData>& getBondDataVector() { return bondData; }
  // For RDKit use only
  std::vector<uint32_t>& getAtomBondStartsVector() { return atomBondStarts; }
  // For RDKit use only
  std::vector<uint32_t>& getOtherAtomIndicesVector() { return otherAtomIndices; }
  // For RDKit use only
  std::vector<uint32_t>& getBondDataIndicesVector() { return bondDataIndices; }

  bool hasProp(const PropToken& name) const;
  bool hasAtomProp(const PropToken& name, uint32_t index) const;
  bool hasBondProp(const PropToken& name, uint32_t index) const;

  std::vector<std::string> getPropList(
      bool includePrivate = true,
      bool includeComputed = true,
      Scope scope = Scope::MOL,
      uint32_t index = PropIterator::anyIndexMarker) const;

private:
  // This is used by ROMol, Atom, and Bond to get the "__computedProps" prop
  void getComputedPropList(
      STR_VECT &res, Scope scope = Scope::MOL,
      uint32_t index = PropIterator::anyIndexMarker) const;

public:

  PropIterator beginProps(bool includeComputed = true, Scope scope = Scope::MOL,
                          uint32_t index = PropIterator::anyIndexMarker) const {
    return PropIterator(properties.cbegin(), properties.cend(), includeComputed,
                        scope, index);
  }
  PropIterator endProps() const {
    return PropIterator(properties.cend(), properties.cend(), true, Scope::MOL,
                        PropIterator::anyIndexMarker);
  }

  template <typename T>
  T getMolProp(const PropToken &name) const {
    const Property *prop = findProp(name, Scope::MOL);
    if (prop == nullptr) {
      throw KeyErrorException(name.getString());
    }
    return from_rdvalue<T>(prop->d_inPlaceData);
  }

  template <typename T>
  T getAtomProp(const PropToken &name, uint32_t index) const {
    const Property *prop = findProp(name, Scope::ATOM);
    if (prop == nullptr || !prop->d_arrayData.isSetMask[index]) {
      throw KeyErrorException(name.getString());
    }
    return prop->d_arrayData.getValueAs<T>(index);
  }

  template <typename T>
  T getBondProp(const PropToken &name, uint32_t index) const {
    const Property *prop = findProp(name, Scope::BOND);
    if (prop == nullptr || !prop->d_arrayData.isSetMask[index]) {
      throw KeyErrorException(name.getString());
    }
    return prop->d_arrayData.getValueAs<T>(index);
  }

  //! Populates res with the property with token name, if present. Returns whether token was found.
  template <typename T>
  bool getMolPropIfPresent(const PropToken& name, T& res) const {
    const Property* prop = findProp(name, Scope::MOL);
    if (prop == nullptr) {
      return false;
    }
    if constexpr (std::is_same_v<T, std::string>) {
      rdvalue_tostring(prop->d_inPlaceData, res);
    } else {
      res = from_rdvalue<T>(prop->d_inPlaceData);
    }
    return true;
  }

  //! Populates res with the property with token name, if present. Returns whether token was found.
  template <typename T>
  bool getAtomPropIfPresent(const PropToken& name, uint32_t index, T& res) const {
    return getArrayPropIfPresent<T>(name, Scope::ATOM, index, res);
  }

  //! Populates res with the property with token name, if present. Returns whether token was found.
  template <typename T>
  bool getBondPropIfPresent(const PropToken& name, uint32_t index, T& res) const {
    return getArrayPropIfPresent<T>(name, Scope::BOND, index, res);
  }

  //! Returns pointer to array of atom properties with token name, or nullptr if not present.
  //! Returns nullptr if the property is not a basic scalar type (use getAtomPropIfPresent).
  template <typename T> T* getAtomPropArrayIfPresent(const PropToken& name) {
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    const Property* prop = findProp(name, Scope::ATOM);
    if (prop == nullptr) {
      return nullptr;
    }
    if (prop->d_arrayData.family != TypeToPropertyType<T>::family) {
      return nullptr;
    }
    return static_cast<T *>(prop->d_arrayData.data);
  }

  //! \overload
  template <typename T> const T* getAtomPropArrayIfPresent(const PropToken& name) const {
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    const Property* prop = findProp(name, Scope::ATOM);
    if (prop == nullptr) {
      return nullptr;
    }
    if (prop->d_arrayData.family != TypeToPropertyType<T>::family) {
      return nullptr;
    }
    return static_cast<T *>(prop->d_arrayData.data);
  }

  //! Returns pointer to array of bond properties with token name, or nullptr if not present.
  //! Returns nullptr if the property is not a basic scalar type (use getBondPropIfPresent).
  template <typename T> T* getBondPropArrayIfPresent(const PropToken& name) {
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    const Property* prop = findProp(name, Scope::BOND);
    if (prop == nullptr) {
      return nullptr;
    }
    if (prop->d_arrayData.family != TypeToPropertyType<T>::family) {
      return nullptr;
    }
    return static_cast<T *>(prop->d_arrayData.data);
  }

  //! \overload
  template <typename T> const T* getBondPropArrayIfPresent(const PropToken& name) const {
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    const Property* prop = findProp(name, Scope::BOND);
    if (prop == nullptr) {
      return nullptr;
    }
    if (prop->d_arrayData.family != TypeToPropertyType<T>::family) {
      return nullptr;
    }
    return static_cast<T *>(prop->d_arrayData.data);
  }

  //! Adds a property with token name and value to the molecule. If the property already exists, it is overwritten.
  template <typename T> void addMolProp(const PropToken& name, const T& value, bool computed = false) {
    if (Property *prop = findProp(name, Scope::MOL); prop != nullptr) {
      prop->setVal(value);
      prop->d_isComputed = computed;
      return;
    }
    Property& newProp = properties.emplace_back();
    newProp.d_name = std::move(name);
    newProp.d_scope = Scope::MOL;
    newProp.d_isComputed = computed;
    newProp.setVal(value);
  }

  //! Alias of addMolProp
  template <typename T> void setMolProp(const PropToken& name, const T& value, bool computed = false) {
    addMolProp(name, value, computed);
  }

  //! Adds a property array with token name and value for all atoms
  //! Returns a pointer to the array if a scalar type, or nullptr if not.
  template <typename T>
  T* addAtomProp(const PropToken& name, T defaultValue, bool computed = false) {
    Property * newProp = addPerElementProp<T>(name, Scope::ATOM, defaultValue, computed);
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    return static_cast<T *>(newProp->d_arrayData.data);
  }

  //! Adds a property array with token name and value for all bonds
  //! Returns a pointer to the array if a scalar type, or nullptr if not.
  template<typename T>
  T* addBondProp(const PropToken& name, T defaultValue, bool computed = false) {
    Property * newProp = addPerElementProp<T>(name, Scope::BOND, defaultValue, computed);
    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      return nullptr;
    }
    return static_cast<T *>(newProp->d_arrayData.data);
  }

  //! Sets an atom property at the given index
  template<typename T>
  void setSingleAtomProp(const PropToken &name, const std::uint32_t index,
                         const T &value, bool computed = false,
                         bool supportTypeMismatch = false) {
    setSingleProp<T>(name, Scope::ATOM, index, value, computed,
                     supportTypeMismatch);
  }

  //! Sets a bond property at the given index
  template<typename T>
  void setSingleBondProp(const PropToken &name, const std::uint32_t index,
                         const T &value, bool computed = false,
                         bool supportTypeMismatch = false) {
    setSingleProp<T>(name, Scope::BOND, index, value, computed,
                     supportTypeMismatch);
  }

  //! Clears a molecule property if present
  void clearMolPropIfPresent(const PropToken& name) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ++it) {
      if (it->name() == name && it->scope() == Scope::MOL) {
        properties.erase(it);
        return;
      }
    }
  }

  //! Clears all atom props with the given name
  void clearAtomPropIfPresent(const PropToken& name) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ++it) {
      if (it->name() == name && it->scope() == Scope::ATOM) {
        properties.erase(it);
        return;
      }
    }
  }

  //! Clears all bond props with the given name
  void clearBondPropIfPresent(const PropToken& name) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ++it) {
      if (it->name() == name && it->scope() == Scope::BOND) {
        properties.erase(it);
        return;
      }
    }
  }

  //! Clears a single atom prop at the given index
  void clearSingleAtomProp(const PropToken& name, const std::uint32_t atomIndex) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ++it) {
      if (it->name() == name && it->scope() == Scope::ATOM &&
          it->d_arrayData.isSetMask[atomIndex]) {
        --it->d_arrayData.numSet;
        if (it->d_arrayData.numSet == 0) {
          properties.erase(it);
          return;
        }
        it->d_arrayData.isSetMask[atomIndex] = false;
        return;
      }
    }
  }

  //! Clears a single bond prop at the given index
  void clearSingleBondProp(const PropToken& name, const std::uint32_t bondIndex) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ++it) {
      if (it->name() == name && it->scope() == Scope::BOND &&
          it->d_arrayData.isSetMask[bondIndex]) {
        --it->d_arrayData.numSet;
        if (it->d_arrayData.numSet == 0) {
          properties.erase(it);
          return;
        }
        it->d_arrayData.isSetMask[bondIndex] = false;
        return;
      }
    }
  }

  //! Clears all props for a given atom
  void clearSingleAtomAllProps(const std::uint32_t atomIndex, const bool onlyComputed = false) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ) {
      bool erased = false;
      if (it->scope() == Scope::ATOM && (!onlyComputed || it->isComputed()) &&
          it->d_arrayData.isSetMask[atomIndex]) {
        --it->d_arrayData.numSet;
        if (it->d_arrayData.numSet == 0) {
          it = properties.erase(it);
          end = properties.end();
          erased = true;
        } else {
          it->d_arrayData.isSetMask[atomIndex] = false;
        }
      }
      if (!erased) {
        ++it;
      }
    }
  }

  //! Clears all props for a given bond
  void clearSingleBondAllProps(const std::uint32_t bondIndex, const bool onlyComputed = false) {
    for (auto it = properties.begin(), end = properties.end(); it != end; ) {
      bool erased = false;
      if (it->scope() == Scope::BOND && (!onlyComputed || it->isComputed()) &&
          it->d_arrayData.isSetMask[bondIndex]) {
        --it->d_arrayData.numSet;
        if (it->d_arrayData.numSet == 0) {
          it = properties.erase(it);
          end = properties.end();
          erased = true;
        } else {
          it->d_arrayData.isSetMask[bondIndex] = false;
        }
      }
      if (!erased) {
        ++it;
      }
    }
  }
  //! Clear all mol, atom, and bond properties.
  void clearProps();

  void copyProp(const PropToken &destinationName, const RDMol &sourceMol,
                const PropToken &sourceName, Scope scope = Scope::MOL);
  void copySingleProp(const PropToken &destinationName,
                           uint32_t destinationIndex, const RDMol &sourceMol,
                           const PropToken &sourceName, uint32_t sourceIndex,
                           Scope scope);

  //! Associates an atom index with a bookmark
  void setAtomBookmark(int atomIdx, int mark);
  //! Associates an atom index with a bookmark, removing any previous bookmark with the same mark
  void replaceAtomBookmark(int atomIdx, int mark);
  //! Returns a vector of atom indices associated with a bookmark
  std::vector<int> getAllAtomsWithBookmarks(int mark);
  //! Returns the first atom index associated with a bookmark. This will be the first atom bookmarked with the given mark.
  int getAtomWithBookmark(int mark);
  //! Removes a bookmark.
  void clearAtomBookmark(int mark);
  //! Removes a bookmark from a specific atom.
  void clearAtomBookmark(int atomIdx, int mark);
  //! Removes all atom bookmarks.
  void clearAllAtomBookmarks();
  //! Returns whether a bookmark exists.
  bool hasAtomBookmark(int mark) const;

  //! Associates a bond index with a bookmark
  void setBondBookmark(int bondIdx, int mark);
  //! Returns a vector of Bond indices associated with a bookmark
  std::vector<int> getAllBondsWithBookmarks(int mark);
  //! Returns the first bond index associated with a bookmark. This will be the first bond bookmarked with the given mark.
  int getBondWithBookmark(int mark);
  //! Removes a bookmark.
  void clearBondBookmark(int mark);
  //! Removes a bookmark from a specific bond.
  void clearBondBookmark(int bondIdx, int mark);
  //! Removes all bond bookmarks.
  void clearAllBondBookmarks();
  //! Returns whether a bookmark exists.
  bool hasBondBookmark(int mark) const;

  const StereoGroups* getStereoGroups() const { return stereoGroups.get(); }
  StereoGroups* getStereoGroups() { return stereoGroups.get(); }
  void setStereoGroups(std::unique_ptr<StereoGroups> &&groups);

  std::vector<SubstanceGroup>& getSubstanceGroups() { return substanceGroups; }
  const std::vector<SubstanceGroup>& getSubstanceGroups() const { return substanceGroups; }

  bool getIsStereochemDone() const { return isStereochemDone; }
  void setIsStereochemDone(bool value) { isStereochemDone = value; }

  //! \name Properties
  //! @{

  //! clears all of our \c computed \c properties
  void clearComputedProps(bool includeRings = true);
  //! calculates any of our lazy \c properties
  /*!
    <b>Notes:</b>
       - this calls \c updatePropertyCache() on each of our Atoms and Bonds
  */
  void updatePropertyCache(bool strict = true);

  void updateAtomPropertyCache(atomindex_t atomIndex, bool strict = true) {
    calcExplicitValence(atomIndex, strict);
    calcImplicitValence(atomIndex, strict);
  }
  void updateBondPropertyCache([[maybe_unused]] uint32_t bondIndex,
                               [[maybe_unused]] bool strict = true) {
    // Currently nothing
  }

  bool needsUpdatePropertyCache() const;
  void clearPropertyCache();

  //! @}

  int getAtomMapNum(uint32_t atomIndex) const {
    const int* atomMap = getAtomPropArrayIfPresent<int>(
        common_properties::molAtomMapNumberToken);
    return (atomMap != nullptr) ? atomMap[atomIndex] : 0;
  }

  const QUERYATOM_QUERY *getAtomQuery(uint32_t atomIndex) const {
    URANGE_CHECK(atomIndex, getNumAtoms());
    PRECONDITION(atomQueries.size() == 0 || atomQueries.size() == getNumAtoms(),
                 "atomQueries invalid size");
    return (atomQueries.size() == 0) ? nullptr : atomQueries[atomIndex].get();
  }
  QUERYATOM_QUERY *getAtomQuery(uint32_t atomIndex) {
    URANGE_CHECK(atomIndex, getNumAtoms());
    PRECONDITION(atomQueries.size() == 0 || atomQueries.size() == getNumAtoms(),
                 "atomQueries invalid size");
    return (atomQueries.size() == 0) ? nullptr : atomQueries[atomIndex].get();
  }
  const QUERYBOND_QUERY *getBondQuery(uint32_t bondIndex) const {
    URANGE_CHECK(bondIndex, getNumBonds());
    PRECONDITION(bondQueries.size() == 0 || bondQueries.size() == getNumBonds(),
                 "bondQueries invalid size");
    return (bondQueries.size() == 0) ? nullptr : bondQueries[bondIndex].get();
  }
  QUERYBOND_QUERY *getBondQuery(uint32_t bondIndex) {
    URANGE_CHECK(bondIndex, getNumBonds());
    PRECONDITION(bondQueries.size() == 0 || bondQueries.size() == getNumBonds(),
                 "bondQueries invalid size");
    return (bondQueries.size() == 0) ? nullptr : bondQueries[bondIndex].get();
  }
  bool hasAtomQuery(uint32_t atomIndex) const {
    return getAtomQuery(atomIndex) != nullptr;
  }
  bool hasBondQuery(uint32_t bondIndex) const {
    return getBondQuery(bondIndex) != nullptr;
  }

  //! Add a new atom to the molecule.
  AtomData& addAtom();
  //! Add a new bond to the molecule between the atoms with the given indices.
  //! If updateAromaticity is true, the atoms will be marked as aromatic if the bond is aromatic.
  BondData& addBond(uint32_t beginAtomIdx, uint32_t endAtomIdx, BondEnums::BondType bondType, bool updateAromaticity = true);
  void removeAtom(atomindex_t atomIndex, bool clearProps = true);
  void removeBond(uint32_t bondIndex);

  //! Begin batch edit. All changes will be made upon commit.
  void beginBatchEdit();
  //! Commit batch edit to apply changes.
  void commitBatchEdit();
  //! Rollback batch edit to discard changes.
  void rollbackBatchEdit();

  //! Returns the number of conformers.
  uint32_t getNumConformers() const;
  //! Return the coordinates of the conformer with the given ID.
  const double* getConformerPositions(int32_t id = -1) const;
  //! Return the coordinates of the conformer with the given ID.
  double* getConformerPositions(int32_t id = -1);
  //! Remove all conformers.
  //! Return whether the conformer with the given ID is 3D.
  bool is3DConformer(uint32_t id) const;
  void clearConformers();
  //! Remove the given conformer ID if it exists;
  void removeConformer(uint32_t id);
  //! Add a conformer, return the index. Assigns an ID if id < 0;
  uint32_t addConformer(const double* positions, int32_t id = -1, bool is3D = true);
  //! Adds an unset conformer, returns the index and data pointer. Assigns an ID if id < 0;
  std::pair<uint32_t, double*> addConformer(int32_t id = -1, bool is3D = true);
  //! Allocates memory for n conformers. Helpful for performance when adding many conformers or atoms.
  void allocateConformers(size_t newNumConformers) {
    allocateConformers(newNumConformers, getNumAtoms());
  }
  //! \overload
  void allocateConformers(size_t newNumConformers, size_t newNumAtoms);

  const RingInfoCache &getRingInfo() const;
  RingInfoCache &getRingInfo();

  void setAtomQuery(atomindex_t atomIndex,
                    std::unique_ptr<QUERYATOM_QUERY> query);
  void setBondQuery(uint32_t bondIndex,
                    std::unique_ptr<QUERYBOND_QUERY> query);

  // Compatibity interface access
  ROMol& asROMol();
  const ROMol& asROMol() const;
  RWMol& asRWMol();
  const RWMol& asRWMol() const;

  Ranges::IndexRange<false, false> atoms();
  Ranges::IndexRange<false, true> atoms() const;
  Ranges::IndexRange<true, false> bonds();
  Ranges::IndexRange<true, true> bonds() const;
  Ranges::IndirectRange<false, false> atomNeighbors(atomindex_t atomIndex);
  Ranges::IndirectRange<false, true> atomNeighbors(atomindex_t atomIndex) const;
  Ranges::IndirectRange<true, false> atomBonds(atomindex_t atomIndex);
  Ranges::IndirectRange<true, true> atomBonds(atomindex_t atomIndex) const;

 private:
  //! Instantiate an RDMol with a stable ROMol pointer.
  RDMol(ROMol* existingPtr);
  //! Copy constructor with stable ROMol pointer.
  RDMol(const RDMol &other, ROMol *existingPtr);
  //! \overload
  RDMol(const RDMol &other, bool quickCopy, int confId, ROMol *existingPtr);

  void initFromOther(const RDMol &other, bool quickCopy, int confId, ROMol *existingPtr);

  //! Used by Bond to initialize single-bond non-owning molecule. Avoids atom index bounds checking.
  void addUnconnectedBond() {bondData.resize(1);}

  //! Returns handle to property if found, else nullptr
  const Property* findProp(const PropToken& name, Scope scope) const {
    for (const auto &property : properties) {
      if (property.name() == name && property.scope() == scope) {
        return &property;
      }
    }
    return nullptr;
  }

  //! Non-const overload
  Property* findProp(const PropToken& name, Scope scope) {
    for (auto &property : properties) {
      if (property.name() == name && property.scope() == scope) {
        return &property;
      }
    }
    return nullptr;
  }

  //! Adds a new array property for bonds or atoms.
  template <typename T>
  Property* addPerElementProp(const PropToken& name, Scope scope, const T& defaultValue, bool computed = false, bool setAll = true) {
    // Throw if not atom or bond scope
    if (scope != Scope::ATOM && scope != Scope::BOND) {
      throw ValueErrorException("Per-element properties must be atom or bond scope");
    }
    Property* existing = findProp(name, scope);
    if (existing != nullptr) {
      if (existing->d_arrayData.family == TypeToPropertyType<T>::family) {
        existing->d_isComputed = computed;
        return existing;
      }
      // If we have an existing prop with the wrong type, clear it.
      if (scope == Scope::ATOM) {
        clearAtomPropIfPresent(name);
      } else {
        clearBondPropIfPresent(name);
      }
    }

    const uint32_t numElements = scope == Scope::ATOM ? getNumAtoms() : getNumBonds();
    Property& newProperty = properties.emplace_back();
    newProperty.d_name = name;
    newProperty.d_scope = scope;
    newProperty.d_isComputed = computed;
    newProperty.d_arrayData.size = numElements;
    newProperty.d_arrayData.family = TypeToPropertyType<T>::family;
    newProperty.d_arrayData.construct(setAll);

    if constexpr (TypeToPropertyType<T>::family == PropertyType::ANY) {
      auto *data = static_cast<RDValue *>(newProperty.d_arrayData.data);
      for (uint32_t i = 0; i < numElements; ++i) {
        if constexpr (std::is_same_v<T, RDValue>) {
          copy_rdvalue(data[i], defaultValue);
        } else if constexpr (std::is_same_v<T, const char*>) {
          data[i] = std::string(defaultValue);
        } else {
          data[i] = defaultValue;
        }
      }
    } else {
      auto *data = static_cast<T *>(newProperty.d_arrayData.data);
      for (uint32_t i = 0; i < numElements; ++i) {
        data[i] = defaultValue;
      }
    }
    return &newProperty;
  }

  // Overload to treat const char* as a string
  Property *addPerElementProp(const PropToken &name, Scope scope,
                              const char *defaultValue, bool computed = false,
                              bool setAll = true) {
    std::string s = defaultValue == nullptr ? "" : defaultValue;
    return addPerElementProp(name, scope, s, computed, setAll);
  }


  template<typename T>
  void setSingleProp(const PropToken &name, const Scope scope,
                     const std::uint32_t index, const T &value,
                     bool computed = false, bool supportTypeMismatch = false) {
    Property* existing = findProp(name, scope);
    if (existing == nullptr) {
      if constexpr (std::is_default_constructible_v<T>) {
        existing = addPerElementProp(name, scope, T(), computed, false);
      } else {
        existing = addPerElementProp(name, scope, value, computed, false);
      }
    }
    existing->d_isComputed = computed;
    auto existingFamily = existing->d_arrayData.family;
    if (existingFamily != TypeToPropertyType<T>::family && existingFamily != PropertyType::ANY) {
      if (!supportTypeMismatch) {
        throw ValueErrorException("Property type mismatch");
      }
      // Convert to RDValue
      existing->d_arrayData.convertToRDValue();
      PRECONDITION(existing->d_arrayData.family == PropertyType::ANY,
                   "convertToRDValue should make family ANY");
      existingFamily = PropertyType::ANY;
    }
    auto &arr = existing->d_arrayData;
    if (TypeToPropertyType<T>::family == PropertyType::ANY || existingFamily == PropertyType::ANY) {
      auto* data = static_cast<RDValue*>(arr.data);
      // Default value may still need cleaning even if not previously set.
      RDValue::cleanup_rdvalue(data[index]);
      if constexpr (std::is_same_v<T, RDValue>) {
        copy_rdvalue(data[index], value);
      } else if constexpr (std::is_same_v<T, const char*>) {
        data[index] = std::string(value);
      } else {
        data[index] = value;
      }
    } else {
      auto *data = static_cast<T *>(arr.data);
      data[index] = value;
    }
    arr.isSetMask[index] = true;
    ++arr.numSet;
  }

  // Overload to treat const char* as a string
  void setSingleProp(const PropToken &name, const Scope scope,
                     const std::uint32_t index, const char *value,
                     bool computed = false, bool supportTypeMismatch = false) {
    std::string s(value);
    setSingleProp(name, scope, index, s, computed, supportTypeMismatch);
  }

  template <typename T>
  bool getArrayPropIfPresent(const PropToken& name, Scope scope, std::uint32_t index, T& res) const {
    const Property* prop = findProp(name, scope);
    if (prop == nullptr) {
      return false;
    }

    if (scope == Scope::ATOM && index >= getNumAtoms()) {
      throw ValueErrorException("Atom index out of bounds");
    }

    if (scope == Scope::BOND && index >= getNumBonds()) {
      throw ValueErrorException("Bond index out of bounds");
    }

    const PropArray &arr = prop->d_arrayData;
    if (arr.isSetMask[index]) {
      res = arr.getValueAs<T>(index);
      return true;
    }
    return false;
  }

  //! Gets the index of a conformer based on its ID, or throws
  //! ConformerException if not found. Negative id returns index 0.
  uint32_t findConformerIndex(int32_t id) const;
  //! Pushes positions and metadata from conformers to compatibility data.
  //! If confId is -1, all conformers are copied.
  //! If confId is positive, only the conformer with that id is copied,
  //! or none, if no conformer exists with that id.
  void copyConformersFromCompatibilityData(const CompatibilityData* compatData, int confId = -1) const;
  //! Pulls positions and metadata from compatibility data to conformers.
  void copyConformersToCompatibilityData(CompatibilityData* compatData) const;
  //! Manually tell rdmol that the compat data conformers may have been modified.
  void markConformersAsCompatModified() const;

  int calculateExplicitValence(atomindex_t atomIndex, bool strict,
                               bool checkIt) const;
  int calculateImplicitValence(atomindex_t atomIndex, bool strict,
                               bool checkIt) const;

  // NOTE: This can change in another thread, so only use it in situations where
  // the caller is planning to modify data or where the result changing would
  // only be relevant conditional upon data being modified.
  bool hasCompatibilityData() const {
    return compatibilityData.load(std::memory_order_relaxed) != nullptr;
  }
  CompatibilityData* getCompatibilityDataIfPresent() {
    return compatibilityData.load(std::memory_order_relaxed);
  }
  const CompatibilityData* getCompatibilityDataIfPresent() const {
    return compatibilityData.load(std::memory_order_relaxed);
  }
  CompatibilityData* ensureCompatInit() const;
  void clearBondStereoAtomsCompat(uint32_t bondIndex);
  bool hasBondStereoAtomsCompat(uint32_t bondIndex) const;
  const INT_VECT* getBondStereoAtomsCompat(uint32_t bondIndex) const;
  INT_VECT* getBondStereoAtomsCompat(uint32_t bondIndex);
  Atom* getAtomCompat(uint32_t atomIndex);
  const Atom* getAtomCompat(uint32_t atomIndex) const;
  Bond* getBondCompat(uint32_t bondIndex);
  const Bond* getBondCompat(uint32_t bondIndex) const;
  void releaseCompatMolOwnership();
  //! Set atom, bond, conformer pointers to be owned by this RDMol. Required for
  //! move operation ownership transfer.
  void setCompatPointersToSelf();
  //! Sets the compatibililty data ROMol pointer.
  void setROMolPointerCompat(ROMol* ptr);

  std::list<Atom *>& getAtomBookmarksCompat(int mark);
  std::list<Bond *>& getBondBookmarksCompat(int mark);

  std::map<int, std::list<Atom *>>* getAtomBookmarksCompat();
  std::map<int, std::list<Bond *>>* getBondBookmarksCompat();


  std::list<boost::shared_ptr<Conformer>>& getConformersCompat();
  const std::list<boost::shared_ptr<Conformer>>& getConformersCompat() const;

  RingInfo &getRingInfoCompat();
  const RingInfo &getRingInfoCompat() const;

  const std::vector<StereoGroup> &getStereoGroupsCompat();

  // Replaces an atom handle in compatibility data. Cleans up previous data mol.
  // Does not copy any data, but replaces pointers in compatibility structures.
  void replaceAtomPointerCompat(Atom* atom, uint32_t index);
  // Replaces a bond handle in compatibility data. Cleans up previous data mol.
  // Does not copy any data, but replaces pointers in compatibility structures.
  void replaceBondPointerCompat(Bond* atom, uint32_t index);

  // source can be owned by a different RDMol, but null means to copy from
  // this RDMol's compatibilityData.
  // If quickCopy is true, bookmarks and conformers are cleared and not copied.
  // If confId is -1, all conformers are copied.
  // If confId is positive, only the conformer with that id is copied,
  // or none, if no conformer exists with that id.
  void copyFromCompatibilityData(const CompatibilityData *source,
                                 bool quickCopy = false, int confId = -1);

  // This only wraps query and sets the wrapper into atomQueries.
  // Called from QueryAtom, which must have already updated its Query
  void setAtomQueryCompat(atomindex_t atomIndex, Queries::Query<int, const Atom *, true> *query);

  // This only wraps query and sets the wrapper into bondQueries.
  // Called from QueryBond, which must have already updated its Query
  void setBondQueryCompat(uint32_t bondIndex, Queries::Query<int, const Bond *, true> *query);

  // This forces the specified atom to be a QueryAtom, sets its query to this, taking ownership,
  // and then also wraps it and sets the wrapper into atomQueries.
  void setAtomQueryCompatFull(atomindex_t atomIndex,
                              Queries::Query<int, const Atom *, true> *query);

  // This forces the specified bond to be a QueryBond, sets its query to this, taking ownership,
  // and then also wraps it and sets the wrapper into bondQueries.
  void setBondQueryCompatFull(uint32_t bondIndex,
                              Queries::Query<int, const Bond *, true> *query);
};


template<bool ISCONST>
class RDMolWrapperBase {
 protected:
  using MolType = std::conditional_t<ISCONST, const RDMol, RDMol>;

 private:
  MolType* d_mol;
  uint32_t d_index;

 protected:
  RDMolWrapperBase() : d_mol(nullptr), d_index(0) {}
  RDMolWrapperBase(MolType* mol, uint32_t index)
      : d_mol(mol), d_index(index) {}
  RDMolWrapperBase(const RDMolWrapperBase&) = default;
  RDMolWrapperBase& operator=(const RDMolWrapperBase&) = default;
 public:
  bool isValid() const { return (d_mol != nullptr); }
  MolType& mol() { return *d_mol; }
  const RDMol& mol() const { return *d_mol; }
  uint32_t index() const { return d_index; }
};

template<bool ISCONST>
class RDMolAtomBase : public RDMolWrapperBase<ISCONST> {
  using MolType = typename RDMolWrapperBase<ISCONST>::MolType;
 protected:
  RDMolAtomBase() : RDMolWrapperBase<ISCONST>() {}
  RDMolAtomBase(MolType* mol, uint32_t atomIndex) : RDMolWrapperBase<ISCONST>(mol, atomIndex) {}
  RDMolAtomBase(const RDMolAtomBase&) = default;
  RDMolAtomBase& operator=(const RDMolAtomBase&) = default;
};

class ConstRDMolAtom : public RDMolAtomBase<true> {
 public:
  ConstRDMolAtom() : RDMolAtomBase() {}
  ConstRDMolAtom(const RDMol* mol, uint32_t atomIndex) : RDMolAtomBase(mol, atomIndex) {}
  ConstRDMolAtom(const ConstRDMolAtom&) = default;
  ConstRDMolAtom& operator=(const ConstRDMolAtom&) = default;
  ConstRDMolAtom(const RDMolAtom&);

  const AtomData& data() const { return mol().getAtom(index()); }
};

class RDMolAtom : public RDMolAtomBase<false> {
 public:
  RDMolAtom() : RDMolAtomBase() {}
  RDMolAtom(RDMol* mol, uint32_t atomIndex) : RDMolAtomBase(mol, atomIndex) {}
  RDMolAtom(const RDMolAtom&) = default;
  RDMolAtom& operator=(const RDMolAtom&) = default;

  AtomData& data() { return mol().getAtom(index()); }
};

inline ConstRDMolAtom::ConstRDMolAtom(const RDMolAtom& other)
    : RDMolAtomBase(&other.mol(), other.index()) {}

template<bool ISCONST>
class RDMolBondBase : public RDMolWrapperBase<ISCONST> {
  using MolType = typename RDMolWrapperBase<ISCONST>::MolType;
 protected:
  RDMolBondBase() : RDMolWrapperBase<ISCONST>() {}
  RDMolBondBase(MolType* mol, uint32_t bondIndex) : RDMolWrapperBase<ISCONST>(mol, bondIndex) {}
  RDMolBondBase(const RDMolBondBase&) = default;
  RDMolBondBase& operator=(const RDMolBondBase&) = default;
};

class ConstRDMolBond : public RDMolBondBase<true> {
 public:
  ConstRDMolBond() : RDMolBondBase() {}
  ConstRDMolBond(const RDMol* mol, uint32_t bondIndex) : RDMolBondBase(mol, bondIndex) {}
  ConstRDMolBond(const ConstRDMolBond&) = default;
  ConstRDMolBond& operator=(const ConstRDMolBond&) = default;
  ConstRDMolBond(const RDMolBond&);

  const BondData& data() const { return mol().getBond(index()); }
};

class RDMolBond : public RDMolBondBase<false> {
 public:
  RDMolBond() : RDMolBondBase() {}
  RDMolBond(RDMol* mol, uint32_t bondIndex) : RDMolBondBase(mol, bondIndex) {}
  RDMolBond(const RDMolBond&) = default;
  RDMolBond& operator=(const RDMolBond&) = default;

  BondData& data() { return mol().getBond(index()); }
};

inline ConstRDMolBond::ConstRDMolBond(const RDMolBond& other)
    : RDMolBondBase(&other.mol(), other.index()) {}

namespace Ranges {
//! This class is only intended to be used by IndexRange
template <bool ISBOND, bool ISCONST>
class IndexIter {
  using MolType = std::conditional_t<ISCONST, const RDMol, RDMol>;
  MolType *d_mol;
  uint32_t d_index;

  using value_type = std::conditional_t<ISBOND,
      std::conditional_t<ISCONST, ConstRDMolBond, RDMolBond>,
      std::conditional_t<ISCONST, ConstRDMolAtom, RDMolAtom>>;
  using iterator_category = std::forward_iterator_tag;

  friend class IndexRange<ISBOND, ISCONST>;

  IndexIter(MolType *mol, uint32_t index) : d_mol(mol), d_index(index) {}

 public:
  IndexIter(const IndexIter &) = default;
  IndexIter &operator=(const IndexIter &) = default;

  value_type operator*() const { return value_type(d_mol, d_index); }

  IndexIter &operator++() {
    ++d_index;
    return *this;
  }
  IndexIter operator++(int) {
    IndexIter copy = *this;
    ++d_index;
    return copy;
  }
  bool operator==(const IndexIter &other) const {
    return d_index == other.d_index && d_mol == other.d_mol;
  }
  bool operator!=(const IndexIter &other) const { return !(*this == other); }
};

//! This class is only intended to be used by range-based for loops via
//! RDMol::atoms() and  RDMol::bonds()
template <bool ISBOND, bool ISCONST>
class IndexRange {
  using MolType = std::conditional_t<ISCONST, const RDMol, RDMol>;
  using Iter = IndexIter<ISBOND, ISCONST>;
  MolType *d_mol;

  friend class RDKit::RDMol;
  IndexRange(MolType *mol) : d_mol(mol) {}

 public:
  IndexRange(const IndexRange<ISBOND, ISCONST> &) = default;
  IndexRange &operator=(const IndexRange<ISBOND, ISCONST> &) = default;
  Iter begin() const { return Iter(d_mol, 0); }
  Iter end() const { return Iter(d_mol, d_mol->getNumAtoms()); }
};

//! This class is only intended to be used by IndirectRange
template <bool ISBOND, bool ISCONST>
class IndirectIter {
  using MolType = std::conditional_t<ISCONST, const RDMol, RDMol>;
  MolType *d_mol;
  const uint32_t *d_indices;

  using value_type = std::conditional_t<ISBOND,
      std::conditional_t<ISCONST, ConstRDMolBond, RDMolBond>,
      std::conditional_t<ISCONST, ConstRDMolAtom, RDMolAtom>>;
  using iterator_category = std::forward_iterator_tag;

  friend class IndirectRange<ISBOND, ISCONST>;

  IndirectIter(MolType *mol, const uint32_t *indices)
      : d_mol(mol), d_indices(indices) {}

 public:
  IndirectIter(const IndirectIter &) = default;
  IndirectIter &operator=(const IndirectIter &) = default;

  value_type operator*() const { return value_type(d_mol, *d_indices); }

  IndirectIter &operator++() {
    ++d_indices;
    return *this;
  }
  IndirectIter operator++(int) {
    IndirectIter copy = *this;
    ++d_indices;
    return copy;
  }
  bool operator==(const IndirectIter &other) const {
    return d_indices == other.d_indices && d_mol == other.d_mol;
  }
  bool operator!=(const IndirectIter &other) const { return !(*this == other); }
};

//! This class is only intended to be used by range-based for loops via
//! RDMol::atomNeighbors() and RDMol::atomBonds()
template <bool ISBOND, bool ISCONST>
class IndirectRange {
  using MolType = std::conditional_t<ISCONST, const RDMol, RDMol>;
  using Iter = IndirectIter<ISBOND, ISCONST>;
  MolType *d_mol;
  const uint32_t *d_begin;
  const uint32_t *d_end;

  friend class RDKit::RDMol;
  IndirectRange(MolType *mol, const uint32_t *beginIndices,
                     const uint32_t *endIndices)
      : d_mol(mol), d_begin(beginIndices), d_end(endIndices) {}

 public:
  IndirectRange(const IndirectRange<ISBOND, ISCONST> &) = default;
  IndirectRange &operator=(const IndirectRange<ISBOND, ISCONST> &) = default;
  Iter begin() const { return Iter(d_mol, d_begin); }
  Iter end() const { return Iter(d_mol, d_end); }
};

}  // namespace Ranges

inline Ranges::IndexRange<false, false> RDMol::atoms() {
  return Ranges::IndexRange<false, false>(this);
}
inline Ranges::IndexRange<false, true> RDMol::atoms() const {
  return Ranges::IndexRange<false, true>(this);
}
inline Ranges::IndexRange<true, false> RDMol::bonds() {
  return Ranges::IndexRange<true, false>(this);
}
inline Ranges::IndexRange<true, true> RDMol::bonds() const {
  return Ranges::IndexRange<true, true>(this);
}
inline Ranges::IndirectRange<false, false> RDMol::atomNeighbors(atomindex_t atomIndex) {
  auto [begin, end] = getAtomNeighbors(atomIndex);
  return Ranges::IndirectRange<false, false>(this, begin, end);
}
inline Ranges::IndirectRange<false, true> RDMol::atomNeighbors(atomindex_t atomIndex) const {
  auto [begin, end] = getAtomNeighbors(atomIndex);
  return Ranges::IndirectRange<false, true>(this, begin, end);
}
inline Ranges::IndirectRange<true, false> RDMol::atomBonds(atomindex_t atomIndex) {
  auto [begin, end] = getAtomBonds(atomIndex);
  return Ranges::IndirectRange<true, false>(this, begin, end);
}
inline Ranges::IndirectRange<true, true> RDMol::atomBonds(atomindex_t atomIndex) const {
  auto [begin, end] = getAtomBonds(atomIndex);
  return Ranges::IndirectRange<true, true>(this, begin, end);
}

}  // namespace RDKit
