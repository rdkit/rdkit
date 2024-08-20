///
//  Copyright (C) 2001-2021 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOLPICKLE_H
#define RD_MOLPICKLE_H

#include <Geometry/point.h>
#include <GraphMol/Atom.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/QueryBond.h>
#include <RDGeneral/StreamOps.h>
#include <boost/utility/binary.hpp>
#include <boost/variant.hpp>
#include <Query/QueryObjects.h>

// Std stuff
#include <iostream>
#include <string>
#include <sstream>
#include <exception>
#ifdef WIN32
#include <ios>
#endif
#include <cstdint>

namespace RDKit {
class ROMol;
class RingInfo;

//! used to indicate exceptions whilst pickling (serializing) molecules
class RDKIT_GRAPHMOL_EXPORT MolPicklerException : public std::exception {
 public:
  MolPicklerException(const char *msg) : _msg(msg) {}
  MolPicklerException(const std::string msg) : _msg(msg) {}
  const char *what() const noexcept override { return _msg.c_str(); }
  ~MolPicklerException() noexcept override = default;

 private:
  std::string _msg;
};

namespace PicklerOps {
typedef enum {
  NoProps = 0,  // no data pickled (default pickling, single-precision coords)
  MolProps = 0x1,  // only public non computed properties
  AtomProps = 0x2,
  BondProps = 0x4,
  QueryAtomData =
      0x2,  // n.b. DEPRECATED and set to AtomProps (does the same work)
  PrivateProps = 0x10,
  ComputedProps = 0x20,
  AllProps = 0x0000FFFF,        // all data pickled
  CoordsAsDouble = 0x00010000,  // save coordinates in double precision
  NoConformers =
      0x00020000  // do not include conformers or associated properties
} PropertyPickleOptions;
}  // namespace PicklerOps

//! handles pickling (serializing) molecules
class RDKIT_GRAPHMOL_EXPORT MolPickler {
 public:
  static const std::int32_t versionMajor;  //!< mark the pickle major version
  static const std::int32_t versionMinor;  //!< mark the pickle minor version
  static const std::int32_t versionPatch;  //!< mark the pickle patch version
  static const std::int32_t endianId;  //!< mark the endian-ness of the pickle

  //! the pickle format is tagged using these tags:
  //! NOTE: if you add to this list, be sure to put new entries AT THE BOTTOM,
  /// otherwise
  //! you will break old pickles.
  typedef enum {
    VERSION = 0,
    BEGINATOM,
    ATOM_INDEX,
    ATOM_NUMBER,
    ATOM_POS,
    ATOM_CHARGE,
    ATOM_NEXPLICIT,
    ATOM_CHIRALTAG,
    ATOM_MASS,
    ATOM_ISAROMATIC,
    ENDATOM,
    BEGINBOND,
    BOND_INDEX,
    BOND_BEGATOMIDX,
    BOND_ENDATOMIDX,
    BOND_TYPE,
    BOND_DIR,
    ENDBOND,
    BEGINPROPS,
    ENDPROPS,
    BEGINSSSR,
    ENDSSSR,
    ENDMOL,
    BEGINCONFS,
    ATOM_MAPNUMBER,
    BEGINQUERY,
    QUERY_VALUE,
    QUERY_ISNEGATED,
    QUERY_NUMCHILDREN,
    QUERY_BOOL,
    QUERY_AND,
    QUERY_OR,
    QUERY_XOR,
    QUERY_EQUALS,
    QUERY_GREATER,
    QUERY_GREATEREQUAL,
    QUERY_LESS,
    QUERY_LESSEQUAL,
    QUERY_RANGE,
    QUERY_SET,
    QUERY_NULL,
    QUERY_ATOMRING,
    QUERY_RECURSIVE,
    ENDQUERY,
    ATOM_DUMMYLABEL,
    BEGIN_ATOM_MONOMER,
    ATOM_PDB_RESIDUE_SERIALNUMBER,
    ATOM_PDB_RESIDUE_ALTLOC,
    ATOM_PDB_RESIDUE_RESIDUENAME,
    ATOM_PDB_RESIDUE_CHAINID,
    ATOM_PDB_RESIDUE_INSERTIONCODE,
    ATOM_PDB_RESIDUE_OCCUPANCY,
    ATOM_PDB_RESIDUE_TEMPFACTOR,
    ATOM_PDB_RESIDUE_ISHETEROATOM,
    ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE,
    ATOM_PDB_RESIDUE_RESIDUENUMBER,
    ATOM_PDB_RESIDUE_SEGMENTNUMBER,
    END_ATOM_MONOMER,
    BEGINATOMPROPS,
    BEGINBONDPROPS,
    BEGINQUERYATOMDATA,
    BEGINSGROUP,
    BEGINSTEREOGROUP,
    BEGINCONFPROPS,
    BEGINCONFS_DOUBLE,
    QUERY_TYPELABEL,
    BEGINSYMMSSSR,
    BEGINFASTFIND,
    BEGINFINDOTHERORUNKNOWN,
    QUERY_PROPERTY,
    // add new entries above here
    INVALID_TAG = 255
  } Tags;

  static unsigned int getDefaultPickleProperties();
  static void setDefaultPickleProperties(unsigned int);

  static const CustomPropHandlerVec &getCustomPropHandlers();
  static void addCustomPropHandler(const CustomPropHandler &handler);

  //! pickles a molecule and sends the results to stream \c ss
  static void pickleMol(const ROMol *mol, std::ostream &ss);
  static void pickleMol(const ROMol *mol, std::ostream &ss,
                        unsigned int propertyFlags);

  static void pickleMol(const ROMol &mol, std::ostream &ss);

  static void pickleMol(const ROMol &mol, std::ostream &ss,
                        unsigned int propertyFlags) {
    MolPickler::pickleMol(&mol, ss, propertyFlags);
  }

  //! pickles a molecule and adds the results to string \c res
  static void pickleMol(const ROMol *mol, std::string &res);
  static void pickleMol(const ROMol *mol, std::string &res,
                        unsigned int propertyFlags);
  static void pickleMol(const ROMol &mol, std::string &res);
  static void pickleMol(const ROMol &mol, std::string &res,
                        unsigned int propertyFlags) {
    MolPickler::pickleMol(&mol, res, propertyFlags);
  }

  //! constructs a molecule from a pickle stored in a string
  static void molFromPickle(const std::string &pickle, ROMol *mol,
                            unsigned int propertyFlags);
  static void molFromPickle(const std::string &pickle, ROMol &mol,
                            unsigned int propertyFlags) {
    MolPickler::molFromPickle(pickle, &mol, propertyFlags);
  }
  static void molFromPickle(const std::string &pickle, ROMol *mol) {
    MolPickler::molFromPickle(pickle, mol, PicklerOps::AllProps);
  }
  static void molFromPickle(const std::string &pickle, ROMol &mol) {
    MolPickler::molFromPickle(pickle, &mol, PicklerOps::AllProps);
  }

  //! constructs a molecule from a pickle stored in a stream
  static void molFromPickle(std::istream &ss, ROMol *mol,
                            unsigned int propertyFlags);
  static void molFromPickle(std::istream &ss, ROMol &mol,
                            unsigned int propertyFlags) {
    MolPickler::molFromPickle(ss, &mol, propertyFlags);
  }
  static void molFromPickle(std::istream &ss, ROMol *mol) {
    MolPickler::molFromPickle(ss, mol, PicklerOps::AllProps);
  }
  static void molFromPickle(std::istream &ss, ROMol &mol) {
    MolPickler::molFromPickle(ss, &mol, PicklerOps::AllProps);
  }

 private:
  //! Pickle nonquery atom data
  static std::int32_t _pickleAtomData(std::ostream &tss, const Atom *atom);
  //! depickle nonquery atom data
  static void _unpickleAtomData(std::istream &tss, Atom *atom, int version);

  static void _pickleQueryAtomData(std::ostream &tss, const Atom *atom);

  //! do the actual work of pickling a molecule
  template <typename T>
  static void _pickle(const ROMol *mol, std::ostream &ss,
                      unsigned int propertyFlags);

  //! do the actual work of pickling an Atom
  template <typename T>
  static void _pickleAtom(std::ostream &ss, const Atom *atom);

  //! do the actual work of pickling a Bond
  template <typename T>
  static void _pickleBond(std::ostream &ss, const Bond *bond,
                          std::map<int, int> &atomIdxMap);

  //! do the actual work of pickling an SSSR structure
  template <typename T>
  static void _pickleSSSR(std::ostream &ss, const RingInfo *ringInfo,
                          std::map<int, int> &atomIdxMap);

  //! do the actual work of pickling a SubstanceGroup
  template <typename T>
  static void _pickleSubstanceGroup(std::ostream &ss,
                                    const SubstanceGroup &sgroup,
                                    std::map<int, int> &atomIdxMap,
                                    std::map<int, int> &bondIdxMap);

  //! do the actual work of pickling Stereo Group data
  template <typename T>
  static void _pickleStereo(std::ostream &ss, std::vector<StereoGroup> groups,
                            std::map<int, int> &atomIdxMap,
                            std::map<int, int> &bondIdxMap);

  //! do the actual work of pickling a Conformer
  template <typename T, typename C>
  static void _pickleConformer(std::ostream &ss, const Conformer *conf);

  //! do the actual work of de-pickling a molecule
  template <typename T>
  static void _depickle(std::istream &ss, ROMol *mol, int version, int numAtoms,
                        unsigned int propertyFlags);

  //! extract atomic data from a pickle and add the resulting Atom to the
  /// molecule
  template <typename T>
  static Atom *_addAtomFromPickle(std::istream &ss, ROMol *mol,
                                  RDGeom::Point3D &pos, int version,
                                  bool directMap = false);

  //! extract bond data from a pickle and add the resulting Bond to the molecule
  template <typename T>
  static Bond *_addBondFromPickle(std::istream &ss, ROMol *mol, int version,
                                  bool directMap = false);

  //! extract ring info from a pickle and add the resulting RingInfo to the
  /// molecule
  template <typename T>
  static void _addRingInfoFromPickle(
      std::istream &ss, ROMol *mol, int version, bool directMap = false,
      FIND_RING_TYPE ringType =
          FIND_RING_TYPE::FIND_RING_TYPE_OTHER_OR_UNKNOWN);

  //! extract a SubstanceGroup from a pickle
  template <typename T>
  static SubstanceGroup _getSubstanceGroupFromPickle(std::istream &ss,
                                                     ROMol *mol, int version);

  template <typename T>
  static void _depickleStereo(std::istream &ss, ROMol *mol, int version);

  //! extract a conformation from a pickle
  template <typename T, typename C>
  static Conformer *_conformerFromPickle(std::istream &ss, int version);

  //! pickle standard properties
  static void _pickleProperties(std::ostream &ss, const RDProps &props,
                                unsigned int pickleFlags);
  //! unpickle standard properties
  static void _unpickleProperties(std::istream &ss, RDProps &props,
                                  int version);

  //! backwards compatibility
  static void _pickleV1(const ROMol *mol, std::ostream &ss);
  //! backwards compatibility
  static void _depickleV1(std::istream &ss, ROMol *mol);
  //! backwards compatibility
  static void _addAtomFromPickleV1(std::istream &ss, ROMol *mol);
  //! backwards compatibility
  static void _addBondFromPickleV1(std::istream &ss, ROMol *mol);
};

namespace PicklerOps {
// clang-format off
using QueryDetails = boost::variant<
    MolPickler::Tags, std::tuple<MolPickler::Tags, int32_t>,
    std::tuple<MolPickler::Tags, int32_t, int32_t>,
    std::tuple<MolPickler::Tags, int32_t, int32_t, int32_t, char>,
    std::tuple<MolPickler::Tags, std::set<int32_t>>,
    std::tuple<MolPickler::Tags, std::string>>;
// clang-format on
template <class T>
QueryDetails getQueryDetails(const Queries::Query<int, T const *, true> *query);

}  // namespace PicklerOps

};  // namespace RDKit

#endif
