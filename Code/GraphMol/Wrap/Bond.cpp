//
//  Copyright (C) 2003-2017 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <string>
#include "props.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit {

void expandQuery(QueryBond *self, const QueryBond *other,
                 Queries::CompositeQueryType how, bool maintainOrder) {
  if (other->hasQuery()) {
    const QueryBond::QUERYBOND_QUERY *qry = other->getQuery();
    self->expandQuery(qry->copy(), how, maintainOrder);
  }
}

void setQuery(QueryBond *self, const QueryBond *other) {
  if (other->hasQuery()) {
    self->setQuery(other->getQuery()->copy());
  }
}

int BondHasProp(const Bond *bond, const char *key) {
  int res = bond->hasProp(key);
  return res;
}

template <class T>
void BondSetProp(const Bond *bond, const char *key, const T &val) {
  bond->setProp<T>(key, val);
}
void BondClearProp(const Bond *bond, const char *key) {
  if (!bond->hasProp(key)) {
    return;
  }
  bond->clearProp(key);
}

bool BondIsInRing(const Bond *bond) {
  if (!bond->getOwningMol().getRingInfo()->isInitialized()) {
    MolOps::findSSSR(bond->getOwningMol());
  }
  return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx()) != 0;
}

bool BondIsInRingSize(const Bond *bond, int size) {
  if (!bond->getOwningMol().getRingInfo()->isInitialized()) {
    MolOps::findSSSR(bond->getOwningMol());
  }
  return bond->getOwningMol().getRingInfo()->isBondInRingOfSize(bond->getIdx(),
                                                                size);
  return false;
}

INT_VECT getBondStereoAtoms(const Bond *bond) { return bond->getStereoAtoms(); }

std::string BondGetSmarts(const Bond *bond, bool allBondsExplicit) {
  std::string res;
  if (bond->hasQuery()) {
    res = SmartsWrite::GetBondSmarts(static_cast<const QueryBond *>(bond));
  } else {
    res = SmilesWrite::GetBondSmiles(bond, -1, false, allBondsExplicit);
  }
  return res;
}

std::string bondClassDoc =
    "The class to store Bonds.\n\
Note: unlike Atoms, is it currently impossible to construct Bonds from\n\
Python.\n";
struct bond_wrapper {
  static void wrap() {
    python::class_<Bond>("Bond", bondClassDoc.c_str(), python::no_init)

        .def("HasOwningMol", &Bond::hasOwningMol,
             "Returns whether or not this instance belongs to a molecule.\n")
        .def("GetOwningMol", &Bond::getOwningMol,
             "Returns the Mol that owns this bond.\n",
             python::return_value_policy<python::reference_existing_object>())

        .def("GetBondType", &Bond::getBondType,
             "Returns the type of the bond as a BondType\n")
        .def("SetBondType", &Bond::setBondType,
             "Set the type of the bond as a BondType\n")
        .def("GetBondTypeAsDouble", &Bond::getBondTypeAsDouble,
             "Returns the type of the bond as a double (i.e. 1.0 for SINGLE, "
             "1.5 for AROMATIC, 2.0 for DOUBLE)\n")

        .def("GetBondDir", &Bond::getBondDir,
             "Returns the type of the bond as a BondDir\n")
        .def("SetBondDir", &Bond::setBondDir,
             "Set the type of the bond as a BondDir\n")

        .def("GetStereo", &Bond::getStereo,
             "Returns the stereo configuration of the bond as a BondStereo\n")
        .def("SetStereo", &Bond::setStereo,
             "Set the stereo configuration of the bond as a BondStereo\n")
        .def("GetStereoAtoms", getBondStereoAtoms,
             "Returns the indices of the atoms setting this bond's "
             "stereochemistry.\n")
        .def("SetStereoAtoms", &Bond::setStereoAtoms,
             "Set the indices of the atoms setting this bond's "
             "stereochemistry.\n")

        .def("GetValenceContrib",
             (double (Bond::*)(const Atom *) const) & Bond::getValenceContrib,
             "Returns the contribution of the bond to the valence of an "
             "Atom.\n\n"
             "  ARGUMENTS:\n\n"
             "    - atom: the Atom to consider.\n")

        .def("GetIsAromatic", &Bond::getIsAromatic)
        .def("SetIsAromatic", &Bond::setIsAromatic)

        .def("GetIsConjugated", &Bond::getIsConjugated,
             "Returns whether or not the bond is considered to be conjugated.")
        .def("SetIsConjugated", &Bond::setIsConjugated)

        .def("GetIdx", &Bond::getIdx,
             "Returns the bond's index (ordering in the molecule)\n")

        .def("GetBeginAtomIdx", &Bond::getBeginAtomIdx,
             "Returns the index of the bond's first atom.\n")
        .def("GetEndAtomIdx", &Bond::getEndAtomIdx,
             "Returns the index of the bond's first atom.\n")
        .def("GetOtherAtomIdx", &Bond::getOtherAtomIdx,
             "Given the index of one of the bond's atoms, returns the\n"
             "index of the other.\n")

        .def("GetBeginAtom", &Bond::getBeginAtom,
             python::return_value_policy<python::reference_existing_object>(),
             "Returns the bond's first atom.\n")
        .def("GetEndAtom", &Bond::getEndAtom,
             python::return_value_policy<python::reference_existing_object>(),
             "Returns the bond's second atom.\n")
        .def("GetOtherAtom", &Bond::getOtherAtom,
             python::return_value_policy<python::reference_existing_object>(),
             "Given one of the bond's atoms, returns the other one.\n")

        // FIX: query stuff
        .def("Match", (bool (Bond::*)(const Bond *) const) & Bond::Match,
             "Returns whether or not this bond matches another Bond.\n\n"
             "  Each Bond (or query Bond) has a query function which is\n"
             "  used for this type of matching.\n\n"
             "  ARGUMENTS:\n"
             "    - other: the other Bond to which to compare\n")

        .def("IsInRingSize", BondIsInRingSize,
             "Returns whether or not the bond is in a ring of a particular "
             "size.\n\n"
             "  ARGUMENTS:\n"
             "    - size: the ring size to look for\n")
        .def("IsInRing", BondIsInRing,
             "Returns whether or not the bond is in a ring of any size.\n\n")

        .def("HasQuery", &Bond::hasQuery,
             "Returns whether or not the bond has an associated query\n\n")

        .def("DescribeQuery", describeQuery,
             "returns a text description of the query. Primarily intended for "
             "debugging purposes.\n\n")

        .def("GetSmarts", BondGetSmarts,
             (python::arg("bond"), python::arg("allBondsExplicit") = false),
             "returns the SMARTS (or SMILES) string for a Bond")

        // properties
        .def("SetProp", BondSetProp<std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets a bond property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n\n")

        .def("GetProp", GetProp<Bond, std::string>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a string\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("SetIntProp", BondSetProp<int>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets a bond property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (an int).\n\n")

        .def("SetUnsignedProp", BondSetProp<unsigned int>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets a bond property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (an int >= 0).\n\n")

        .def("GetIntProp", GetProp<Bond, int>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an int).\n\n"
             "  RETURNS: an int\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("GetUnsignedProp", GetProp<Bond, unsigned int>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an unsigned "
             "integer).\n\n"
             "  RETURNS: an int (Python has no unsigned type)\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("SetDoubleProp", BondSetProp<double>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets a bond property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a double).\n\n")

        .def("GetDoubleProp", GetProp<Bond, double>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a double).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("SetBoolProp", BondSetProp<bool>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets a bond property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a boolean).\n\n")

        .def("GetBoolProp", GetProp<Bond, bool>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a boolean).\n\n"
             "  RETURNS: a boolean\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("HasProp", BondHasProp,
             "Queries a Bond to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")

        .def("ClearProp", BondClearProp,
             "Removes a particular property from an Bond (does nothing if not "
             "already set).\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be removed.\n")

        .def("GetPropNames", &Bond::getPropList,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
             "Returns a list of the properties set on the Bond.\n\n")

        .def("GetPropsAsDict", GetPropsAsDict<Bond>,
             (python::arg("self"), python::arg("includePrivate") = true,
              python::arg("includeComputed") = true),
             "Returns a dictionary of the properties set on the Bond.\n"
             " n.b. some properties cannot be converted to python types.\n")

        ;

    python::enum_<Bond::BondType>("BondType")
        .value("UNSPECIFIED", Bond::UNSPECIFIED)
        .value("SINGLE", Bond::SINGLE)
        .value("DOUBLE", Bond::DOUBLE)
        .value("TRIPLE", Bond::TRIPLE)
        .value("QUADRUPLE", Bond::QUADRUPLE)
        .value("QUINTUPLE", Bond::QUINTUPLE)
        .value("HEXTUPLE", Bond::HEXTUPLE)
        .value("ONEANDAHALF", Bond::ONEANDAHALF)
        .value("TWOANDAHALF", Bond::TWOANDAHALF)
        .value("THREEANDAHALF", Bond::THREEANDAHALF)
        .value("FOURANDAHALF", Bond::FOURANDAHALF)
        .value("FIVEANDAHALF", Bond::FIVEANDAHALF)
        .value("AROMATIC", Bond::AROMATIC)
        .value("IONIC", Bond::IONIC)
        .value("HYDROGEN", Bond::HYDROGEN)
        .value("THREECENTER", Bond::THREECENTER)
        .value("DATIVEONE", Bond::DATIVEONE)
        .value("DATIVE", Bond::DATIVE)
        .value("DATIVEL", Bond::DATIVEL)
        .value("DATIVER", Bond::DATIVER)
        .value("OTHER", Bond::OTHER)
        .value("ZERO", Bond::ZERO);
    python::enum_<Bond::BondDir>("BondDir")
        .value("NONE", Bond::NONE)
        .value("BEGINWEDGE", Bond::BEGINWEDGE)
        .value("BEGINDASH", Bond::BEGINDASH)
        .value("ENDDOWNRIGHT", Bond::ENDDOWNRIGHT)
        .value("ENDUPRIGHT", Bond::ENDUPRIGHT)
        .value("EITHERDOUBLE", Bond::EITHERDOUBLE)
        .value("UNKNOWN", Bond::UNKNOWN);
    python::enum_<Bond::BondStereo>("BondStereo")
        .value("STEREONONE", Bond::STEREONONE)
        .value("STEREOANY", Bond::STEREOANY)
        .value("STEREOZ", Bond::STEREOZ)
        .value("STEREOE", Bond::STEREOE)
        .value("STEREOCIS", Bond::STEREOCIS)
        .value("STEREOTRANS", Bond::STEREOTRANS);

    bondClassDoc =
        "The class to store QueryBonds.\n\
These cannot currently be constructed directly from Python\n";
    python::class_<QueryBond, python::bases<Bond>>(
        "QueryBond", bondClassDoc.c_str(), python::no_init)
        .def("ExpandQuery", expandQuery,
             (python::arg("self"), python::arg("other"),
              python::arg("how") = Queries::COMPOSITE_AND,
              python::arg("maintainOrder") = true),
             "combines the query from other with ours")
        .def("SetQuery", setQuery, (python::arg("self"), python::arg("other")),
             "Replace our query with a copy of the other query");
  };
};
}  // namespace RDKit
void wrap_bond() { RDKit::bond_wrapper::wrap(); }
