//
//  Copyright (C) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <string>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDBoost/Wrap.h>

// #include "seqs.hpp"
// #include "props.hpp"
#include <algorithm>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
void expandQuery(QueryAtom *self, const QueryAtom *other,
                 Queries::CompositeQueryType how, bool maintainOrder) {
  if (!other) {
    throw ValueErrorException("other Atom is null");
  }
  if (other->hasQuery()) {
    const QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
    self->expandQuery(qry->copy(), how, maintainOrder);
  }
}

void setQuery(QueryAtom *self, const QueryAtom *other) {
  if (!other) {
    throw ValueErrorException("other Atom is null");
  }
  if (other->hasQuery()) {
    self->setQuery(other->getQuery()->copy());
  }
}

template <class T>
void AtomSetProp(const Atom *atom, const std::string &key, const T &val) {
  // std::cerr<<"asp: "<<atom<<" " << key<<" - " << val << std::endl;
  atom->setProp<T>(key, val);
}

int AtomHasProp(const Atom *atom, const std::string &key) {
  // std::cerr<<"ahp: "<<atom<<" " << key<< std::endl;
  int res = atom->hasProp(key);
  return res;
}

void AtomClearProp(const Atom *atom, const std::string &key) {
  if (!atom->hasProp(key)) {
    return;
  }
  atom->clearProp(key);
}

// python::tuple AtomGetNeighbors(Atom *atom) {
//   python::list res;
//   for (auto nbr : atom->getOwningMol().atomNeighbors(atom)) {
//     res.append(python::ptr(nbr));
//   }
//   return python::tuple(res);
// }

// python::tuple AtomGetBonds(Atom *atom) {
//   python::list res;
//   for (auto bond : atom->getOwningMol().atomBonds(atom)) {
//     res.append(python::ptr(bond));
//   }
//   return python::tuple(res);
// }

bool AtomIsInRing(const Atom *atom) {
  if (!atom->getOwningMol().getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(atom->getOwningMol());
  }
  return atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) != 0;
}
bool AtomIsInRingSize(const Atom *atom, int size) {
  if (!atom->getOwningMol().getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(atom->getOwningMol());
  }
  return atom->getOwningMol().getRingInfo()->isAtomInRingOfSize(atom->getIdx(),
                                                                size);
}

std::string AtomGetSmarts(const Atom *atom, bool doKekule, bool allHsExplicit,
                          bool isomericSmiles) {
  std::string res;
  if (atom->hasQuery()) {
    res = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
  } else {
    // FIX: this should not be necessary
    res = SmilesWrite::GetAtomSmiles(atom, doKekule, nullptr, allHsExplicit,
                                     isomericSmiles);
  }
  return res;
}

void SetAtomMonomerInfo(Atom *atom, const AtomMonomerInfo *info) {
  if (!info) {
    atom->setMonomerInfo(nullptr);
  } else {
    atom->setMonomerInfo(info->copy());
  }
}

AtomMonomerInfo *AtomGetMonomerInfo(Atom *atom) {
  return atom->getMonomerInfo();
}

void AtomSetPDBResidueInfo(Atom *atom, const AtomMonomerInfo *info) {
  if (!info) {
    // This clears out the monomer info
    atom->setMonomerInfo(nullptr);
    return;
  }

  if (info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    throw ValueErrorException("MonomerInfo is not a PDB Residue");
  }
  atom->setMonomerInfo(info->copy());
}

AtomPDBResidueInfo *AtomGetPDBResidueInfo(Atom *atom) {
  AtomMonomerInfo *res = atom->getMonomerInfo();
  if (!res) {
    return nullptr;
  }
  if (res->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    throw ValueErrorException("MonomerInfo is not a PDB Residue");
  }
  return (AtomPDBResidueInfo *)res;
}

namespace {
int getExplicitValenceHelper(const Atom *atom) {
  RDLog::deprecationWarning("please use GetValence(which=)");
  return atom->getValence(Atom::ValenceType::EXPLICIT);
};
int getImplicitValenceHelper(const Atom *atom) {
  RDLog::deprecationWarning("please use GetValence(getExplicit=False)");
  return atom->getValence(Atom::ValenceType::IMPLICIT);
};
}  // namespace

struct MDLDummy {};
struct DaylightDummy {};

// FIX: is there any reason at all to not just prevent the construction of
// Atoms?
std::string atomClassDoc =
    "The class to store Atoms.\n\
Note that, though it is possible to create one, having an Atom on its own\n\
(i.e not associated with a molecule) is not particularly useful.\n";
struct atom_wrapper {
  static void wrap(nb::module_ &m) {
#if defined(__GNUC__) or defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    nb::class_<Atom>(m, "Atom")
        .def(nb::init<std::string>(), "what"_a)
        .def(nb::init<const Atom &>(), "other"_a)
        .def(nb::init<unsigned int>(), "num"_a,
             "Constructor, takes the atomic number")

        .def("__copy__", &Atom::copy, nb::rv_policy::take_ownership,
             "Create a copy of the atom")

        .def("GetAtomicNum", &Atom::getAtomicNum, "Returns the atomic number.")

        .def("SetAtomicNum", &Atom::setAtomicNum, "newNum"_a,
             "Sets the atomic number, takes an integer value as an argument")

        .def("GetSymbol", &Atom::getSymbol,
             "Returns the atomic symbol (a string)\n")

        .def("GetIdx", &Atom::getIdx,
             "Returns the atom's index (ordering in the molecule)\n")

        .def("GetDegree", &Atom::getDegree,
             "Returns the degree of the atom in the molecule.\n\n"
             "  The degree of an atom is defined to be its number of\n"
             "  directly-bonded neighbors.\n"
             "  The degree is independent of bond orders, but is dependent\n"
             "    on whether or not Hs are explicit in the graph.\n")
        .def("GetTotalDegree", &Atom::getTotalDegree,
             "Returns the degree of the atom in the molecule including Hs.\n\n"
             "  The degree of an atom is defined to be its number of\n"
             "  directly-bonded neighbors.\n"
             "  The degree is independent of bond orders.\n")

        .def("GetTotalNumHs", &Atom::getTotalNumHs,
             "includeNeighbors"_a = false,
             "Returns the total number of Hs (explicit and implicit) on the "
             "atom.\n\n"
             "  ARGUMENTS:\n\n"
             "    - includeNeighbors: (optional) toggles inclusion of "
             "neighboring H atoms in the sum.\n"
             "      Defaults to 0.\n")
        .def("GetNumImplicitHs", &Atom::getNumImplicitHs,
             "Returns the total number of implicit Hs on the atom.\n")
        .def(
            "GetExplicitValence", &getExplicitValenceHelper,
            "DEPRECATED, please use GetValence(Chem.ValenceType,EXPLICIT) instead.\nReturns the explicit valence of the atom.\n")
        .def(
            "GetImplicitValence", &getImplicitValenceHelper,
            "DEPRECATED, please use getValence(Chem.ValenceType,IMPLICIT) instead.\nReturns the number of implicit Hs on the atom.\n")
        .def("GetValence", &Atom::getValence, "which"_a,
             "Returns the valence (explicit or implicit) of the atom.\n")
        .def("GetTotalValence", &Atom::getTotalValence,
             "Returns the total valence (explicit + implicit) of the atom.\n\n")
        .def("HasValenceViolation", &Atom::hasValenceViolation,
             "Returns whether the atom has a valence violation or not.\n\n")

        .def("GetFormalCharge", &Atom::getFormalCharge)
        .def("SetFormalCharge", &Atom::setFormalCharge, "what"_a)

        .def("SetNoImplicit", &Atom::setNoImplicit,
             "what"_a
             "Sets a marker on the atom that *disallows* implicit Hs.\n"
             "  This holds even if the atom would otherwise have implicit Hs "
             "added.\n")
        .def("GetNoImplicit", &Atom::getNoImplicit,
             "Returns whether or not the atom is *allowed* to have implicit "
             "Hs.\n")

        .def("SetNumExplicitHs", &Atom::setNumExplicitHs, "what"_a)
        .def("GetNumExplicitHs", &Atom::getNumExplicitHs)
        .def("SetIsAromatic", &Atom::setIsAromatic, "what"_a)
        .def("GetIsAromatic", &Atom::getIsAromatic)
        .def("GetMass", &Atom::getMass)
        .def("SetIsotope", &Atom::setIsotope, "what"_a)
        .def("GetIsotope", &Atom::getIsotope)
        .def("SetNumRadicalElectrons", &Atom::setNumRadicalElectrons, "num"_a)
        .def("GetNumRadicalElectrons", &Atom::getNumRadicalElectrons)
        .def("GetQueryType", &Atom::getQueryType)

        .def("SetChiralTag", &Atom::setChiralTag, "what"_a)
        .def("InvertChirality", &Atom::invertChirality)
        .def("GetChiralTag", &Atom::getChiralTag)

        .def("SetHybridization", &Atom::setHybridization, "what"_a,
             "Sets the hybridization of the atom.\n"
             "  The argument should be a HybridizationType\n")
        .def("GetHybridization", &Atom::getHybridization,
             "Returns the atom's hybridization.\n")

        .def("HasOwningMol", &Atom::hasOwningMol,
             "Returns whether or not this instance belongs to a molecule.\n")
        .def("GetOwningMol", &Atom::getOwningMol,
             "Returns the Mol that owns this atom.\n",
             nb::rv_policy::reference_internal)
        //    .def("GetNeighbors", AtomGetNeighbors, python::args("self"),
        //         "Returns a read-only sequence of the atom's neighbors\n")

        //    .def("GetBonds", AtomGetBonds, python::args("self"),
        //         "Returns a read-only sequence of the atom's bonds\n")

        .def("Match", (bool(Atom::*)(const Atom *) const) & Atom::Match,
             "other"_a,
             "Returns whether or not this atom matches another Atom.\n\n"
             "  Each Atom (or query Atom) has a query function which is\n"
             "  used for this type of matching.\n\n"
             "  ARGUMENTS:\n"
             "    - other: the other Atom to which to compare\n")

        .def("IsInRingSize", AtomIsInRingSize, "size"_a,
             "Returns whether or not the atom is in a ring of a particular "
             "size.\n\n"
             "  ARGUMENTS:\n"
             "    - size: the ring size to look for\n")

        .def("IsInRing", AtomIsInRing,
             "Returns whether or not the atom is in a ring\n\n")

        .def("HasQuery", &Atom::hasQuery,
             "Returns whether or not the atom has an associated query\n\n")

        .def("DescribeQuery", describeQuery,
             "returns a text description of the query. Primarily intended for "
             "debugging purposes.\n\n")

        .def("GetSmarts", AtomGetSmarts, "doKekule"_a = false,
             "allHsExplicit"_a = false, "isomericSmiles"_a = true,
             "returns the SMARTS (or SMILES) string for an Atom\n\n")
#if 0
        // properties
        .def("SetProp", AtomSetProp<std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n\n")

        .def(
            "GetProp", GetPyProp<Atom>,
            (python::arg("self"), python::arg("key"),
             python::arg("autoConvert") = false),
            "Returns the value of the property.\n\n"
            "  ARGUMENTS:\n"
            "    - key: the name of the property to return (a string).\n\n"
            "    - autoConvert: if True attempt to convert the property into a python object\n\n"
            "  RETURNS: a string\n\n"
            "  NOTE:\n"
            "    - If the property has not been set, a KeyError exception "
            "will be raised.\n",
            boost::python::return_value_policy<return_pyobject_passthrough>())

        .def("SetIntProp", AtomSetProp<int>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a int).\n"
             "    - value: the property value (a int).\n\n")

        .def("SetUnsignedProp", AtomSetProp<unsigned>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an unsigned "
             "integer).\n"
             "    - value: the property value (a int >= 0).\n\n")

        .def("GetIntProp", GetProp<Atom, int>, python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an int).\n\n"
             "  RETURNS: an int\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n",
             boost::python::return_value_policy<return_pyobject_passthrough>())

        .def("GetUnsignedProp", GetProp<Atom, unsigned>,
             python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an unsigned "
             "integer).\n\n"
             "  RETURNS: an integer (Python has no unsigned type)\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n",
             boost::python::return_value_policy<return_pyobject_passthrough>())
        .def("SetDoubleProp", AtomSetProp<double>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a double).\n"
             "    - value: the property value (a double).\n\n")

        .def("GetDoubleProp", GetProp<Atom, double>,
             python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a double).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n",
             boost::python::return_value_policy<return_pyobject_passthrough>())

        .def("SetBoolProp", AtomSetProp<bool>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a bool).\n"
             "    - value: the property value (a bool).\n\n")

        .def("GetBoolProp", GetProp<Atom, bool>, python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a bool).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n",
             boost::python::return_value_policy<return_pyobject_passthrough>())

        .def("SetExplicitBitVectProp", AtomSetProp<ExplicitBitVect>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an "
             "ExplicitBitVect).\n"
             "    - value: the property value (an ExplicitBitVect).\n\n")

        .def("GetExplicitBitVectProp", GetProp<Atom, ExplicitBitVect>,
             python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a "
             "ExplicitBitVect).\n\n"
             "  RETURNS: an ExplicitBitVect \n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n",
             boost::python::return_value_policy<return_pyobject_passthrough>())

        .def("HasProp", AtomHasProp, python::args("self", "key"),
             "Queries a Atom to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")

        .def("ClearProp", AtomClearProp, python::args("self", "key"),
             "Removes a particular property from an Atom (does nothing if not "
             "already set).\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be removed.\n")

        .def("GetPropNames", &Atom::getPropList,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
             "Returns a list of the properties set on the Atom.\n\n")

        .def("GetPropsAsDict", GetPropsAsDict<Atom>,
             (python::arg("self"), python::arg("includePrivate") = true,
              python::arg("includeComputed") = true,
              python::arg("autoConvertStrings") = true),
             getPropsAsDictDocString.c_str())
#endif
        .def("UpdatePropertyCache", &Atom::updatePropertyCache,
             "strict"_a = true,
             "Regenerates computed properties like implicit valence and ring "
             "information.\n\n")

        .def("NeedsUpdatePropertyCache", &Atom::needsUpdatePropertyCache,
             "Returns true or false depending on whether implicit and explicit "
             "valence of the molecule have already been calculated.\n\n")

        .def("ClearPropertyCache", &Atom::clearPropertyCache,
             "Clears implicit and explicit valence information.\n\n")

        .def("GetMonomerInfo", AtomGetMonomerInfo,
             nb::rv_policy::reference_internal,
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("GetPDBResidueInfo", AtomGetPDBResidueInfo,
             nb::rv_policy::reference_internal,
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("SetMonomerInfo", SetAtomMonomerInfo, "info"_a,
             "Sets the atom's MonomerInfo object.\n\n")
        .def("SetPDBResidueInfo", AtomSetPDBResidueInfo, "info"_a,
             "Sets the atom's MonomerInfo object.\n\n")
        .def("GetAtomMapNum", &Atom::getAtomMapNum,
             "Gets the atoms map number, returns 0 if not set")
        .def("SetAtomMapNum", &Atom::setAtomMapNum, "mapno"_a,
             "strict"_a = false,
             "Sets the atoms map number, a value of 0 clears the atom map")
        .doc() = atomClassDoc.c_str();
#if defined(__GNUC__) or defined(__clang__)
#pragma GCC diagnostic pop
#endif

    nb::enum_<Atom::HybridizationType>(m, "HybridizationType")
        .value("UNSPECIFIED", Atom::UNSPECIFIED)
        .value("S", Atom::S)
        .value("SP", Atom::SP)
        .value("SP2", Atom::SP2)
        .value("SP3", Atom::SP3)
        .value("SP2D", Atom::SP2D)
        .value("SP3D", Atom::SP3D)
        .value("SP3D2", Atom::SP3D2)
        .value("OTHER", Atom::OTHER);
    nb::enum_<Atom::ChiralType>(m, "ChiralType")
        .value("CHI_UNSPECIFIED", Atom::CHI_UNSPECIFIED)
        .value("CHI_TETRAHEDRAL_CW", Atom::CHI_TETRAHEDRAL_CW)
        .value("CHI_TETRAHEDRAL_CCW", Atom::CHI_TETRAHEDRAL_CCW)
        .value("CHI_OTHER", Atom::CHI_OTHER)
        .value("CHI_TETRAHEDRAL", Atom::CHI_TETRAHEDRAL)
        .value("CHI_ALLENE", Atom::CHI_ALLENE)
        .value("CHI_SQUAREPLANAR", Atom::CHI_SQUAREPLANAR)
        .value("CHI_TRIGONALBIPYRAMIDAL", Atom::CHI_TRIGONALBIPYRAMIDAL)
        .value("CHI_OCTAHEDRAL", Atom::CHI_OCTAHEDRAL)
        .export_values();

    nb::enum_<Atom::ValenceType>(m, "ValenceType")
        .value("IMPLICIT", Atom::ValenceType::IMPLICIT)
        .value("EXPLICIT", Atom::ValenceType::EXPLICIT)
        .export_values();

    nb::enum_<Queries::CompositeQueryType>(m, "CompositeQueryType")
        .value("COMPOSITE_AND", Queries::COMPOSITE_AND)
        .value("COMPOSITE_OR", Queries::COMPOSITE_OR)
        .value("COMPOSITE_XOR", Queries::COMPOSITE_XOR)
        .export_values();

    atomClassDoc =
        "The class to store QueryAtoms.\n\
These cannot currently be constructed directly from Python\n";
    nb::class_<QueryAtom, Atom>(m, "QueryAtom")
        .def("ExpandQuery", expandQuery, "other"_a,
             "how"_a = Queries::COMPOSITE_AND, "maintainOrder"_a = true,
             "combines the query from other with ours")
        .def("SetQuery", setQuery, "other"_a,
             "Replace our query with a copy of the other query")
        .doc() = atomClassDoc.c_str();

    m.def(
        "GetAtomRLabel", getAtomRLabel, "atom"_a,
        "Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)");
    m.def("SetAtomRLabel", setAtomRLabel, "atom"_a, "rlabel"_a,
          "Sets the atom's MDL RLabel (this is an integer from 0 to "
          "99).\nSetting to 0 clears the rlabel.");

    m.def("GetAtomAlias", getAtomAlias, "atom"_a,
          "Returns the atom's MDL alias text");
    m.def("SetAtomAlias", setAtomAlias, "atom"_a, "rlabel"_a,
          "Sets the atom's MDL alias text.\nSetting to an empty string "
          "clears the alias.");
    m.def("GetAtomValue", getAtomValue, "atom"_a,
          "Returns the atom's MDL alias text");
    m.def("SetAtomValue", setAtomValue, "atom"_a, "rlabel"_a,
          "Sets the atom's MDL alias text.\nSetting to an empty string "
          "clears the alias.");

    m.def("GetSupplementalSmilesLabel", getSupplementalSmilesLabel, "atom"_a,
          "Gets the supplemental smiles label on an atom, returns an "
          "empty string if not present.");
    m.def("SetSupplementalSmilesLabel", setSupplementalSmilesLabel, "atom"_a,
          "label"_a,
          "Sets a supplemental label on an atom that is written to the smiles "
          "string.\n\n"
          ">>> m = Chem.MolFromSmiles(\"C\")\n"
          ">>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')\n"
          ">>> Chem.MolToSmiles(m)\n"
          "'C<xxx>'\n");
    m.def("GetNumPiElectrons", numPiElectrons, "atom"_a,
          "Returns the number of electrons an atom is using for pi bonding");
  }
};
}  // namespace RDKit
void wrap_atom(nb::module_ &m) { RDKit::atom_wrapper::wrap(m); }
