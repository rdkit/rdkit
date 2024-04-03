//
//  Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDBoost/Wrap.h>

#include "seqs.hpp"
#include "props.hpp"
#include <algorithm>

namespace python = boost::python;
namespace RDKit {
void expandQuery(QueryAtom *self, const QueryAtom *other,
                 Queries::CompositeQueryType how, bool maintainOrder) {
  if (other->hasQuery()) {
    const QueryAtom::QUERYATOM_QUERY *qry = other->getQuery();
    self->expandQuery(qry->copy(), how, maintainOrder);
  }
}

void setQuery(QueryAtom *self, const QueryAtom *other) {
  if (other->hasQuery()) {
    self->setQuery(other->getQuery()->copy());
  }
}

template <class T>
void AtomSetProp(const Atom *atom, const char *key, const T &val) {
  // std::cerr<<"asp: "<<atom<<" " << key<<" - " << val << std::endl;
  atom->setProp<T>(key, val);
}

int AtomHasProp(const Atom *atom, const char *key) {
  // std::cerr<<"ahp: "<<atom<<" " << key<< std::endl;
  int res = atom->hasProp(key);
  return res;
}

void AtomClearProp(const Atom *atom, const char *key) {
  if (!atom->hasProp(key)) {
    return;
  }
  atom->clearProp(key);
}

python::tuple AtomGetNeighbors(Atom *atom) {
  python::list res;
  for (auto nbr : atom->getOwningMol().atomNeighbors(atom)) {
    res.append(python::ptr(nbr));
  }
  return python::tuple(res);
}

python::tuple AtomGetBonds(Atom *atom) {
  python::list res;
  for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    res.append(python::ptr(bond));
  }
  return python::tuple(res);
}

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
  atom->setMonomerInfo(info->copy());
}

AtomMonomerInfo *AtomGetMonomerInfo(Atom *atom) {
  return atom->getMonomerInfo();
}

void AtomSetPDBResidueInfo(Atom *atom, const AtomMonomerInfo *info) {
  if (info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    throw_value_error("MonomerInfo is not a PDB Residue");
  }
  atom->setMonomerInfo(info->copy());
}

AtomPDBResidueInfo *AtomGetPDBResidueInfo(Atom *atom) {
  AtomMonomerInfo *res = atom->getMonomerInfo();
  if (!res) {
    return nullptr;
  }
  if (res->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    throw_value_error("MonomerInfo is not a PDB Residue");
  }
  return (AtomPDBResidueInfo *)res;
}

struct MDLDummy {};
struct DaylightDummy {};

// FIX: is there any reason at all to not just prevent the construction of
// Atoms?
std::string atomClassDoc =
    "The class to store Atoms.\n\
Note that, though it is possible to create one, having an Atom on its own\n\
(i.e not associated with a molecule) is not particularly useful.\n";
struct atom_wrapper {
  static void wrap() {
    python::class_<Atom>(
        "Atom", atomClassDoc.c_str(),
        python::init<std::string>(python::args("self", "what")))

        .def(python::init<const Atom &>(python::args("self", "other")))
        .def(python::init<unsigned int>(python::args("self", "num"),
                                        "Constructor, takes the atomic number"))

        .def("__copy__", &Atom::copy,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             python::args("self"), "Create a copy of the atom")

        .def("GetAtomicNum", &Atom::getAtomicNum, python::args("self"),
             "Returns the atomic number.")

        .def("SetAtomicNum", &Atom::setAtomicNum,
             python::args("self", "newNum"),
             "Sets the atomic number, takes an integer value as an argument")

        .def("GetSymbol", &Atom::getSymbol, python::args("self"),
             "Returns the atomic symbol (a string)\n")

        .def("GetIdx", &Atom::getIdx, python::args("self"),
             "Returns the atom's index (ordering in the molecule)\n")

        .def("GetDegree", &Atom::getDegree, python::args("self"),
             "Returns the degree of the atom in the molecule.\n\n"
             "  The degree of an atom is defined to be its number of\n"
             "  directly-bonded neighbors.\n"
             "  The degree is independent of bond orders, but is dependent\n"
             "    on whether or not Hs are explicit in the graph.\n")
        .def("GetTotalDegree", &Atom::getTotalDegree, python::args("self"),
             "Returns the degree of the atom in the molecule including Hs.\n\n"
             "  The degree of an atom is defined to be its number of\n"
             "  directly-bonded neighbors.\n"
             "  The degree is independent of bond orders.\n")

        .def("GetTotalNumHs", &Atom::getTotalNumHs,
             (python::arg("self"), python::arg("includeNeighbors") = false),
             "Returns the total number of Hs (explicit and implicit) on the "
             "atom.\n\n"
             "  ARGUMENTS:\n\n"
             "    - includeNeighbors: (optional) toggles inclusion of "
             "neighboring H atoms in the sum.\n"
             "      Defaults to 0.\n")
        .def("GetNumImplicitHs", &Atom::getNumImplicitHs, python::args("self"),
             "Returns the total number of implicit Hs on the atom.\n")

        .def("GetExplicitValence", &Atom::getExplicitValence,
             python::args("self"),
             "Returns the explicit valence of the atom.\n")
        .def("GetImplicitValence", &Atom::getImplicitValence,
             python::args("self"),
             "Returns the number of implicit Hs on the atom.\n")
        .def("GetTotalValence", &Atom::getTotalValence, python::args("self"),
             "Returns the total valence (explicit + implicit) of the atom.\n\n")
        .def("HasValenceViolation", &Atom::hasValenceViolation,
             "Returns whether the atom has a valence violation or not.\n\n")

        .def("GetFormalCharge", &Atom::getFormalCharge, python::args("self"))
        .def("SetFormalCharge", &Atom::setFormalCharge,
             python::args("self", "what"))

        .def("SetNoImplicit", &Atom::setNoImplicit,
             python::args("self", "what"),
             "Sets a marker on the atom that *disallows* implicit Hs.\n"
             "  This holds even if the atom would otherwise have implicit Hs "
             "added.\n")
        .def("GetNoImplicit", &Atom::getNoImplicit, python::args("self"),
             "Returns whether or not the atom is *allowed* to have implicit "
             "Hs.\n")

        .def("SetNumExplicitHs", &Atom::setNumExplicitHs,
             python::args("self", "what"))
        .def("GetNumExplicitHs", &Atom::getNumExplicitHs, python::args("self"))
        .def("SetIsAromatic", &Atom::setIsAromatic,
             python::args("self", "what"))
        .def("GetIsAromatic", &Atom::getIsAromatic, python::args("self"))
        .def("GetMass", &Atom::getMass, python::args("self"))
        .def("SetIsotope", &Atom::setIsotope, python::args("self", "what"))
        .def("GetIsotope", &Atom::getIsotope, python::args("self"))
        .def("SetNumRadicalElectrons", &Atom::setNumRadicalElectrons,
             python::args("self", "num"))
        .def("GetNumRadicalElectrons", &Atom::getNumRadicalElectrons,
             python::args("self"))
        .def("GetQueryType", &Atom::getQueryType, python::args("self"))

        .def("SetChiralTag", &Atom::setChiralTag, python::args("self", "what"))
        .def("InvertChirality", &Atom::invertChirality, python::args("self"))
        .def("GetChiralTag", &Atom::getChiralTag, python::args("self"))

        .def("SetHybridization", &Atom::setHybridization,
             python::args("self", "what"),
             "Sets the hybridization of the atom.\n"
             "  The argument should be a HybridizationType\n")
        .def("GetHybridization", &Atom::getHybridization, python::args("self"),
             "Returns the atom's hybridization.\n")

        .def("HasOwningMol", &Atom::hasOwningMol, python::args("self"),
             "Returns whether or not this instance belongs to a molecule.\n")
        .def("GetOwningMol", &Atom::getOwningMol,
             "Returns the Mol that owns this atom.\n",
             python::return_internal_reference<>(), python::args("self"))

        .def("GetNeighbors", AtomGetNeighbors, python::args("self"),
             "Returns a read-only sequence of the atom's neighbors\n")

        .def("GetBonds", AtomGetBonds, python::args("self"),
             "Returns a read-only sequence of the atom's bonds\n")

        .def("Match", (bool(Atom::*)(const Atom *) const) & Atom::Match,
             python::args("self", "what"),
             "Returns whether or not this atom matches another Atom.\n\n"
             "  Each Atom (or query Atom) has a query function which is\n"
             "  used for this type of matching.\n\n"
             "  ARGUMENTS:\n"
             "    - other: the other Atom to which to compare\n")

        .def("IsInRingSize", AtomIsInRingSize, python::args("self", "size"),
             "Returns whether or not the atom is in a ring of a particular "
             "size.\n\n"
             "  ARGUMENTS:\n"
             "    - size: the ring size to look for\n")

        .def("IsInRing", AtomIsInRing, python::args("self"),
             "Returns whether or not the atom is in a ring\n\n")

        .def("HasQuery", &Atom::hasQuery, python::args("self"),
             "Returns whether or not the atom has an associated query\n\n")

        .def("DescribeQuery", describeQuery, python::args("self"),
             "returns a text description of the query. Primarily intended for "
             "debugging purposes.\n\n")

        .def("GetSmarts", AtomGetSmarts,
             (python::arg("self"), python::arg("doKekule") = false,
              python::arg("allHsExplicit") = false,
              python::arg("isomericSmiles") = true),
             "returns the SMARTS (or SMILES) string for an Atom\n\n")

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
            "will be raised.\n")

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
             "will be raised.\n")

        .def("GetUnsignedProp", GetProp<Atom, unsigned>,
             python::args("self", "key"),
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an unsigned "
             "integer).\n\n"
             "  RETURNS: an integer (Python has no unsigned type)\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

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
             "will be raised.\n")

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
             "will be raised.\n")

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
             "will be raised.\n")

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
             "Returns a dictionary of the properties set on the Atom.\n"
             " n.b. some properties cannot be converted to python types.\n")

        .def("UpdatePropertyCache", &Atom::updatePropertyCache,
             (python::arg("self"), python::arg("strict") = true),
             "Regenerates computed properties like implicit valence and ring "
             "information.\n\n")

        .def("NeedsUpdatePropertyCache", &Atom::needsUpdatePropertyCache,
             (python::arg("self")),
             "Returns true or false depending on whether implicit and explicit "
             "valence of the molecule have already been calculated.\n\n")

        .def("GetMonomerInfo", AtomGetMonomerInfo,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>(),
             python::args("self"),
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("GetPDBResidueInfo", AtomGetPDBResidueInfo,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>(),
             python::args("self"),
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("SetMonomerInfo", SetAtomMonomerInfo, python::args("self", "info"),
             "Sets the atom's MonomerInfo object.\n\n")
        .def("SetPDBResidueInfo", AtomSetPDBResidueInfo,
             python::args("self", "info"),
             "Sets the atom's MonomerInfo object.\n\n")
        .def("GetAtomMapNum", &Atom::getAtomMapNum, python::args("self"),
             "Gets the atoms map number, returns 0 if not set")
        .def("SetAtomMapNum", &Atom::setAtomMapNum,
             (python::arg("self"), python::arg("mapno"),
              python::arg("strict") = false),
             "Sets the atoms map number, a value of 0 clears the atom map");

    python::enum_<Atom::HybridizationType>("HybridizationType")
        .value("UNSPECIFIED", Atom::UNSPECIFIED)
        .value("S", Atom::S)
        .value("SP", Atom::SP)
        .value("SP2", Atom::SP2)
        .value("SP3", Atom::SP3)
        .value("SP2D", Atom::SP2D)
        .value("SP3D", Atom::SP3D)
        .value("SP3D2", Atom::SP3D2)
        .value("OTHER", Atom::OTHER);
    python::enum_<Atom::ChiralType>("ChiralType")
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
    ;

    python::enum_<Queries::CompositeQueryType>("CompositeQueryType")
        .value("COMPOSITE_AND", Queries::COMPOSITE_AND)
        .value("COMPOSITE_OR", Queries::COMPOSITE_OR)
        .value("COMPOSITE_XOR", Queries::COMPOSITE_XOR)
        .export_values();
    ;

    atomClassDoc =
        "The class to store QueryAtoms.\n\
These cannot currently be constructed directly from Python\n";
    python::class_<QueryAtom, python::bases<Atom>>(
        "QueryAtom", atomClassDoc.c_str(), python::no_init)
        .def("ExpandQuery", expandQuery,
             (python::arg("self"), python::arg("other"),
              python::arg("how") = Queries::COMPOSITE_AND,
              python::arg("maintainOrder") = true),
             "combines the query from other with ours")
        .def("SetQuery", setQuery, (python::arg("self"), python::arg("other")),
             "Replace our query with a copy of the other query");

    python::def(
        "GetAtomRLabel", getAtomRLabel, (python::arg("atom")),
        "Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)");
    python::def("SetAtomRLabel", setAtomRLabel,
                (python::arg("atom"), python::arg("rlabel")),
                "Sets the atom's MDL RLabel (this is an integer from 0 to "
                "99).\nSetting to 0 clears the rlabel.");

    python::def("GetAtomAlias", getAtomAlias, (python::arg("atom")),
                "Returns the atom's MDL alias text");
    python::def("SetAtomAlias", setAtomAlias,
                (python::arg("atom"), python::arg("rlabel")),
                "Sets the atom's MDL alias text.\nSetting to an empty string "
                "clears the alias.");
    python::def("GetAtomValue", getAtomValue, (python::arg("atom")),
                "Returns the atom's MDL alias text");
    python::def("SetAtomValue", setAtomValue,
                (python::arg("atom"), python::arg("rlabel")),
                "Sets the atom's MDL alias text.\nSetting to an empty string "
                "clears the alias.");

    python::def("GetSupplementalSmilesLabel", getSupplementalSmilesLabel,
                (python::arg("atom")),
                "Gets the supplemental smiles label on an atom, returns an "
                "empty string if not present.");
    python::def(
        "SetSupplementalSmilesLabel", setSupplementalSmilesLabel,
        (python::arg("atom"), python::arg("label")),
        "Sets a supplemental label on an atom that is written to the smiles "
        "string.\n\n"
        ">>> m = Chem.MolFromSmiles(\"C\")\n"
        ">>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')\n"
        ">>> Chem.MolToSmiles(m)\n"
        "'C<xxx>'\n");
    python::def(
        "GetNumPiElectrons", numPiElectrons, (python::arg("atom")),
        "Returns the number of electrons an atom is using for pi bonding");
  }
};
}  // namespace RDKit
void wrap_atom() { RDKit::atom_wrapper::wrap(); }
