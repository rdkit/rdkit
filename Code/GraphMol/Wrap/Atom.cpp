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
  const ROMol *parent = &atom->getOwningMol();
  ROMol::ADJ_ITER begin, end;
  boost::tie(begin, end) = parent->getAtomNeighbors(atom);
  while (begin != end) {
    res.append(python::ptr(parent->getAtomWithIdx(*begin)));
    begin++;
  }
  return python::tuple(res);
}

python::tuple AtomGetBonds(Atom *atom) {
  python::list res;
  const ROMol *parent = &atom->getOwningMol();
  ROMol::OEDGE_ITER begin, end;
  boost::tie(begin, end) = parent->getAtomBonds(atom);
  while (begin != end) {
    Bond *tmpB = (*parent)[*begin].get();
    res.append(python::ptr(tmpB));
    begin++;
  }
  return python::tuple(res);
}

bool AtomIsInRing(const Atom *atom) {
  if (!atom->getOwningMol().getRingInfo()->isInitialized()) {
    MolOps::findSSSR(atom->getOwningMol());
  }
  return atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) != 0;
}
bool AtomIsInRingSize(const Atom *atom, int size) {
  if (!atom->getOwningMol().getRingInfo()->isInitialized()) {
    MolOps::findSSSR(atom->getOwningMol());
  }
  return atom->getOwningMol().getRingInfo()->isAtomInRingOfSize(atom->getIdx(),
                                                                size);
}

std::string AtomGetSmarts(const Atom *atom) {
  std::string res;
  if (atom->hasQuery()) {
    res = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
  } else {
    res = SmilesWrite::GetAtomSmiles(atom);
  }
  return res;
}

void SetAtomMonomerInfo(Atom *atom, const AtomMonomerInfo *info) {
  atom->setMonomerInfo(info->copy());
}

AtomMonomerInfo *AtomGetMonomerInfo(Atom *atom) {
  return atom->getMonomerInfo();
}
AtomPDBResidueInfo *AtomGetPDBResidueInfo(Atom *atom) {
  AtomMonomerInfo *res = atom->getMonomerInfo();
  if (!res) return NULL;
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
    python::class_<Atom>("Atom", atomClassDoc.c_str(),
                         python::init<std::string>())

        .def(python::init<unsigned int>(
            "Constructor, takes either an int (atomic number) or a string "
            "(atomic symbol).\n"))

        .def("GetAtomicNum", &Atom::getAtomicNum, "Returns the atomic number.")

        .def("SetAtomicNum", &Atom::setAtomicNum,
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
             (python::arg("self"), python::arg("includeNeighbors") = false),
             "Returns the total number of Hs (explicit and implicit) on the "
             "atom.\n\n"
             "  ARGUMENTS:\n\n"
             "    - includeNeighbors: (optional) toggles inclusion of "
             "neighboring H atoms in the sum.\n"
             "      Defaults to 0.\n")
        .def("GetNumImplicitHs", &Atom::getNumImplicitHs,
             "Returns the total number of implicit Hs on the atom.\n")

        .def("GetExplicitValence", &Atom::getExplicitValence,
             "Returns the explicit valence of the atom.\n")
        .def("GetImplicitValence", &Atom::getImplicitValence,
             "Returns the number of implicit Hs on the atom.\n")
        .def("GetTotalValence", &Atom::getTotalValence,
             "Returns the total valence (explicit + implicit) of the atom.\n\n")

        .def("GetFormalCharge", &Atom::getFormalCharge)
        .def("SetFormalCharge", &Atom::setFormalCharge)

        .def("SetNoImplicit", &Atom::setNoImplicit,
             "Sets a marker on the atom that *disallows* implicit Hs.\n"
             "  This holds even if the atom would otherwise have implicit Hs "
             "added.\n")
        .def("GetNoImplicit", &Atom::getNoImplicit,
             "Returns whether or not the atom is *allowed* to have implicit "
             "Hs.\n")

        .def("SetNumExplicitHs", &Atom::setNumExplicitHs)
        .def("GetNumExplicitHs", &Atom::getNumExplicitHs)
        .def("SetIsAromatic", &Atom::setIsAromatic)
        .def("GetIsAromatic", &Atom::getIsAromatic)
        .def("GetMass", &Atom::getMass)
        .def("SetIsotope", &Atom::setIsotope)
        .def("GetIsotope", &Atom::getIsotope)
        .def("SetNumRadicalElectrons", &Atom::setNumRadicalElectrons)
        .def("GetNumRadicalElectrons", &Atom::getNumRadicalElectrons)

        .def("SetChiralTag", &Atom::setChiralTag)
        .def("InvertChirality", &Atom::invertChirality)
        .def("GetChiralTag", &Atom::getChiralTag)

        .def("SetHybridization", &Atom::setHybridization,
             "Sets the hybridization of the atom.\n"
             "  The argument should be a HybridizationType\n")
        .def("GetHybridization", &Atom::getHybridization,
             "Returns the atom's hybridization.\n")

        .def("GetOwningMol", &Atom::getOwningMol,
             "Returns the Mol that owns this atom.\n",
             python::return_internal_reference<>())

        .def("GetNeighbors", AtomGetNeighbors,
             "Returns a read-only sequence of the atom's neighbors\n")

        .def("GetBonds", AtomGetBonds,
             "Returns a read-only sequence of the atom's bonds\n")

        .def("Match", (bool (Atom::*)(const Atom *) const) & Atom::Match,
             "Returns whether or not this atom matches another Atom.\n\n"
             "  Each Atom (or query Atom) has a query function which is\n"
             "  used for this type of matching.\n\n"
             "  ARGUMENTS:\n"
             "    - other: the other Atom to which to compare\n")

        .def("IsInRingSize", AtomIsInRingSize,
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

        .def("GetSmarts", AtomGetSmarts,
             "returns the SMARTS (or SMILES) string for an Atom\n\n")

        // properties
        .def("SetProp", AtomSetProp<std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val")),
             "Sets an atomic property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n\n")

        .def("GetProp", GetProp<Atom, std::string>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
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

        .def("GetIntProp", GetProp<Atom, int>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (an int).\n\n"
             "  RETURNS: an int\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("GetUnsignedProp", GetProp<Atom, unsigned>,
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

        .def("GetBoolProp", GetProp<Atom, bool>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a bool).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")

        .def("HasProp", AtomHasProp,
             "Queries a Atom to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")

        .def("ClearProp", AtomClearProp,
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
              python::arg("includeComputed") = true),
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
                 1, python::with_custodian_and_ward_postcall<0, 1> >(),
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("GetPDBResidueInfo", AtomGetPDBResidueInfo,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1> >(),
             "Returns the atom's MonomerInfo object, if there is one.\n\n")
        .def("SetMonomerInfo", SetAtomMonomerInfo,
             "Sets the atom's MonomerInfo object.\n\n")
        .def("GetAtomMapNum", &Atom::getAtomMapNum,
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
        .value("SP3D", Atom::SP3D)
        .value("SP3D2", Atom::SP3D2)
        .value("OTHER", Atom::OTHER);
    python::enum_<Atom::ChiralType>("ChiralType")
        .value("CHI_UNSPECIFIED", Atom::CHI_UNSPECIFIED)
        .value("CHI_TETRAHEDRAL_CW", Atom::CHI_TETRAHEDRAL_CW)
        .value("CHI_TETRAHEDRAL_CCW", Atom::CHI_TETRAHEDRAL_CCW)
        .value("CHI_OTHER", Atom::CHI_OTHER)
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
    python::class_<QueryAtom, python::bases<Atom> >(
        "QueryAtom", atomClassDoc.c_str(), python::no_init)
        .def("ExpandQuery", expandQuery,
             (python::arg("self"), python::arg("other"),
              python::arg("how") = Queries::COMPOSITE_AND,
              python::arg("maintainOrder") = true),
             "combines the query from other with ours");

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
        "string.\n"
        ">>> m = Chem.MolFromSmiles(\"C\")\n"
        ">>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')\n"
        ">>> Chem.MolToSmiles(m)\n"
        "'C<xxx>'\n");
  }
};
}  // end of namespace
void wrap_atom() { RDKit::atom_wrapper::wrap(); }
