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
#include <span>
#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/pair.h>
#include <nanobind/make_iterator.h>
#include "props.hpp"
// #include "rdchem.h"
// #include "seqs.hpp"
#include "substructmethods.h"

// ours
#include <RDBoost/Wrap_nb.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
nb::bytes MolToBinary(const ROMol &self) {
  std::string res;
  {
    NOGIL gil;
    MolPickler::pickleMol(self, res);
  }
  nb::bytes retval = nb::bytes(res.c_str(), res.length());
  return retval;
}
nb::bytes MolToBinaryWithProps(const ROMol &self, unsigned int props) {
  std::string res;
  {
    NOGIL gil;
    MolPickler::pickleMol(self, res, props);
  }
  nb::bytes retval = nb::bytes(res.c_str(), res.length());
  return retval;
}
#if 0
void MolClearComputedPropsHelper(const ROMol &mol, bool includeRings) {
  mol.clearComputedProps(includeRings);
}


//
// allows molecules to be pickled.
//  since molecules have a constructor that takes a binary string
//  we only need to provide getinitargs()
//
struct mol_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const ROMol &self) {
    return python::make_tuple(MolToBinary(self));
  };
};

bool HasSubstructMatchStr(std::string pkl, const ROMol &query,
                          bool recursionPossible = true,
                          bool useChirality = false,
                          bool useQueryQueryMatches = false) {
  NOGIL gil;
  ROMol *mol;
  try {
    mol = new ROMol(pkl);
  } catch (...) {
    mol = nullptr;
  }
  if (!mol) {
    throw ValueErrorException("Null Molecule");
  }
  MatchVectType res;
  bool hasM = SubstructMatch(*mol, query, res, recursionPossible, useChirality,
                             useQueryQueryMatches);
  delete mol;
  return hasM;
}

unsigned int AddMolConformer(ROMol &mol, Conformer *conf,
                             bool assignId = false) {
  auto *nconf = new Conformer(*conf);
  return mol.addConformer(nconf, assignId);
}

Conformer *GetMolConformer(ROMol &mol, int id = -1) {
  return &(mol.getConformer(id));
}

// FIX: we should eventually figure out how to do iterators properly
QueryAtomIterSeq *MolGetAromaticAtoms(const ROMOL_SPTR &mol) {
  auto *qa = new QueryAtom();
  qa->setQuery(makeAtomAromaticQuery());
  QueryAtomIterSeq *res =
      new QueryAtomIterSeq(mol, mol->beginQueryAtoms(qa), mol->endQueryAtoms(),
                           AtomCountFunctor(mol));
  return res;
}
QueryAtomIterSeq *MolGetQueryAtoms(const ROMOL_SPTR &mol, QueryAtom *qa) {
  QueryAtomIterSeq *res =
      new QueryAtomIterSeq(mol, mol->beginQueryAtoms(qa), mol->endQueryAtoms(),
                           AtomCountFunctor(mol));
  return res;
}

ConformerIterSeq *GetMolConformers(const ROMOL_SPTR &mol) {
  ConformerIterSeq *res =
      new ConformerIterSeq(mol, mol->beginConformers(), mol->endConformers(),
                           ConformerCountFunctor(mol));
  return res;
}

int getMolNumAtoms(const ROMol &mol, int onlyHeavy, bool onlyExplicit) {
  if (onlyHeavy > -1) {
    BOOST_LOG(rdWarningLog)
        << "WARNING: the onlyHeavy argument to mol.GetNumAtoms() has been "
           "deprecated. Please use the onlyExplicit argument instead or "
           "mol.GetNumHeavyAtoms() if you want the heavy atom count."
        << std::endl;
    return mol.getNumAtoms(onlyHeavy);
  }
  return mol.getNumAtoms(onlyExplicit);
}
#endif

namespace {

void setSubstructMatchFinalCheck(SubstructMatchParameters &ps,
                                 nb::object func) {
  ps.extraFinalCheck = pyFinalMatchFunctor(func);
}

void setExtraAtomCheckFunc(
    SubstructMatchParameters &ps,
    std::function<bool(const Atom &queryAtom, const Atom &molAtom)> &func) {
  ps.extraAtomCheck = func;
}
void setExtraAtomCheckFunc2(SubstructMatchParameters &ps,
                            const AtomCoordsMatchFunctor &ftor) {
  ps.extraAtomCheck = std::bind(&AtomCoordsMatchFunctor::operator(), &ftor,
                                std::placeholders::_1, std::placeholders::_2);
}
void setExtraBondCheckFunc(
    SubstructMatchParameters &ps,
    std::function<bool(const Bond &queryBond, const Bond &molBond)> &func) {
  ps.extraBondCheck = func;
}

}  // namespace
namespace {
void MolDebug(const ROMol &mol, bool useStdout) {
  if (useStdout) {
    mol.debugMol(std::cout);
  } else {
    std::ostream *dest = &std::cerr;
    if (rdInfoLog != nullptr) {
      if (rdInfoLog->teestream) {
        dest = rdInfoLog->teestream;
      } else if (rdInfoLog->dp_dest) {
        dest = rdInfoLog->dp_dest;
      }
      mol.debugMol(*dest);
    }
  }
}

struct AtomSeqHolder {
  ROMol &d_mol;
  AtomSeqHolder(ROMol &mol) : d_mol(mol) {};
  size_t size() const { return d_mol.getNumAtoms(); }
  ROMol::AtomIterator begin() const { return d_mol.beginAtoms(); }
  ROMol::AtomIterator end() const { return d_mol.endAtoms(); }
};

}  // namespace

#if 0
class ReadWriteMol : public RWMol {
 public:
  ReadWriteMol() {};
  ReadWriteMol(const ROMol &m, bool quickCopy = false, int confId = -1)
      : RWMol(m, quickCopy, confId) {};

  void RemoveAtom(unsigned int idx) { removeAtom(idx); };
  void RemoveBond(unsigned int idx1, unsigned int idx2) {
    removeBond(idx1, idx2);
  };
  int AddBond(unsigned int begAtomIdx, unsigned int endAtomIdx,
              Bond::BondType order = Bond::UNSPECIFIED) {
    return addBond(begAtomIdx, endAtomIdx, order);
  };
  int AddAtom(Atom *atom) {
    PRECONDITION(atom, "bad atom");
    return addAtom(atom, true, false);
  };
  void ReplaceAtom(unsigned int idx, Atom *atom, bool updateLabel,
                   bool preserveProps) {
    PRECONDITION(atom, "bad atom");
    replaceAtom(idx, atom, updateLabel, preserveProps);
  };
  void ReplaceBond(unsigned int idx, Bond *bond, bool preserveProps,
                   bool keepSGroups) {
    PRECONDITION(bond, "bad bond");
    replaceBond(idx, bond, preserveProps, keepSGroups);
  };
  void SetStereoGroups(python::list &stereo_groups) {
    std::vector<StereoGroup> groups;
    pythonObjectToVect<StereoGroup>(stereo_groups, groups);
    for (const auto &group : groups) {
      for (const auto atom : group.getAtoms()) {
        if (!atom) {
          throw_value_error("NULL atom in StereoGroup");
        }
        if (&atom->getOwningMol() != this) {
          throw_value_error(
              "atom in StereoGroup does not belong to this molecule.");
        }
      }
    }
    setStereoGroups(std::move(groups));
  }
  ROMol *GetMol() const {
    auto *res = new ROMol(*this);
    return res;
  }

  ReadWriteMol *enter() {
    beginBatchEdit();
    return this;
  }
  bool exit(python::object exc_type, python::object exc_val,
            python::object traceback) {
    RDUNUSED_PARAM(exc_val);
    RDUNUSED_PARAM(traceback);
    if (exc_type != python::object()) {
      // exception thrown, abort the edits
      rollbackBatchEdit();
    } else {
      commitBatchEdit();
    }
    // we haven't handled any possible exceptions (and shouldn't do so),
    // so just return false;
    return false;
  }

 private:
  boost::shared_ptr<RWMol> dp_mol;
};
#endif

std::string molClassDoc =
    "The Molecule class.\n\n\
  In addition to the expected Atoms and Bonds, molecules contain:\n\
    - a collection of Atom and Bond bookmarks indexed with integers\n\
        that can be used to flag and retrieve particular Atoms or Bonds\n\
        using the {get|set}{Atom|Bond}Bookmark() methods.\n\n\
    - a set of string-valued properties. These can have arbitrary string\n\
        labels and can be set and retrieved using the {set|get}Prop() methods\n\
        Molecular properties can be tagged as being *computed*, in which case\n\
          they will be automatically cleared under certain circumstances (when the\n\
          molecule itself is modified, for example).\n\
        Molecules also have the concept of *private* properties, which are tagged\n\
          by beginning the property name with an underscore (_).\n";
std::string rwmolClassDoc =
    "The RW molecule class (read/write)\n\n\
  This class is a more-performant version of the EditableMolecule class in that\n\
  it is a 'live' molecule and shares the interface from the Mol class.\n\
  All changes are performed without the need to create a copy of the\n\
  molecule using GetMol() (this is still available, however).\n\
  \n\
  n.b. Eventually this class may become a direct replacement for EditableMol";

struct mol_wrapper {
  static void wrap(nb::module_ &m) {
    nb::exception<ConformerException>(m, "ConformerException",
                                      PyExc_ValueError);

    nb::enum_<RDKit::PicklerOps::PropertyPickleOptions>(m,
                                                        "PropertyPickleOptions")
        .value("NoProps", RDKit::PicklerOps::NoProps)
        .value("MolProps", RDKit::PicklerOps::MolProps)
        .value("AtomProps", RDKit::PicklerOps::AtomProps)
        .value("BondProps", RDKit::PicklerOps::BondProps)
        .value("QueryAtomData", RDKit::PicklerOps::QueryAtomData)
        .value("PrivateProps", RDKit::PicklerOps::PrivateProps)
        .value("ComputedProps", RDKit::PicklerOps::ComputedProps)
        .value("AllProps", RDKit::PicklerOps::AllProps)
        .value("CoordsAsDouble", RDKit::PicklerOps::CoordsAsDouble)
        .value("NoConformers", RDKit::PicklerOps::NoConformers)
        .export_values();

    //     RegisterVectorConverter<StereoGroup>("StereoGroup_vect");

    m.def("GetDefaultPickleProperties", MolPickler::getDefaultPickleProperties,
          "Get the current global mol pickler options.");
    m.def("SetDefaultPickleProperties", MolPickler::setDefaultPickleProperties,
          "arg1"_a, "Set the current global mol pickler options.");

    // REVIEW: There's probably a better place for the next few definitions
    nb::class_<RDKit::AtomCoordsMatchFunctor>(
        m, "AtomCoordsMatcher",
        "Allows using atom coordinates as part of substructure matching")
        .def(nb::init<>())
        .def(
            nb::init<int, int, double>(), "refConfId"_a = -1,
            "queryConfId"_a = -1, "tol"_a = 1e-4,
            "constructor taking reference and query conformer IDs and a distance tolerance")
        .def("__call__", &RDKit::AtomCoordsMatchFunctor::operator())
        .def_rw("refConfId", &RDKit::AtomCoordsMatchFunctor::d_refConfId,
                "reference conformer ID")
        .def_rw("queryConfId", &RDKit::AtomCoordsMatchFunctor::d_queryConfId,
                "query conformer ID")
        .def_rw("tol2", &RDKit::AtomCoordsMatchFunctor::d_tol2,
                "squared distance tolerance");
    nb::class_<RDKit::SubstructMatchParameters>(
        m, "SubstructMatchParameters",
        "Parameters controlling substructure matching")
        .def(nb::init<>(), "Constructor")
        .def_rw("useChirality", &RDKit::SubstructMatchParameters::useChirality,
                "Use chirality in determining whether or not atoms/bonds match")
        .def_rw(
            "useEnhancedStereo",
            &RDKit::SubstructMatchParameters::useEnhancedStereo,
            "take enhanced stereochemistry into account while doing the match. "
            "This only has an effect if useChirality is also True.")
        .def_rw("aromaticMatchesConjugated",
                &RDKit::SubstructMatchParameters::aromaticMatchesConjugated,
                "aromatic and conjugated bonds match each other")
        .def_rw("aromaticMatchesSingleOrDouble",
                &RDKit::SubstructMatchParameters::aromaticMatchesSingleOrDouble,
                "aromatic and single or double bonds match each other")
        .def_rw(
            "useGenericMatchers",
            &RDKit::SubstructMatchParameters::useGenericMatchers,
            "use generic groups (=homology groups) as a post-filtering step "
            "(if any are present in the molecule)")
        .def_rw("useQueryQueryMatches",
                &RDKit::SubstructMatchParameters::useQueryQueryMatches,
                "Consider query-query matches, not just simple matches")
        .def_rw("recursionPossible",
                &RDKit::SubstructMatchParameters::recursionPossible,
                "Allow recursive queries")
        .def_rw("uniquify", &RDKit::SubstructMatchParameters::uniquify,
                "uniquify (by atom index) match results")
        .def_rw("maxMatches", &RDKit::SubstructMatchParameters::maxMatches,
                "maximum number of matches to return")
        .def_rw("maxRecursiveMatches",
                &RDKit::SubstructMatchParameters::maxRecursiveMatches,
                "maximum number of recursive matches to find")
        .def_rw(
            "numThreads", &RDKit::SubstructMatchParameters::numThreads,
            "number of threads to use when multi-threading is possible."
            "0 selects the number of concurrent threads supported by the"
            "hardware. negative values are added to the number of concurrent"
            "threads supported by the hardware.")
        .def_rw("bondProperties",
                &RDKit::SubstructMatchParameters::bondProperties,
                "bond properties that must be equivalent in order to match.")
        .def_rw("atomProperties",
                &RDKit::SubstructMatchParameters::atomProperties,
                "atom properties that must be equivalent in order to match.")
        .def_rw(
            "specifiedStereoQueryMatchesUnspecified",
            &RDKit::SubstructMatchParameters::
                specifiedStereoQueryMatchesUnspecified,
            "If set, query atoms and bonds with specified stereochemistry will match atoms and bonds with unspecified stereochemistry.")
        .def("setExtraFinalCheck", setSubstructMatchFinalCheck,
             // FIX: This probably doesn't have the right custodian/ward
             "func"_a,
             R"DOC(allows you to provide a function that will be called
               with the molecule
           and a vector of atom IDs containing a potential match.
           The function should return true or false indicating whether or not
           that match should be accepted.)DOC")
        .def("setExtraAtomCheckFunc", setExtraAtomCheckFunc,
             // FIX: This probably doesn't have the right custodian/ward
             "func"_a,
             R"DOC(allows you to provide a function that will be called
           for each atom pair that matches during substructure searching,
           after all other comparisons have passed.
           The function should return true or false indicating whether or not
           that atom-match should be accepted.)DOC")
        .def(
            "setExtraAtomCheckFunc", setExtraAtomCheckFunc2,
            // FIX: This probably doesn't have the right custodian/ward
            "atomCoordsMatcher"_a,
            R"DOC(allows you to provide an AtomCoordsMatcher that will be called
           for each atom pair that matches during substructure searching,
           after all other comparisons have passed.)DOC")
        .def_rw(
            "extraAtomCheckOverridesDefaultCheck",
            &RDKit::SubstructMatchParameters::
                extraAtomCheckOverridesDefaultCheck,
            "if set, only the extraAtomCheck will be used to determine whether or not atoms match")
        .def("setExtraBondCheckFunc", setExtraBondCheckFunc,
             // FIX: This probably doesn't have the right custodian/ward
             "func"_a,
             R"DOC(allows you to provide a function that will be called
           for each bond pair that matches during substructure searching,
           after all other comparisons have passed.
           The function should return true or false indicating whether or not
           that bond-match should be accepted.)DOC")
        .def_rw(
            "extraBondCheckOverridesDefaultCheck",
            &RDKit::SubstructMatchParameters::
                extraBondCheckOverridesDefaultCheck,
            "if set, only the extraBondCheck will be used to determine whether or not bonds match")
        .def("__setattr__", &safeSetattr);

    nb::class_<AtomSeqHolder>(m, "_AtomSeqHolder",
                              "A sequence-like holder of a molecule's atoms")
        .def(nb::init<ROMol &>(), "mol"_a)
        .def("__len__", &AtomSeqHolder::size)
        .def(
            "__iter__",
            [](const AtomSeqHolder &a) {
              return nb::make_iterator(nb::type<AtomSeqHolder>(), "iterator",
                                       a.begin(), a.end());
            },
            nb::keep_alive<0, 1>());
    nb::class_<ROMol>(m, "Mol", nb::dynamic_attr())
        .def(nb::init<>(), "Constructor, takes no arguments")
        .def(
            "__init__",
            [](ROMol *t, nb::bytes b) {
              new (t) ROMol(std::string(static_cast<const char *>(b.data()),
                                        static_cast<size_t>(b.size())));
            },
            "Constructor from a binary string", "pklString"_a)
        .def(
            "__init__",
            [](ROMol *t, nb::bytes b, unsigned int propertyFlags) {
              new (t) ROMol(std::string(static_cast<const char *>(b.data()),
                                        static_cast<size_t>(b.size())),
                            propertyFlags);
            },
            "Constructor from a binary string", "pklString"_a,
            "propertyFlags"_a)
        .def(nb::init<const ROMol &, bool, int>(), "mol"_a,
             "quickCopy"_a = false, "confId"_a = -1)
        .def("__copy__", &generic__copy__<ROMol>)
        .def("__deepcopy__", &generic__deepcopy__<ROMol>, "memo"_a)
        .def(
            "GetAtoms", [](ROMol &mol) { return AtomSeqHolder(mol); },
            nb::keep_alive<0, 1>(),
            "Returns a sequence-like object of the molecule's atoms.")
        .def("GetNumAtoms",
             nb::overload_cast<>(&ROMol::getNumAtoms, nb::const_),
             "Returns the number of atoms in the molecule.")
        .def("GetNumHeavyAtoms", &ROMol::getNumHeavyAtoms,
             "Returns the number of heavy atoms (atomic number >1) in the "
             "molecule.")
        .def("GetAtomWithIdx",
             // nb::overload_cast doesn't seem to work with this one:
             (Atom * (ROMol::*)(unsigned int)) & ROMol::getAtomWithIdx,
             nb::rv_policy::reference_internal, "idx"_a,
             "Returns a particular Atom.\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which Atom to return\n\n"
             "  NOTE: atom indices start at 0")

        .def("GetNumBonds", &ROMol::getNumBonds, "onlyHeavy"_a = true,
             "Returns the number of Bonds in the molecule.\n\n"
             "  ARGUMENTS:\n"
             "    - onlyHeavy: (optional) include only bonds to heavy atoms "
             "(not Hs)\n"
             "                  defaults to True.")

        .def("GetBondWithIdx",
             (Bond * (ROMol::*)(unsigned int)) & ROMol::getBondWithIdx,
             nb::rv_policy::reference_internal, "idx"_a,
             "Returns a particular Bond.\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which Bond to return\n\n"
             "  NOTE: bond indices start at 0")

        .def("GetNumConformers", &ROMol::getNumConformers,
             "Return the number of conformations on the molecule")

        .def("AddConformer", &ROMol::addConformer, "conf"_a,
             "assignId"_a = false,
             "Add a conformer to the molecule and return the conformer ID")

        .def("GetConformer",
             (Conformer & (ROMol::*)(int)) & ROMol::getConformer, "id"_a = -1,
             "Get the conformer with a specified ID",
             nb::rv_policy::reference_internal)

        //    .def("GetConformers", GetMolConformers,
        //         python::return_value_policy<
        //             python::manage_new_object,
        //             python::with_custodian_and_ward_postcall<0, 1>>(),
        //         python::args("self"),
        //         "Returns a read-only sequence containing all of the
        //         molecule's " "Conformers.")

        .def("RemoveAllConformers", &ROMol::clearConformers,
             "Remove all the conformations on the molecule")

        .def("RemoveConformer", &ROMol::removeConformer, "id"_a,
             "Remove the conformer with the specified ID")
        .def("GetBondBetweenAtoms",
             (Bond * (ROMol::*)(unsigned int, unsigned int)) &
                 ROMol::getBondBetweenAtoms,
             nb::rv_policy::reference_internal, "idx1"_a, "idx2"_a,
             "Returns the bond between two atoms, if there is one.\n\n"
             "  ARGUMENTS:\n"
             "    - idx1,idx2: the Atom indices\n\n"
             "  Returns:\n"
             "    The Bond between the two atoms, if such a bond exists.\n"
             "    If there is no Bond between the atoms, None is returned "
             "instead.\n\n"
             "  NOTE: atom indices start at 0\n")

        .def("HasQuery", &ROMol::hasQuery,
             "Returns if any atom or bond in molecule has a query")

        // substructures
        .def("HasSubstructMatch",
             (bool (*)(const ROMol &m, const ROMol &query, bool, bool,
                       bool))HasSubstructMatch,
             "query"_a, "recursionPossible"_a = true, "useChirality"_a = false,
             "useQueryQueryMatches"_a = false,
             "Queries whether or not the molecule contains a particular "
             "substructure.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - recursionPossible: (optional)\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (nb::object(*)(const ROMol &m, const ROMol &query, bool,
                            bool))GetSubstructMatch,
             "query"_a, "useChirality"_a = false,
             "useQueryQueryMatches"_a = false,
             "Returns the indices of the molecule's atoms that match a "
             "substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "  RETURNS: a tuple of integers\n\n"
             "  NOTES:\n"
             "     - only a single match is returned\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        .def("GetSubstructMatches",
             (nb::object(*)(const ROMol &m, const ROMol &query, bool, bool,
                            bool, unsigned int))GetSubstructMatches,
             "query"_a, "uniquify"_a = true, "useChirality"_a = false,
             "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
             "Returns tuples of the indices of the molecule's atoms that "
             "match "
             "a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule.\n"
             "    - uniquify: (optional) determines whether or not the "
             "matches "
             "are uniquified.\n"
             "                Defaults to 1.\n\n"
             "    - useChirality: enables the use of stereochemistry in the "
             "matching\n\n"
             "    - useQueryQueryMatches: use query-query matching logic\n\n"
             "    - maxMatches: The maximum number of matches that will be "
             "returned.\n"
             "                  In high-symmetry cases with medium-sized "
             "molecules, it is\n"
             "                  very easy to end up with a combinatorial "
             "explosion in the\n"
             "                  number of possible matches. This argument "
             "prevents that from\n"
             "                  having unintended consequences\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        //    .def("HasSubstructMatch",
        //         (bool (*)(const ROMol &m, const MolBundle &query, bool, bool,
        //                   bool))HasSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //          python::arg("recursionPossible") = true,
        //          python::arg("useChirality") = false,
        //          python::arg("useQueryQueryMatches") = false))
        //    .def("GetSubstructMatch",
        //         (PyObject * (*)(const ROMol &m, const MolBundle &query, bool,
        //                         bool)) GetSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //          python::arg("useChirality") = false,
        //          python::arg("useQueryQueryMatches") = false))

        //    .def("GetSubstructMatches",
        //         (PyObject * (*)(const ROMol &m, const MolBundle &query, bool,
        //         bool,
        //                         bool, unsigned int)) GetSubstructMatches,
        //         (python::arg("self"), python::arg("query"),
        //          python::arg("uniquify") = true,
        //          python::arg("useChirality") = false,
        //          python::arg("useQueryQueryMatches") = false,
        //          python::arg("maxMatches") = 1000))

        //    //--------------------------------------------
        //    .def("HasSubstructMatch",
        //         (bool (*)(const ROMol &m, const ROMol &query,
        //                   const SubstructMatchParameters
        //                   &))helpHasSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //         python::arg("params")), "Queries whether or not the molecule
        //         contains a particular " "substructure.\n\n" "  ARGUMENTS:\n"
        //         " - query: a Molecule\n\n" "    - params: parameters
        //         controlling the substructure match\n\n" "  RETURNS: True or
        //         False\n")
        //    .def("GetSubstructMatch",
        //         (PyObject * (*)(const ROMol &m, const ROMol &query,
        //                         const SubstructMatchParameters &params))
        //             helpGetSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //         python::arg("params")), "Returns the indices of the
        //         molecule's atoms that match a " "substructure query.\n\n" "
        //         ARGUMENTS:\n" "    - query: a Molecule\n\n" "    - params:
        //         parameters controlling the substructure match\n\n" " RETURNS:
        //         a tuple of integers\n\n" "  NOTES:\n" "     - only a single
        //         match is returned\n" "     - the ordering of the indices
        //         corresponds to the atom " "ordering\n" "         in the
        //         query. For example, the first index is for the " "atom in\n"
        //         "         this molecule that matches the first atom in the "
        //         "query.\n")

        //    .def("GetSubstructMatches",
        //         (PyObject * (*)(const ROMol &m, const ROMol &query,
        //                         const SubstructMatchParameters &))
        //             helpGetSubstructMatches,
        //         (python::arg("self"), python::arg("query"),
        //         python::arg("params")), "Returns tuples of the indices of the
        //         molecule's atoms that " "match " "a substructure query.\n\n"
        //         " ARGUMENTS:\n" "    - query: a Molecule.\n" "    - params:
        //         parameters controlling the substructure match\n\n" " RETURNS:
        //         a tuple of tuples of integers\n\n" "  NOTE:\n" "     - the
        //         ordering of the indices corresponds to the atom "
        //         "ordering\n" " in the query. For example, the first index is
        //         for the " "atom in\n" " this molecule that matches the first
        //         atom in the " "query.\n")

        //    .def("HasSubstructMatch",
        //         (bool (*)(const ROMol &m, const MolBundle &query,
        //                   const SubstructMatchParameters
        //                   &))helpHasSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //          python::arg("params") = true))
        //    .def("GetSubstructMatch",
        //         (PyObject * (*)(const ROMol &m, const MolBundle &query,
        //                         const SubstructMatchParameters &))
        //             helpGetSubstructMatch,
        //         (python::arg("self"), python::arg("query"),
        //         python::arg("params")))
        //    .def("GetSubstructMatches",
        //         (PyObject * (*)(const ROMol &m, const MolBundle &query,
        //                         const SubstructMatchParameters &))
        //             helpGetSubstructMatches,
        //         (python::arg("self"), python::arg("query"),
        //         python::arg("params")))

        // properties
        .def("SetProp", MolSetProp<ROMol, std::string>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetDoubleProp", MolSetProp<ROMol, double>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a double valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a double.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to 0.\n\n")
        .def("SetIntProp", MolSetProp<ROMol, int>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets an integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an unsigned "
             "number).\n"
             "    - value: the property value as an integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetUnsignedProp", MolSetProp<ROMol, unsigned int>, "key"_a,
             "val"_a, "computed"_a = false,
             "Sets an unsigned integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as an unsigned integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetBoolProp", MolSetProp<ROMol, bool>, "key"_a, "val"_a,
             "computed"_a = false,
             "Sets a boolean valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a bool.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("HasProp", MolHasProp<ROMol>, "key"_a,
             "Queries a molecule to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")
        .def(
            "GetProp", GetPyProp<ROMol>, "key"_a, "autoConvert"_a = false,
            "Returns the value of the property.\n\n"
            "  ARGUMENTS:\n"
            "    - key: the name of the property to return (a string).\n\n"
            "    - autoConvert: if True attempt to convert the property into a python object\n\n"
            "  RETURNS: a string\n\n"
            "  NOTE:\n"
            "    - If the property has not been set, a KeyError exception will be raised.\n",
            boost::python::return_value_policy<return_pyobject_passthrough>())
        .def("GetDoubleProp", GetProp<ROMol, double>, "key"_a,
             "Returns the double value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetIntProp", GetProp<ROMol, int>, "key"_a,
             "Returns the integer value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetUnsignedProp", GetProp<ROMol, unsigned int>,
             python::args("self", "key"),
             "Returns the unsigned int value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an unsigned integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetBoolProp", GetProp<ROMol, bool>, "key"_a,
             "Returns the Bool value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n", )
        .def("ClearProp", MolClearProp<ROMol>, "key"_a,
             "Removes a property from the molecule.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to clear (a string).\n")

        .def("GetPropNames", &ROMol::getPropList, "includePrivate"_a = false,
             "includeComputed"_a = false,
             "Returns a tuple with all property names for this molecule.\n\n"
             "  ARGUMENTS:\n"
             "    - includePrivate: (optional) toggles inclusion of private "
             "properties in the result set.\n"
             "                      Defaults to 0.\n"
             "    - includeComputed: (optional) toggles inclusion of computed "
             "properties in the result set.\n"
             "                      Defaults to 0.\n\n"
             "  RETURNS: a tuple of strings\n")

        .def("GetPropsAsDict", GetPropsAsDict<ROMol>,
             "includePrivate"_a = false, "includeComputed"_a = false,
             "autoConvertStrings"_a = true, getPropsAsDictDocString.c_str())
        .def("GetStereoGroups", &ROMol::getStereoGroups,
             "Returns a list of StereoGroups defining the relative "
             "stereochemistry "
             "of the atoms.\n",
             nb::keep_alive<0, 1>(), nb::rv_policy::reference_internal)
#if 0
        .def("GetAromaticAtoms", MolGetAromaticAtoms,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             python::args("self"),
             "Returns a read-only sequence containing all of the molecule's "
             "aromatic Atoms.\n")
        .def(
            "GetAtomsMatchingQuery", MolGetQueryAtoms,
            python::return_value_policy<
                python::manage_new_object,
                python::with_custodian_and_ward_postcall<0, 1>>(),
            python::args("self", "qa"),
            "Returns a read-only sequence containing all of the atoms in a "
            "molecule that match the query atom. "
            "Atom query options are defined in the rdkit.Chem.rdqueries module.\n")

#endif

        .def("ClearComputedProps", &ROMol::clearComputedProps,
             "includeRings"_a = true,
             "Removes all computed properties from the molecule.\n\n")

        .def("UpdatePropertyCache", &ROMol::updatePropertyCache,
             "strict"_a = true,
             "Regenerates computed properties like implicit valence and ring "
             "information.\n\n")

        .def("NeedsUpdatePropertyCache", &ROMol::needsUpdatePropertyCache,
             "Returns true or false depending on whether implicit and "
             "explicit "
             "valence of the molecule have already been calculated.\n\n")

        .def(
            "ClearPropertyCache", &ROMol::clearPropertyCache,
            (python::arg("self")),
            "Clears implicit and explicit valence information from all atoms.\n\n")

        .def("Debug", MolDebug, "useStdout"_a = true,
             "Prints debugging information about the molecule.")

        .def("ToBinary", MolToBinary,
             "Returns a binary string representation of the molecule.")
        .def("ToBinary", MolToBinaryWithProps, "propertyFlags"_a,
             "Returns a binary string representation of the molecule pickling "
             "the "
             "specified properties.\n")

        .def("GetRingInfo", &ROMol::getRingInfo,
             nb::rv_policy::reference_internal,
             "Returns the number of molecule's RingInfo object.\n\n")
        .def("__getstate__",
             [](const ROMol &mol) {
               const auto pkl = MolToBinary(mol);
               return std::make_tuple(pkl);
             })
        .def("__setstate__",
             [](ROMol &mol, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&mol) ROMol(pkl);
             })
        .doc() = molClassDoc.c_str();

#if 0        
    // ---------------------------------------------------------------------------------------------
    python::def("_HasSubstructMatchStr", HasSubstructMatchStr,
                (python::arg("pkl"), python::arg("query"),
                 python::arg("recursionPossible") = true,
                 python::arg("useChirality") = false,
                 python::arg("useQueryQueryMatches") = false),
                "This function is included to speed substructure queries from "
                "databases, \n"
                "it's probably not of\n"
                "general interest.\n\n"
                "  ARGUMENTS:\n"
                "    - pkl: a Molecule pickle\n\n"
                "    - query: a Molecule\n\n"
                "    - recursionPossible: (optional)\n\n"
                "    - useChirality: (optional)\n\n"
                "    - useQueryQueryMatches: use query-query matching logic\n\n"
                "  RETURNS: True or False\n");

    python::class_<ReadWriteMol, python::bases<ROMol>>(
        "RWMol", rwmolClassDoc.c_str(),
        python::init<const ROMol &>(python::args("self", "m"),
                                    "Construct from a Mol"))
        .def(python::init<>(python::args("self")))
        .def(python::init<const std::string &>(
            python::args("self", "pklString")))
        .def(python::init<const std::string &, unsigned int>(
            (python::args("self", "pklString", "propertyFlags"))))
        .def(python::init<const ROMol &, bool, int>(
            (python::arg("self"), python::arg("mol"),
             python::arg("quickCopy") = false, python::arg("confId") = -1)))
        .def("__copy__", &generic__copy__<ReadWriteMol>, python::args("self"))
        .def("__deepcopy__", &generic__deepcopy__<ReadWriteMol>,
             python::args("self", "memo"))
        .def("__enter__", &ReadWriteMol::enter,
             python::return_internal_reference<>())
        .def("__exit__", &ReadWriteMol::exit)

        .def("RemoveAtom", &ReadWriteMol::RemoveAtom,
             python::args("self", "idx"),
             "Remove the specified atom from the molecule")
        .def("RemoveBond", &ReadWriteMol::RemoveBond,
             python::args("self", "idx1", "idx2"),
             "Remove the specified bond from the molecule")

        .def("AddBond", &ReadWriteMol::AddBond,
             ((python::arg("self"), python::arg("beginAtomIdx")),
              python::arg("endAtomIdx"),
              python::arg("order") = Bond::UNSPECIFIED),
             "add a bond, returns the new number of bonds")

        .def("AddAtom", &ReadWriteMol::AddAtom,
             ((python::arg("self"), python::arg("atom"))),
             "add an atom, returns the index of the newly added atom")
        .def("ReplaceAtom", &ReadWriteMol::ReplaceAtom,
             ((python::arg("self"), python::arg("index")),
              python::arg("newAtom"), python::arg("updateLabel") = false,
              python::arg("preserveProps") = false),
             "replaces the specified atom with the provided one\n"
             "If updateLabel is True, the new atom becomes the active atom\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new atom")
        .def("ReplaceBond", &ReadWriteMol::ReplaceBond,
             ((python::arg("self"), python::arg("index")),
              python::arg("newBond"), python::arg("preserveProps") = false,
              python::arg("keepSGroups") = true),
             "replaces the specified bond with the provided one.\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new bond. If keepSGroups is False, all"
             "Substance Groups referencing the bond will be dropped.")
        .def("GetMol", &ReadWriteMol::GetMol,
             "Returns a Mol (a normal molecule)",
             python::return_value_policy<python::manage_new_object>(),
             python::args("self"))

        .def("SetStereoGroups", &ReadWriteMol::SetStereoGroups,
             ((python::arg("self"), python::arg("stereo_groups"))),
             "Set the stereo groups")

        .def("InsertMol", &ReadWriteMol::insertMol,
             ((python::arg("self"), python::arg("mol"))),
             "Insert (add) the given molecule into this one")

        .def("BeginBatchEdit", &RWMol::beginBatchEdit, python::args("self"),
             "starts batch editing")
        .def("RollbackBatchEdit", &RWMol::rollbackBatchEdit,
             python::args("self"), "cancels batch editing")
        .def("CommitBatchEdit", &RWMol::commitBatchEdit, python::args("self"),
             "finishes batch editing and makes the actual changes")

        // enable pickle support
        .def_pickle(mol_pickle_suite());
#endif
  }
};
}  // namespace RDKit
void wrap_mol(nb::module_ &m) { RDKit::mol_wrapper::wrap(m); }
