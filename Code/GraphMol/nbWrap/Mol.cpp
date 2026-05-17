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
#include "seqholders.hpp"
#include "substructmethods.h"

// ours
#include <RDBoost/Wrap_nb.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/RDKitBase.h>
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
bool HasSubstructMatchStr(nb::bytes pkl, const ROMol &query,
                          bool recursionPossible = true,
                          bool useChirality = false,
                          bool useQueryQueryMatches = false) {
  NOGIL gil;
  std::unique_ptr<ROMol> mol;
  try {
    mol = std::make_unique<ROMol>(
        std::string(static_cast<const char *>(pkl.data()), pkl.size()));
  } catch (...) {
    mol = nullptr;
  }
  if (!mol) {
    throw ValueErrorException("Null Molecule");
  }
  MatchVectType res;
  bool hasM = SubstructMatch(*mol, query, res, recursionPossible, useChirality,
                             useQueryQueryMatches);
  return hasM;
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

}  // namespace

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
  void SetStereoGroups(nb::list &stereo_groups) {
    std::vector<StereoGroup> groups;
    pythonObjectToVect<StereoGroup>(stereo_groups, groups);
    for (const auto &group : groups) {
      for (const auto atom : group.getAtoms()) {
        if (!atom) {
          throw ValueErrorException("NULL atom in StereoGroup");
        }
        if (&atom->getOwningMol() != this) {
          throw ValueErrorException(
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
  bool exit(nb::object exc_type, nb::object *exc_val, nb::object traceback) {
    RDUNUSED_PARAM(exc_val);
    RDUNUSED_PARAM(traceback);
    if (!exc_type.is_none()) {
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
  std::shared_ptr<RWMol> dp_mol;
};

std::string molClassDoc = R"DOC(The Molecule class.

     In addition to the expected Atoms and Bonds, molecules contain:
          - a collection of Atom and Bond bookmarks indexed with integers
                    that can be used to flag and retrieve particular Atoms or Bonds
                    using the {get|set}{Atom|Bond}Bookmark() methods.

          - a set of string-valued properties. These can have arbitrary string
                    labels and can be set and retrieved using the {set|get}Prop() methods
                    Molecular properties can be tagged as being *computed*, in which case
                         they will be automatically cleared under certain circumstances (when the
                         molecule itself is modified, for example).
                    Molecules also have the concept of *private* properties, which are tagged
                         by beginning the property name with an underscore (_).
)DOC";
std::string rwmolClassDoc = R"DOC(The RW molecule class (read/write)

     This class is a more-performant version of the EditableMolecule class in that
     it is a 'live' molecule and shares the interface from the Mol class.
     All changes are performed without the need to create a copy of the
     molecule using GetMol() (this is still available, however).
  
     n.b. Eventually this class may become a direct replacement for EditableMol)DOC";

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

    nb::class_<AtomSeqHolder<>>(m, "_AtomSeqHolder1",
                                "A sequence-like holder of a molecule's atoms")
        .def(nb::init<ROMol &>(), "mol"_a)
        .def("__len__", &AtomSeqHolder<>::size)
        .def(
            "__iter__",
            [](AtomSeqHolder<> &a) {
              return nb::make_iterator(nb::type<AtomSeqHolder<>>(), "iterator",
                                       a.begin(), a.end());
            },
            nb::keep_alive<0, 1>())
        .def("__getitem__", &AtomSeqHolder<>::operator[],
             nb::rv_policy::reference_internal, "idx"_a);
    nb::class_<AtomSeqHolder<AtomNeighborsIterator>>(
        m, "_AtomSeqHolder2", "A sequence-like holder of an atom's neighbors")
        //    .def(nb::init<ROMol &>(), "mol"_a)
        .def("__len__", &AtomSeqHolder<AtomNeighborsIterator>::size)
        .def(
            "__iter__",
            [](AtomSeqHolder<AtomNeighborsIterator> &a) {
              return nb::make_iterator(
                  nb::type<AtomSeqHolder<AtomNeighborsIterator>>(), "iterator",
                  a.begin(), a.end());
            },
            nb::keep_alive<0, 1>())
        .def("__getitem__", &AtomSeqHolder<AtomNeighborsIterator>::operator[],
             nb::rv_policy::reference_internal, "idx"_a);
    nb::class_<BondSeqHolder<>>(m, "_BondSeqHolder1",
                                "A sequence-like holder of a molecule's bonds")
        .def(nb::init<ROMol &>(), "mol"_a)
        .def("__len__", &BondSeqHolder<>::size)
        .def(
            "__iter__",
            [](BondSeqHolder<> &a) {
              return nb::make_iterator(nb::type<BondSeqHolder<>>(), "iterator",
                                       a.begin(), a.end());
            },
            nb::keep_alive<0, 1>())
        .def("__getitem__", &BondSeqHolder<>::operator[],
             nb::rv_policy::reference_internal, "idx"_a);
    nb::class_<BondSeqHolder<AtomBondsIterator>>(
        m, "_BondSeqHolder2", "A sequence-like holder of an atom's bonds")
        .def("__len__", &BondSeqHolder<AtomBondsIterator>::size)
        .def(
            "__iter__",
            [](BondSeqHolder<AtomBondsIterator> &a) {
              return nb::make_iterator(
                  nb::type<BondSeqHolder<AtomBondsIterator>>(), "iterator",
                  a.begin(), a.end());
            },
            nb::keep_alive<0, 1>())
        .def("__getitem__", &BondSeqHolder<AtomBondsIterator>::operator[],
             nb::rv_policy::reference_internal, "idx"_a);
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
            "GetAtoms", [](ROMol &mol) { return AtomSeqHolder<>(mol); },
            nb::keep_alive<0, 1>(),
            "Returns a sequence-like object of the molecule's atoms.")
        .def(
            "GetBonds", [](ROMol &mol) { return BondSeqHolder<>(mol); },
            nb::keep_alive<0, 1>(),
            "Returns a sequence-like object of the molecule's bonds.")
        .def("GetNumAtoms",
             nb::overload_cast<>(&ROMol::getNumAtoms, nb::const_),
             "Returns the number of atoms in the molecule.")
        .def(
            "GetNumAtoms",
            nb::overload_cast<bool>(&ROMol::getNumAtoms, nb::const_),
            "Returns the number of atoms in the molecule. Optionally, only count explicit atoms.",
            "onlyExplicit"_a = true)
        .def("GetNumHeavyAtoms", &ROMol::getNumHeavyAtoms,
             "Returns the number of heavy atoms (atomic number >1) in the "
             "molecule.")
        .def("GetAtomWithIdx",
             // nb::overload_cast doesn't seem to work with this one:
             (Atom * (ROMol::*)(unsigned int)) & ROMol::getAtomWithIdx,
             nb::rv_policy::reference_internal, "idx"_a,
             R"DOC(Returns a particular Atom.

     ARGUMENTS:
          - idx: which Atom to return

     NOTE: atom indices start at 0)DOC")

        .def("GetNumBonds", &ROMol::getNumBonds, "onlyHeavy"_a = true,
             R"DOC(Returns the number of Bonds in the molecule.

     ARGUMENTS:
          - onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)
                                             defaults to True.)DOC")

        .def("GetBondWithIdx",
             (Bond * (ROMol::*)(unsigned int)) & ROMol::getBondWithIdx,
             nb::rv_policy::reference_internal, "idx"_a,
             R"DOC(Returns a particular Bond.

     ARGUMENTS:
          - idx: which Bond to return

     NOTE: bond indices start at 0)DOC")

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
             R"DOC(Returns the bond between two atoms, if there is one.

     ARGUMENTS:
          - idx1,idx2: the Atom indices

     Returns:
          The Bond between the two atoms, if such a bond exists.
          If there is no Bond between the atoms, None is returned instead.

     NOTE: atom indices start at 0
)DOC")

        .def("HasQuery", &ROMol::hasQuery,
             "Returns if any atom or bond in molecule has a query")

        // substructures
        .def(
            "HasSubstructMatch",
            (bool (*)(const ROMol &m, const ROMol &query, bool, bool,
                      bool))HasSubstructMatch,
            "query"_a, "recursionPossible"_a = true, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Queries whether or not the molecule contains a particular substructure.

     ARGUMENTS:
          - query: a Molecule

          - recursionPossible: (optional)

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: True or False
)DOC")
        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(const ROMol &m, const ROMol &query, bool,
                                  bool))GetSubstructMatch,
            "query"_a, "useChirality"_a = false,
            "useQueryQueryMatches"_a = false,
            R"DOC(Returns the indices of the molecule's atoms that match a substructure query.

     ARGUMENTS:
          - query: a Molecule

          - useChirality: enables the use of stereochemistry in the matching

          - useQueryQueryMatches: use query-query matching logic

     RETURNS: a list of integers

     NOTES:
           - only a single match is returned
           - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
)DOC")

        .def("GetSubstructMatches",
             (std::vector<std::vector<int>> (*)(
                 const ROMol &m, const ROMol &query, bool, bool, bool,
                 unsigned int))GetSubstructMatches,
             "query"_a, "uniquify"_a = true, "useChirality"_a = false,
             "useQueryQueryMatches"_a = false, "maxMatches"_a = 1000,
             R"DOC(Returns lists of the indices of the molecule's
                    atoms that match a substructure query.

             ARGUMENTS:
                  - query: a Molecule.
                  - uniquify: (optional) determines whether or not the
                  matches are uniquified.
                                                Defaults to 1.

                  - useChirality: enables the use of stereochemistry in the
                  matching

                  - useQueryQueryMatches: use query-query matching logic

                  - maxMatches: The maximum number of matches that will be
                  returned.
                                                     In high-symmetry cases
                                                     with medium-sized
                                                     molecules, it is very
                                                     easy to end up with a
                                                     combinatorial explosion
                                                     in the number of
                                                     possible matches. This
                                                     argument prevents that
                                                     from having unintended
                                                     consequences

             RETURNS: a list of lists of integers

             NOTE:
                   - the ordering of the indices corresponds to the atom
                   ordering
                             in the query. For example, the first index is
                             for the atom in this molecule that matches the
                             first atom in the query.
        )DOC")

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

        //--------------------------------------------
        .def(
            "HasSubstructMatch",
            (bool (*)(const ROMol &m, const ROMol &query,
                      const std::optional<SubstructMatchParameters>))
                helpHasSubstructMatch,
            "query"_a, "params"_a = nb::none(),
            R"DOC(Queries whether or not the molecule contains a particular substructure.
          
          ARGUMENTS:
          - query: a Molecule
          
          - params: parameters controlling the substructure match
          
          RETURNS: True or False
          )DOC")

        .def(
            "GetSubstructMatch",
            (std::vector<int> (*)(
                const ROMol &m, const ROMol &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatch,
            "query"_a, "params"_a = nb::none(),
            R"DOC(Returns the indices of the molecule's atoms that match a substructure query.

  ARGUMENTS:
    - query: a Molecule

    - params: parameters controlling the substructure match

  RETURNS: a list of integers

  NOTES:
     - only a single match is returned
     - the ordering of the indices corresponds to the atom ordering
         in the query. For example, the first index is for the atom in
         this molecule that matches the first atom in the query.
)DOC")

        .def(
            "GetSubstructMatches",
            (std::vector<std::vector<int>> (*)(
                const ROMol &m, const ROMol &query,
                const std::optional<SubstructMatchParameters>))
                helpGetSubstructMatches,
            "query"_a, "params"_a = nb::none(),
            R"DOC(Returns lists of the indices of the molecule's atoms that match a substructure query.

  ARGUMENTS:
    - query: a Molecule.

    - params: parameters controlling the substructure match

  RETURNS: a list of lists of integers

  NOTE:
     - the ordering of the indices corresponds to the atom ordering
         in the query. For example, the first index is for the atom in
         this molecule that matches the first atom in the query.
)DOC")
        //    .def("HasSubstructMatch",
        //         (bool (*)(const ROMol &m, const MolBundle &query,
        //                   const SubstructMatchParameters
        //                   &))helpHasSubstructMatch,
        //         "query"_a, "params"_a = SubstructMatchParameters())
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
             R"DOC(Sets a molecular property

     ARGUMENTS:
          - key: the name of the property to be set (a string).
          - value: the property value (a string).
          - computed: (optional) marks the property as being computed.
                                        Defaults to False.

)DOC")
        .def("SetDoubleProp", MolSetProp<ROMol, double>, "key"_a, "val"_a,
             "computed"_a = false,
             R"DOC(Sets a double valued molecular property

     ARGUMENTS:
          - key: the name of the property to be set (a string).
          - value: the property value as a double.
          - computed: (optional) marks the property as being computed.
                                        Defaults to 0.

)DOC")
        .def("SetIntProp", MolSetProp<ROMol, int>, "key"_a, "val"_a,
             "computed"_a = false,
             R"DOC(Sets an integer valued molecular property

     ARGUMENTS:
          - key: the name of the property to be set (an unsigned number).
          - value: the property value as an integer.
          - computed: (optional) marks the property as being computed.
                                        Defaults to False.

)DOC")
        .def("SetUnsignedProp", MolSetProp<ROMol, unsigned int>, "key"_a,
             "val"_a, "computed"_a = false,
             R"DOC(Sets an unsigned integer valued molecular property

     ARGUMENTS:
          - key: the name of the property to be set (a string).
          - value: the property value as an unsigned integer.
          - computed: (optional) marks the property as being computed.
                                        Defaults to False.

)DOC")
        .def("SetBoolProp", MolSetProp<ROMol, bool>, "key"_a, "val"_a,
             "computed"_a = false,
             R"DOC(Sets a boolean valued molecular property

     ARGUMENTS:
          - key: the name of the property to be set (a string).
          - value: the property value as a bool.
          - computed: (optional) marks the property as being computed.
                                        Defaults to False.

)DOC")
        .def(
            "HasProp", MolHasProp<ROMol>, "key"_a,
            R"DOC(Queries a molecule to see if a particular property has been assigned.

     ARGUMENTS:
          - key: the name of the property to check for (a string).
)DOC")
        .def("GetProp", GetPyProp<ROMol>, "key"_a, "autoConvert"_a = false,
             R"DOC(Returns the value of the property.

     ARGUMENTS:
          - key: the name of the property to return (a string).

          - autoConvert: if True attempt to convert the property into a python object

     RETURNS: a string

     NOTE:
          - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("GetDoubleProp", GetProp<ROMol, double>, "key"_a,
             R"DOC(Returns the double value of the property if possible.

     ARGUMENTS:
          - key: the name of the property to return (a string).

     RETURNS: a double

     NOTE:
          - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("GetIntProp", GetProp<ROMol, int>, "key"_a,
             R"DOC(Returns the integer value of the property if possible.

     ARGUMENTS:
          - key: the name of the property to return (a string).

     RETURNS: an integer

     NOTE:
          - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("GetUnsignedProp", GetProp<ROMol, unsigned int>, "key"_a,
             R"DOC(Returns the unsigned int value of the property if possible.

     ARGUMENTS:
          - key: the name of the property to return (a string).

     RETURNS: an unsigned integer

     NOTE:
          - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("GetBoolProp", GetProp<ROMol, bool>, "key"_a,
             R"DOC(Returns the Bool value of the property if possible.

     ARGUMENTS:
          - key: the name of the property to return (a string).

     RETURNS: a bool

     NOTE:
          - If the property has not been set, a KeyError exception will be raised.
)DOC")
        .def("ClearProp", MolClearProp<ROMol>, "key"_a,
             R"DOC(Removes a property from the molecule.

     ARGUMENTS:
          - key: the name of the property to clear (a string).
)DOC")

        .def("GetPropNames", &ROMol::getPropList, "includePrivate"_a = false,
             "includeComputed"_a = false,
             R"DOC(Returns a tuple with all property names for this molecule.

     ARGUMENTS:
          - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                                       Defaults to 0.
          - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                                       Defaults to 0.

     RETURNS: a tuple of strings
)DOC")

        .def("GetPropsAsDict", GetPropsAsDict<ROMol>,
             "includePrivate"_a = false, "includeComputed"_a = false,
             "autoConvertStrings"_a = true, getPropsAsDictDocString.c_str())
        .def(
            "GetStereoGroups", &ROMol::getStereoGroups,
            R"DOC(Returns a list of StereoGroups defining the relative stereochemistry of the atoms.
)DOC",
            nb::keep_alive<0, 1>(), nb::rv_policy::reference_internal)
#if 0
        .def("GetAromaticAtoms", MolGetAromaticAtoms, nb::keep_alive<0, 1>(),
             "Returns a read-only sequence containing all of the molecule's "
             "aromatic Atoms.\n")
             .def(
                  "GetAtomsMatchingQuery", MolGetQueryAtoms, "qa"_a,
                  nb::rv_policy::reference_internal,
                  R"DOC(Returns a read-only sequence containing all of the atoms in a molecule that match the query atom. 
                  Atom query options are defined in the rdkit.Chem.rdqueries module.
                  )DOC")
#endif

        .def("ClearComputedProps", &ROMol::clearComputedProps,
             "includeRings"_a = true,
             R"DOC(Removes all computed properties from the molecule.

)DOC")

        .def(
            "UpdatePropertyCache", &ROMol::updatePropertyCache,
            "strict"_a = true,
            R"DOC(Regenerates computed properties like implicit valence and ring information.

)DOC")

        .def(
            "NeedsUpdatePropertyCache", &ROMol::needsUpdatePropertyCache,
            R"DOC(Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.

)DOC")

        .def(
            "ClearPropertyCache", &ROMol::clearPropertyCache,
            R"DOC(Clears implicit and explicit valence information from all atoms.

)DOC")

        .def("Debug", MolDebug, "useStdout"_a = true,
             "Prints debugging information about the molecule.")

        .def("ToBinary", MolToBinary,
             "Returns a binary string representation of the molecule.")
        .def(
            "ToBinary", MolToBinaryWithProps, "propertyFlags"_a,
            R"DOC(Returns a binary string representation of the molecule pickling the specified properties.
)DOC")

        .def("GetRingInfo", &ROMol::getRingInfo,
             nb::rv_policy::reference_internal,
             R"DOC(Returns the number of molecule's RingInfo object.

)DOC")
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

    m.def(
        "_HasSubstructMatchStr", HasSubstructMatchStr, "pkl"_a, "query"_a,
        "recursionPossible"_a = true, "useChirality"_a = false,
        "useQueryQueryMatches"_a = false,
        R"DOC(This function is included to speed substructure queries from databases,
it's probably not of general interest.

     ARGUMENTS:
          - pkl: a Molecule pickle
          - query: a Molecule
          - recursionPossible: (optional)
          - useChirality: (optional)
          - useQueryQueryMatches: use query-query matching logic

     RETURNS: True or False
)DOC");
    // ---------------------------------------------------------------------------------------------

    nb::class_<ReadWriteMol, ROMol>(m, "RWMol", nb::dynamic_attr())
        .def(nb::init<>(), "default constructor, takes no arguments")
        .def(
            "__init__",
            [](ReadWriteMol *t, nb::bytes b) {
              new (t)
                  ReadWriteMol(std::string(static_cast<const char *>(b.data()),
                                           static_cast<size_t>(b.size())));
            },
            "pklString"_a, "Constructor from a binary string")

        .def(
            "__init__",
            [](ReadWriteMol *t, nb::bytes b, unsigned int propertyFlags) {
              new (t)
                  ReadWriteMol(std::string(static_cast<const char *>(b.data()),
                                           static_cast<size_t>(b.size())),
                               propertyFlags);
            },
            "pklString"_a, "propertyFlags"_a,
            "Constructor from a binary string")

        .def(
            "__init__",
            [](ReadWriteMol *t, const ROMol &m, bool quickCopy, int confId) {
              new (t) ReadWriteMol(m, quickCopy, confId);
            },
            "mol"_a, "quickCopy"_a = false, "confId"_a = -1,
            "Constructor from a ROMol")

        .def("__copy__", &generic__copy__<ReadWriteMol>)
        .def("__deepcopy__", &generic__deepcopy__<ReadWriteMol>, "memo"_a)
        .def("__enter__", &ReadWriteMol::enter,
             nb::rv_policy::reference_internal)
        .def("__exit__", &ReadWriteMol::exit, "excType"_a = nb::none(),
             "excValue"_a = nb::none(), "traceback"_a = nb::none())
        .def("RemoveAtom", &ReadWriteMol::RemoveAtom, "idx"_a,
             "Remove the specified atom from the molecule")
        .def("RemoveBond", &ReadWriteMol::RemoveBond, "idx1"_a, "idx2"_a,
             "Remove the specified bond from the molecule")

        .def("AddBond", &ReadWriteMol::AddBond, "beginAtomIdx"_a,
             "endAtomIdx"_a, "order"_a = Bond::UNSPECIFIED,
             "add a bond, returns the new number of bonds")
        .def("AddAtom", &ReadWriteMol::AddAtom, "atom"_a,
             "add an atom, returns the index of the newly added atom")
        .def("ReplaceAtom", &ReadWriteMol::ReplaceAtom, "index"_a, "newAtom"_a,
             "updateLabel"_a = false, "preserveProps"_a = false,
             "replaces the specified atom with the provided one\n"
             "If updateLabel is True, the new atom becomes the active atom\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new atom")
        .def("ReplaceBond", &ReadWriteMol::ReplaceBond, "index"_a, "newBond"_a,
             "preserveProps"_a = false, "keepSGroups"_a = true,
             "replaces the specified bond with the provided one.\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new bond. If keepSGroups is False, all"
             "Substance Groups referencing the bond will be dropped.")
        .def("GetMol", &ReadWriteMol::GetMol,
             "Returns a Mol (a normal molecule)", nb::rv_policy::take_ownership)

        .def("SetStereoGroups", &ReadWriteMol::SetStereoGroups,
             "stereo_groups"_a, "Set the stereo groups")

        .def("InsertMol", &ReadWriteMol::insertMol, "mol"_a,
             "Insert (add) the given molecule into this one")

        .def("BeginBatchEdit", &RWMol::beginBatchEdit, "starts batch editing")
        .def("RollbackBatchEdit", &RWMol::rollbackBatchEdit,
             "cancels batch editing")
        .def("CommitBatchEdit", &RWMol::commitBatchEdit,
             "finishes batch editing and makes the actual changes")

        .def("__getstate__",
             [](const ReadWriteMol &mol) {
               const auto pkl = MolToBinary(mol);
               return std::make_tuple(pkl);
             })
        .def("__setstate__",
             [](ReadWriteMol &mol, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&mol) ReadWriteMol(pkl);
             });
  }
};
}  // namespace RDKit
void wrap_mol(nb::module_ &m) { RDKit::mol_wrapper::wrap(m); }
