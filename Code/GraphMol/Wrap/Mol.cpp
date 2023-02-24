
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

#include "rdchem.h"
#include "seqs.hpp"
#include "props.hpp"
#include "substructmethods.h"

// ours
#include <RDBoost/pyint_api.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/python/iterator.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace python = boost::python;

namespace RDKit {

void MolClearComputedPropsHelper(const ROMol &mol, bool includeRings) {
  mol.clearComputedProps(includeRings);
}

python::object MolToBinary(const ROMol &self) {
  std::string res;
  {
    NOGIL gil;
    MolPickler::pickleMol(self, res);
  }
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object MolToBinaryWithProps(const ROMol &self, unsigned int props) {
  std::string res;
  {
    NOGIL gil;
    MolPickler::pickleMol(self, res, props);
  }
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
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

// FIX: we should eventually figure out how to do iterators properly
AtomIterSeq *MolGetAtoms(const ROMOL_SPTR &mol) {
  AtomIterSeq *res = new AtomIterSeq(mol, mol->beginAtoms(), mol->endAtoms(),
                                     AtomCountFunctor(mol));
  return res;
}
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

// AtomIterSeq *MolGetHeteros(ROMol *mol){
//  AtomIterSeq *res = new AtomIterSeq(mol->beginHeteros(),
//                                     mol->endHeteros());
//  return res;
//}
BondIterSeq *MolGetBonds(const ROMOL_SPTR &mol) {
  BondIterSeq *res = new BondIterSeq(mol, mol->beginBonds(), mol->endBonds(),
                                     BondCountFunctor(mol));
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

namespace {
class pyobjFunctor {
 public:
  pyobjFunctor(python::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() = default;
  bool operator()(const ROMol &m, const std::vector<unsigned int> &match) {
    return python::extract<bool>(dp_obj(boost::ref(m), boost::ref(match)));
  }

 private:
  python::object dp_obj;
};
void setSubstructMatchFinalCheck(SubstructMatchParameters &ps,
                                 python::object func) {
  ps.extraFinalCheck = pyobjFunctor(func);
}
}  // namespace

class ReadWriteMol : public RWMol {
 public:
  ReadWriteMol(){};
  ReadWriteMol(const ROMol &m, bool quickCopy = false, int confId = -1)
      : RWMol(m, quickCopy, confId){};

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
  static void wrap() {
    python::register_exception_translator<ConformerException>(
        &rdExceptionTranslator);

    python::enum_<RDKit::PicklerOps::PropertyPickleOptions>(
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
    ;

    RegisterVectorConverter<StereoGroup>("StereoGroup_vect");

    python::def("GetDefaultPickleProperties",
                MolPickler::getDefaultPickleProperties,
                "Get the current global mol pickler options.");
    python::def("SetDefaultPickleProperties",
                MolPickler::setDefaultPickleProperties,
                "Set the current global mol pickler options.");

    // REVIEW: There's probably a better place for this definition
    python::class_<RDKit::SubstructMatchParameters, boost::noncopyable>(
        "SubstructMatchParameters",
        "Parameters controlling substructure matching")
        .def_readwrite(
            "useChirality", &RDKit::SubstructMatchParameters::useChirality,
            "Use chirality in determining whether or not atoms/bonds match")
        .def_readwrite(
            "useEnhancedStereo",
            &RDKit::SubstructMatchParameters::useEnhancedStereo,
            "take enhanced stereochemistry into account while doing the match. "
            "This only has an effect if useChirality is also True.")
        .def_readwrite(
            "aromaticMatchesConjugated",
            &RDKit::SubstructMatchParameters::aromaticMatchesConjugated,
            "aromatic and conjugated bonds match each other")
        .def_readwrite(
            "useGenericMatchers",
            &RDKit::SubstructMatchParameters::useGenericMatchers,
            "use generic groups (=homology groups) as a post-filtering step "
            "(if any are present in the molecule)")
        .def_readwrite("useQueryQueryMatches",
                       &RDKit::SubstructMatchParameters::useQueryQueryMatches,
                       "Consider query-query matches, not just simple matches")
        .def_readwrite("recursionPossible",
                       &RDKit::SubstructMatchParameters::recursionPossible,
                       "Allow recursive queries")
        .def_readwrite("uniquify", &RDKit::SubstructMatchParameters::uniquify,
                       "uniquify (by atom index) match results")
        .def_readwrite("maxMatches",
                       &RDKit::SubstructMatchParameters::maxMatches,
                       "maximum number of matches to return")
        .def_readwrite(
            "numThreads", &RDKit::SubstructMatchParameters::numThreads,
            "number of threads to use when multi-threading is possible."
            "0 selects the number of concurrent threads supported by the"
            "hardware. negative values are added to the number of concurrent"
            "threads supported by the hardware.")
        .def("setExtraFinalCheck", setSubstructMatchFinalCheck,
             python::with_custodian_and_ward<1, 2>(),
             R"DOC(allows you to provide a function that will be called
               with the molecule
           and a vector of atom IDs containing a potential match.
           The function should return true or false indicating whether or not
           that match should be accepted.)DOC");

    python::class_<ROMol, ROMOL_SPTR, boost::noncopyable>(
        "Mol", molClassDoc.c_str(),
        python::init<>("Constructor, takes no arguments"))
        .def(python::init<const std::string &>(python::args("pklString")))
        .def(python::init<const std::string &, unsigned int>(
            (python::args("pklString", "propertyFlags"))))
        .def(python::init<const ROMol &, bool, int>(
            (python::arg("mol"), python::arg("quickCopy") = false,
             python::arg("confId") = -1)))
        .def("__copy__", &generic__copy__<ROMol>)
        .def("__deepcopy__", &generic__deepcopy__<ROMol>)
        .def(
            "GetNumAtoms", getMolNumAtoms,
            (python::arg("onlyHeavy") = -1, python::arg("onlyExplicit") = true),
            "Returns the number of atoms in the molecule.\n\n"
            "  ARGUMENTS:\n"
            "    - onlyExplicit: (optional) include only explicit atoms "
            "(atoms in the molecular graph)\n"
            "                    defaults to 1.\n"
            "  NOTE: the onlyHeavy argument is deprecated\n"

            )
        .def("GetNumHeavyAtoms", &ROMol::getNumHeavyAtoms,
             "Returns the number of heavy atoms (atomic number >1) in the "
             "molecule.\n\n")
        .def("GetAtomWithIdx",
             (Atom * (ROMol::*)(unsigned int)) & ROMol::getAtomWithIdx,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a particular Atom.\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which Atom to return\n\n"
             "  NOTE: atom indices start at 0\n")

        .def("GetNumBonds", &ROMol::getNumBonds,
             (python::arg("onlyHeavy") = true),
             "Returns the number of Bonds in the molecule.\n\n"
             "  ARGUMENTS:\n"
             "    - onlyHeavy: (optional) include only bonds to heavy atoms "
             "(not Hs)\n"
             "                  defaults to 1.\n")

        .def("GetBondWithIdx",
             (Bond * (ROMol::*)(unsigned int)) & ROMol::getBondWithIdx,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a particular Bond.\n\n"
             "  ARGUMENTS:\n"
             "    - idx: which Bond to return\n\n"
             "  NOTE: bond indices start at 0\n")

        .def("GetNumConformers", &ROMol::getNumConformers,
             "Return the number of conformations on the molecule")

        .def("AddConformer", AddMolConformer,
             (python::arg("self"), python::arg("conf"),
              python::arg("assignId") = false),
             "Add a conformer to the molecule and return the conformer ID")

#if 0
        .def("AddConformersFromTrajectory", &ROMol::addConformersFromTrajectory,
             (python::arg("self"), python::arg("traj"),
              python::arg("nConf") = -1),
             "Add conformers from a Trajectory to the molecule and return\n"
             "the number of conformations that were added")
#endif

        .def("GetConformer", GetMolConformer,
             (python::arg("self"), python::arg("id") = -1),
             "Get the conformer with a specified ID",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())

        .def("GetConformers", GetMolConformers,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a read-only sequence containing all of the molecule's "
             "Conformers.")

        .def("RemoveAllConformers", &ROMol::clearConformers,
             "Remove all the conformations on the molecule")

        .def("RemoveConformer", &ROMol::removeConformer,
             "Remove the conformer with the specified ID")
        .def("GetBondBetweenAtoms",
             (Bond * (ROMol::*)(unsigned int, unsigned int)) &
                 ROMol::getBondBetweenAtoms,
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns the bond between two atoms, if there is one.\n\n"
             "  ARGUMENTS:\n"
             "    - idx1,idx2: the Atom indices\n\n"
             "  Returns:\n"
             "    The Bond between the two atoms, if such a bond exists.\n"
             "    If there is no Bond between the atoms, None is returned "
             "instead.\n\n"
             "  NOTE: bond indices start at 0\n")

        // substructures
        .def("HasSubstructMatch",
             (bool (*)(const ROMol &m, const ROMol &query, bool, bool,
                       bool))HasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
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
             (PyObject * (*)(const ROMol &m, const ROMol &query, bool, bool))
                 GetSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false),
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
             (PyObject * (*)(const ROMol &m, const ROMol &query, bool, bool,
                             bool, unsigned int)) GetSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000),
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

        .def("HasSubstructMatch",
             (bool (*)(const ROMol &m, const MolBundle &query, bool, bool,
                       bool))HasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("recursionPossible") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))
        .def("GetSubstructMatch",
             (PyObject * (*)(const ROMol &m, const MolBundle &query, bool,
                             bool)) GetSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false))

        .def("GetSubstructMatches",
             (PyObject * (*)(const ROMol &m, const MolBundle &query, bool, bool,
                             bool, unsigned int)) GetSubstructMatches,
             (python::arg("self"), python::arg("query"),
              python::arg("uniquify") = true,
              python::arg("useChirality") = false,
              python::arg("useQueryQueryMatches") = false,
              python::arg("maxMatches") = 1000))

        //--------------------------------------------
        .def("HasSubstructMatch",
             (bool (*)(const ROMol &m, const ROMol &query,
                       const SubstructMatchParameters &))helpHasSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Queries whether or not the molecule contains a particular "
             "substructure.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: True or False\n")
        .def("GetSubstructMatch",
             (PyObject * (*)(const ROMol &m, const ROMol &query,
                             const SubstructMatchParameters &params))
                 helpGetSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns the indices of the molecule's atoms that match a "
             "substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule\n\n"
             "    - params: parameters controlling the substructure match\n\n"
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
             (PyObject * (*)(const ROMol &m, const ROMol &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatches,
             (python::arg("self"), python::arg("query"), python::arg("params")),
             "Returns tuples of the indices of the molecule's atoms that "
             "match "
             "a substructure query.\n\n"
             "  ARGUMENTS:\n"
             "    - query: a Molecule.\n"
             "    - params: parameters controlling the substructure match\n\n"
             "  RETURNS: a tuple of tuples of integers\n\n"
             "  NOTE:\n"
             "     - the ordering of the indices corresponds to the atom "
             "ordering\n"
             "         in the query. For example, the first index is for the "
             "atom in\n"
             "         this molecule that matches the first atom in the "
             "query.\n")

        .def("HasSubstructMatch",
             (bool (*)(const ROMol &m, const MolBundle &query,
                       const SubstructMatchParameters &))helpHasSubstructMatch,
             (python::arg("self"), python::arg("query"),
              python::arg("params") = true))
        .def("GetSubstructMatch",
             (PyObject * (*)(const ROMol &m, const MolBundle &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatch,
             (python::arg("self"), python::arg("query"), python::arg("params")))
        .def("GetSubstructMatches",
             (PyObject * (*)(const ROMol &m, const MolBundle &query,
                             const SubstructMatchParameters &))
                 helpGetSubstructMatches,
             (python::arg("self"), python::arg("query"), python::arg("params")))

        // properties
        .def("SetProp", MolSetProp<ROMol, std::string>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value (a string).\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetDoubleProp", MolSetProp<ROMol, double>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a double valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a double.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to 0.\n\n")
        .def("SetIntProp", MolSetProp<ROMol, int>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets an integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (an unsigned "
             "number).\n"
             "    - value: the property value as an integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetUnsignedProp", MolSetProp<ROMol, unsigned int>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets an unsigned integer valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as an unsigned integer.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("SetBoolProp", MolSetProp<ROMol, bool>,
             (python::arg("self"), python::arg("key"), python::arg("val"),
              python::arg("computed") = false),
             "Sets a boolean valued molecular property\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to be set (a string).\n"
             "    - value: the property value as a bool.\n"
             "    - computed: (optional) marks the property as being "
             "computed.\n"
             "                Defaults to False.\n\n")
        .def("HasProp", MolHasProp<ROMol>,
             "Queries a molecule to see if a particular property has been "
             "assigned.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to check for (a string).\n")
        .def("GetProp", GetProp<ROMol, std::string>,
             "Returns the value of the property.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a string\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetDoubleProp", GetProp<ROMol, double>,
             "Returns the double value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a double\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetIntProp", GetProp<ROMol, int>,
             "Returns the integer value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetUnsignedProp", GetProp<ROMol, unsigned int>,
             "Returns the unsigned int value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: an unsigned integer\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("GetBoolProp", GetProp<ROMol, bool>,
             "Returns the Bool value of the property if possible.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to return (a string).\n\n"
             "  RETURNS: a bool\n\n"
             "  NOTE:\n"
             "    - If the property has not been set, a KeyError exception "
             "will be raised.\n")
        .def("ClearProp", MolClearProp<ROMol>,
             "Removes a property from the molecule.\n\n"
             "  ARGUMENTS:\n"
             "    - key: the name of the property to clear (a string).\n")

        .def("ClearComputedProps", MolClearComputedPropsHelper,
             (python::arg("self"), python::arg("includeRings") = true),
             "Removes all computed properties from the molecule.\n\n")

        .def("UpdatePropertyCache", &ROMol::updatePropertyCache,
             (python::arg("self"), python::arg("strict") = true),
             "Regenerates computed properties like implicit valence and ring "
             "information.\n\n")

        .def("NeedsUpdatePropertyCache", &ROMol::needsUpdatePropertyCache,
             (python::arg("self")),
             "Returns true or false depending on whether implicit and "
             "explicit "
             "valence of the molecule have already been calculated.\n\n")

        .def("GetStereoGroups", &ROMol::getStereoGroups,
             "Returns a list of StereoGroups defining the relative "
             "stereochemistry "
             "of the atoms.\n",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())

        .def("GetPropNames", &ROMol::getPropList,
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
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
             (python::arg("self"), python::arg("includePrivate") = false,
              python::arg("includeComputed") = false),
             "Returns a dictionary populated with the molecules properties.\n"
             " n.b. Some properties are not able to be converted to python "
             "types.\n\n"
             "  ARGUMENTS:\n"
             "    - includePrivate: (optional) toggles inclusion of private "
             "properties in the result set.\n"
             "                      Defaults to False.\n"
             "    - includeComputed: (optional) toggles inclusion of computed "
             "properties in the result set.\n"
             "                      Defaults to False.\n\n"
             "  RETURNS: a dictionary\n")

        .def("GetAtoms", MolGetAtoms,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a read-only sequence containing all of the molecule's "
             "Atoms.\n")
        .def("GetAromaticAtoms", MolGetAromaticAtoms,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a read-only sequence containing all of the molecule's "
             "aromatic Atoms.\n")
        .def("GetAtomsMatchingQuery", MolGetQueryAtoms,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a read-only sequence containing all of the atoms in a "
             "molecule that match the query atom.\n")

        .def("GetBonds", MolGetBonds,
             python::return_value_policy<
                 python::manage_new_object,
                 python::with_custodian_and_ward_postcall<0, 1>>(),
             "Returns a read-only sequence containing all of the molecule's "
             "Bonds.\n")

        // enable pickle support
        .def_pickle(mol_pickle_suite())

        .def("Debug", MolDebug,
             (python::arg("mol"), python::arg("useStdout") = true),
             "Prints debugging information about the molecule.\n")

        .def("ToBinary", MolToBinary,
             "Returns a binary string representation of the molecule.\n")
        .def("ToBinary", MolToBinaryWithProps,
             (python::arg("mol"), python::arg("propertyFlags")),
             "Returns a binary string representation of the molecule pickling "
             "the "
             "specified properties.\n")

        .def("GetRingInfo", &ROMol::getRingInfo,
             python::return_value_policy<python::reference_existing_object>(),
             "Returns the number of molecule's RingInfo object.\n\n");
    python::register_ptr_to_python<std::shared_ptr<ROMol>>();

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
        python::init<const ROMol &>("Construct from a Mol"))
        .def(python::init<>())
        .def(python::init<const std::string &>(python::args("pklString")))
        .def(python::init<const std::string &, unsigned int>(
            (python::args("pklString", "propertyFlags"))))
        .def(python::init<const ROMol &, bool, int>(
            (python::arg("mol"), python::arg("quickCopy") = false,
             python::arg("confId") = -1)))
        .def("__copy__", &generic__copy__<ReadWriteMol>)
        .def("__deepcopy__", &generic__deepcopy__<ReadWriteMol>)
        .def("__enter__", &ReadWriteMol::enter,
             python::return_internal_reference<>())
        .def("__exit__", &ReadWriteMol::exit)

        .def("RemoveAtom", &ReadWriteMol::RemoveAtom,
             "Remove the specified atom from the molecule")
        .def("RemoveBond", &ReadWriteMol::RemoveBond,
             "Remove the specified bond from the molecule")

        .def("AddBond", &ReadWriteMol::AddBond,
             (python::arg("beginAtomIdx"), python::arg("endAtomIdx"),
              python::arg("order") = Bond::UNSPECIFIED),
             "add a bond, returns the new number of bonds")

        .def("AddAtom", &ReadWriteMol::AddAtom, (python::arg("atom")),
             "add an atom, returns the index of the newly added atom")
        .def("ReplaceAtom", &ReadWriteMol::ReplaceAtom,
             (python::arg("index"), python::arg("newAtom"),
              python::arg("updateLabel") = false,
              python::arg("preserveProps") = false),
             "replaces the specified atom with the provided one\n"
             "If updateLabel is True, the new atom becomes the active atom\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new atom")
        .def("ReplaceBond", &ReadWriteMol::ReplaceBond,
             (python::arg("index"), python::arg("newBond"),
              python::arg("preserveProps") = false,
              python::arg("keepSGroups") = true),
             "replaces the specified bond with the provided one.\n"
             "If preserveProps is True preserve keep the existing props unless "
             "explicit set on the new bond. If keepSGroups is False, all"
             "Substance Groups referencing the bond will be dropped.")
        .def("GetMol", &ReadWriteMol::GetMol,
             "Returns a Mol (a normal molecule)",
             python::return_value_policy<python::manage_new_object>())

        .def("SetStereoGroups", &ReadWriteMol::SetStereoGroups,
             (python::arg("stereo_groups")), "Set the stereo groups")

        .def("InsertMol", &ReadWriteMol::insertMol, (python::arg("mol")),
             "Insert (add) the given molecule into this one")

        .def("BeginBatchEdit", &RWMol::beginBatchEdit, "starts batch editing")
        .def("RollbackBatchEdit", &RWMol::rollbackBatchEdit,
             "cancels batch editing")
        .def("CommitBatchEdit", &RWMol::commitBatchEdit,
             "finishes batch editing and makes the actual changes")

        // enable pickle support
        .def_pickle(mol_pickle_suite());
  };
};
}  // namespace RDKit
void wrap_mol() { RDKit::mol_wrapper::wrap(); }
