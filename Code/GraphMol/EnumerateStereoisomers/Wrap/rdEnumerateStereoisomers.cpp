//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>

namespace python = boost::python;

namespace {
class LocalStereoEnumerator {
 public:
  LocalStereoEnumerator() = delete;
  LocalStereoEnumerator(python::object &py_mol, python::object &py_options,
                        const bool verbose) {
    RDKit::EnumerateStereoisomers::StereoEnumerationOptions opts;
    if (!py_options.is_none()) {
      opts = python::extract<
          RDKit::EnumerateStereoisomers::StereoEnumerationOptions>(py_options);
    }
    RDKit::ROMol mol = python::extract<RDKit::ROMol>(py_mol);
    dp_enumerator.reset(
        new RDKit::EnumerateStereoisomers::StereoisomerEnumerator(mol, opts,
                                                                  verbose));
  }
  LocalStereoEnumerator(const LocalStereoEnumerator &other) = delete;
  LocalStereoEnumerator(LocalStereoEnumerator &&other) = delete;
  LocalStereoEnumerator &operator=(const LocalStereoEnumerator &other) = delete;
  LocalStereoEnumerator &operator=(LocalStereoEnumerator &&other) = delete;
  ~LocalStereoEnumerator() = default;

  boost::shared_ptr<RDKit::ROMol> next() {
    auto iso = dp_enumerator->next();
    if (!iso) {
      return boost::shared_ptr<RDKit::ROMol>();
    }
    return boost::shared_ptr<RDKit::ROMol>(iso.release());
  }

  boost::shared_ptr<RDKit::ROMol> nextForIter() {
    auto iso = next();
    if (!iso) {
      PyErr_SetString(PyExc_StopIteration,
                      "Stereoisomer enumeration exhausted");
      boost::python::throw_error_already_set();
    }
    return iso;
  }

  uint64_t GetStereoisomerCount() {
    return dp_enumerator->getStereoisomerCount();
  }

 private:
  std::unique_ptr<RDKit::EnumerateStereoisomers::StereoisomerEnumerator>
      dp_enumerator;
};

}  // namespace
namespace RDKit {

BOOST_PYTHON_MODULE(rdEnumerateStereoisomers) {
  python::scope().attr("__doc__") =
      "Module containing functions to enumerate stereoisomers of a molecule."
      "  Chiral centers and double bonds will be enumerated if unassigned, or,"
      " if the appropriate option is set, if assigned.  Atropisomers will only"
      " be enumerated if assigned.  There is, as yet, no means of finding "
      " unassigned atropisomers.";

  std::string docString = "EnumerateSteroisomers options.";
  python::class_<EnumerateStereoisomers::StereoEnumerationOptions,
                 boost::noncopyable>("StereoEnumerationOptions",
                                     docString.c_str())
      .def(python::init<bool, bool, bool, bool, std::uint64_t, int, bool>(
          (python::args("self"), python::args("tryEmbedding") = false,
           python::args("onlyUnassigned") = true,
           python::args("onlyStereoGroups") = false,
           python::args("unique") = true, python::args("maxIsomers") = 16384,
           python::args("randomSeed") = 0xC0FFEE,
           python::args("useRingSystemFilter") = true)))
      .def_readwrite(
          "tryEmbedding",
          &EnumerateStereoisomers::StereoEnumerationOptions::tryEmbedding,
          "If true, the process attempts to generate a standard RDKit distance geometry"
          " conformation for the stereoisomer.  If this fails, we assume that the stereoisomer is"
          " non-physical and don't return it.  NOTE that this is computationally expensive and is"
          " just a heuristic that could result in stereoisomers being lost.  Default=False")
      .def_readwrite(
          "onlyUnassigned",
          &EnumerateStereoisomers::StereoEnumerationOptions::onlyUnassigned,
          "If true, stereocenters which have a specified stereochemistry will not be"
          " perturbed unless they are part of a relative stereo group.  Default=True.")
      .def_readwrite(
          "onlyStereoGroups",
          &EnumerateStereoisomers::StereoEnumerationOptions::onlyStereoGroups,
          "If true, only find stereoisomers that differ at the StereoGroups associated with"
          " the molecule.  Default=False.")
      .def_readwrite(
          "unique", &EnumerateStereoisomers::StereoEnumerationOptions::unique,
          "If true, only stereoisomers that differ in canonical CXSmiles will be"
          " returned.  Default=True.")
      .def_readwrite(
          "maxIsomers",
          &EnumerateStereoisomers::StereoEnumerationOptions::maxIsomers,
          "The maximum number of isomers to yield.  If the number of possible isomers"
          " is greater than maxIsomers, a random subset will be yielded.  If 0, there"
          " is no maximum.  Since every additional stereocenter doubles the number of"
          " results (and execution time) it's important to keep an eye on this."
          " The default value is 16384.")
      .def_readwrite(
          "randomSeed",
          &EnumerateStereoisomers::StereoEnumerationOptions::randomSeed,
          "Seed for random number generator. -1 means no seed. Default=0xC0FFEE")
      .def_readwrite(
          "useRingSystemFilter",
          &EnumerateStereoisomers::StereoEnumerationOptions::
              useRingSystemFilter,
          "filter out stereoisomers containing ring systems that are known to be"
          " impossible in 3D space.  Default=True.")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(returns a generator that yields possible stereoisomers for a molecule

    Arguments:
      - m: the molecule to work with
      - options: parameters controlling the enumeration

    If m has stereogroups, they will be expanded

    A small example with 3 chiral atoms and 1 chiral bond (16 theoretical stereoisomers):

    >>> from rdkit import Chem
    >>> from rdkit.Chem.rdEnumerateStereoisomers import StereoisomerEnumerator, StereoEnumerationOptions
    >>> m = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    >>> isomers = tuple(StereoisomerEnumerator(m))
    >>> len(isomers)
    8
    >>> for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2

    Because the molecule is constrained, not all of those isomers can
    actually exist. We can check that:

    >>> opts = StereoEnumerationOptions(tryEmbedding=True)
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
    >>> len(isomers)
    8
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2

    Or we can force the output to only give us unique isomers:

    >>> m = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    >>> opts = StereoEnumerationOptions(unique=True)
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
    >>> len(isomers)
    10
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl
    F[C@@H](Cl)/C=C\C=C/[C@@H](F)Cl
    F[C@@H](Cl)/C=C\C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl
    F[C@H](Cl)/C=C/C=C/[C@H](F)Cl
    F[C@H](Cl)/C=C/C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C/[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C/[C@H](F)Cl
    F[C@H](Cl)/C=C\C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C\[C@H](F)Cl

    By default the code only expands unspecified stereocenters:

    >>> m = Chem.MolFromSmiles('BrC=C[C@H]1OC(C2)(F)C2(Cl)C1')
    >>> isomers = tuple(StereoisomerEnumerator(m))
    >>> len(isomers)
    4
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2

    But we can change that behavior:

    >>> opts = StereoEnumerationOptions(onlyUnassigned=False)
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
    >>> len(isomers)
    8

    Since the result is a generator, we can allow exploring at least parts of very
    large result sets:

    >>> m = Chem.MolFromSmiles('Br' + '[CH](Cl)' * 20 + 'F')
    >>> opts = StereoEnumerationOptions(maxIsomers=0)
    >>> isomers = StereoisomerEnumerator(m, options=opts)
    >>> for x in range(5):
    ...   print(Chem.MolToSmiles(next(isomers),isomericSmiles=True))
    F[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)Br

    Or randomly sample a small subset. Note that if we want that sampling to be consistent
    across runs we need to provide a random number seed:

    >>> m = Chem.MolFromSmiles('Br' + '[CH](Cl)' * 20 + 'F')
    >>> opts = StereoEnumerationOptions(maxIsomers=3,randomSeed=0xf00d)
    >>> isomers = StereoisomerEnumerator(m, options=opts)
    >>> for smi in isomers: #sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
    ...     print(Chem.MolToSmiles(smi))
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)Br
  )DOC";
  python::class_<LocalStereoEnumerator, boost::noncopyable>(
      "StereoisomerEnumerator", docString.c_str(), python::no_init)
      .def(python::init<python::object &, python::object &, bool>(
          (python::arg("options") = python::object(),
           python::arg("verbose") = true)))
      .def(
          "__iter__", +[](python::object const &self) { return self; },
          python::args("self"))
      .def(
          "__next__", &LocalStereoEnumerator::nextForIter,
          "Get the next isomer in the sequence. Raises StopIteration at the end.")
      .def("next", &LocalStereoEnumerator::next,
           "Get next isomer in the sequence, or None if at the end.")
      .def("GetStereoisomerCount", &LocalStereoEnumerator::GetStereoisomerCount,
           "Get the number of stereoisomers.");
}

}  // namespace RDKit
