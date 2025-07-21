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

  unsigned int GetStereoisomerCount() {
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
          " results (and execution time) it's important to keep an eye on this.")
      .def_readwrite(
          "randomSeed",
          &EnumerateStereoisomers::StereoEnumerationOptions::randomSeed,
          "Seed for random number generator.  Default=-1 means no seed.")
      .def("__setattr__", &safeSetattr);

  docString = "Stereoisomer enumerator.";
  python::class_<LocalStereoEnumerator, boost::noncopyable>(
      "StereoisomerEnumerator", docString.c_str(), python::no_init)
      .def(python::init<python::object &, python::object &, bool>(
          (python::arg("options") = python::object(),
           python::arg("verbose") = true)))
      .def("next", &LocalStereoEnumerator::next,
           "Get next isomer in the sequence, or None if at the end.")
      .def("GetStereoisomerCount", &LocalStereoEnumerator::GetStereoisomerCount,
           "Get the number of stereoisomers.");
}

}  // namespace RDKit