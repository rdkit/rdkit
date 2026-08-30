//
//  Copyright (C) 2020-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/trampoline.h>
#include <nanobind/make_iterator.h>

#include <optional>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <RDBoost/boost_shared_ptr.h>

NB_MAKE_OPAQUE(std::vector<
               RDKit::MolStandardize::TautomerScoringFunctions::SubstructTerm>);

#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

std::shared_ptr<RDKit::ROMol> toStd(const RDKit::ROMOL_SPTR &bptr) {
  return {bptr.get(), [b = bptr](RDKit::ROMol *) {}};
}

std::vector<std::shared_ptr<RDKit::ROMol>> tERTautomersGetterHelper(
    const RDKit::MolStandardize::TautomerEnumeratorResult &self) {
  auto tautomers = self.tautomers();
  std::vector<std::shared_ptr<RDKit::ROMol>> out;
  out.reserve(tautomers.size());
  for (const auto &sptr : tautomers) {
    out.push_back(toStd(sptr));
  }
  return out;
}

template <typename T>
nb::tuple bitsetToTuple(const boost::dynamic_bitset<T> &bs) {
  nb::list atList;
  for (auto i = bs.find_first(); i != boost::dynamic_bitset<T>::npos;
       i = bs.find_next(i)) {
    atList.append(i);
  }
  return nb::tuple(atList);
}

struct TautomerEnumeratorCallbackTrampoline
    : public MolStandardize::TautomerEnumeratorCallback {
  NB_TRAMPOLINE(MolStandardize::TautomerEnumeratorCallback, 1);

  bool operator()(
      const ROMol &mol,
      const MolStandardize::TautomerEnumeratorResult &res) override {
    NB_OVERRIDE_PURE_NAME("__call__", operator(), mol, res);
  }
};

nb::tuple smilesTautomerMapKeysHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list keys;
  for (const auto &pair : self) {
    keys.append(pair.first);
  }
  return nb::tuple(keys);
}

nb::tuple smilesTautomerMapValuesHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list values;
  for (const auto &pair : self) {
    auto *t = new MolStandardize::Tautomer(pair.second);
    values.append(nb::cast(t, nb::rv_policy::take_ownership));
  }
  return nb::tuple(values);
}

nb::tuple smilesTautomerMapItemsHelper(
    const MolStandardize::SmilesTautomerMap &self) {
  nb::list items;
  for (const auto &pair : self) {
    auto *t = new MolStandardize::Tautomer(pair.second);
    items.append(
        nb::make_tuple(pair.first, nb::cast(t, nb::rv_policy::take_ownership)));
  }
  return nb::tuple(items);
}

nb::object getCallbackHelper(const MolStandardize::TautomerEnumerator &te) {
  auto cppCallback = te.getCallback();
  if (!cppCallback) {
    return nb::none();
  }

  auto pyCallback = nb::cast(cppCallback);
  if (!pyCallback.is_valid()) {
    return nb::none();
  }

  return pyCallback;
}

void setCallbackHelper(MolStandardize::TautomerEnumerator &te,
                       nb::object callback) {
  if (callback.is_none()) {
    te.setCallback(nullptr);
    return;
  }

  std::shared_ptr<MolStandardize::TautomerEnumeratorCallback> cppObj = nullptr;
  if (!nb::try_cast(callback, cppObj)) {
    throw nb::type_error(
        "Expected an instance of a rdMolStandardize.TautomerEnumeratorCallback subclass");
  }

  // Verify that the Python subclass has a properly overridden __call__:
  // - it must be defined in the immediate class dict (not just inherited)
  // - it must be callable
  nb::handle cls(PyObject_Type(callback.ptr()));
  nb::object cls_dict = cls.attr("__dict__");
  cls.dec_ref();

  if (!cls_dict.is_valid() ||
      !nb::cast<bool>(cls_dict.attr("__contains__")("__call__"))) {
    throw nb::attribute_error(
        "TautomerEnumeratorCallback subclass must override __call__");
  }

  nb::object call_attr = cls_dict.attr("__getitem__")("__call__");
  if (!PyCallable_Check(call_attr.ptr())) {
    throw nb::attribute_error(
        "TautomerEnumeratorCallback.__call__ must be callable");
  }

  te.setCallback(cppObj);
}

class pyobjFunctor {
 public:
  pyobjFunctor(nb::object obj) : dp_obj(std::move(obj)) {}
  ~pyobjFunctor() = default;
  int operator()(const ROMol &m) {
    nb::object result = dp_obj(m);
    try {
      return nb::cast<int>(result);
    } catch (const nb::cast_error &) {
      throw nb::type_error("scoring function must return a numeric value");
    }
  }

 private:
  nb::object dp_obj;
};

ROMol *canonicalizeHelper(const MolStandardize::TautomerEnumerator &self,
                          const ROMol &mol) {
  return self.canonicalize(mol);
}

ROMol *canonicalizeHelper2(const MolStandardize::TautomerEnumerator &self,
                           const ROMol &mol, nb::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  return self.canonicalize(mol, ftor);
}

inline std::vector<ROMOL_SPTR> extractPythonIterable(const nb::object &o) {
  if (!nb::hasattr(o, "__iter__")) {
    throw nb::type_error(
        "the passed object should be an iterable of Mol objects");
  }
  std::vector<ROMOL_SPTR> result;
  for (nb::handle h : o) {
    RDKit::ROMol *molPtr = nullptr;
    try {
      molPtr = nb::cast<RDKit::ROMol *>(h);
    } catch (const nb::cast_error &) {
      throw nb::type_error(
          "the passed object should be an iterable of Mol objects");
    }
    // Copy the molecule so the vector owns its lifetime independently of
    // Python
    result.push_back(ROMOL_SPTR(new RDKit::ROMol(*molPtr)));
  }
  return result;
}

ROMol *pickCanonicalHelper(const MolStandardize::TautomerEnumerator &self,
                           const nb::object &o) {
  try {
    auto e = nb::cast<MolStandardize::TautomerEnumeratorResult *>(o);
    return self.pickCanonical(*e);
  } catch (const nb::cast_error &) {
    return self.pickCanonical(extractPythonIterable(o));
  }
}

ROMol *pickCanonicalHelper2(const MolStandardize::TautomerEnumerator &self,
                            const nb::object &o, nb::object scoreFunc) {
  pyobjFunctor ftor(scoreFunc);
  try {
    auto e = nb::cast<MolStandardize::TautomerEnumeratorResult *>(o);
    return self.pickCanonical(*e, ftor);
  } catch (const nb::cast_error &) {
    return self.pickCanonical(extractPythonIterable(o), ftor);
  }
}

MolStandardize::TautomerEnumeratorResult enumerateHelper(
    const MolStandardize::TautomerEnumerator &self, const ROMol &mol) {
  return self.enumerate(mol);
}

std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm>
GetDefaultTautomerSubstructsHelper() {
  std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm> terms;
  for (auto term : MolStandardize::TautomerScoringFunctions::
           getDefaultTautomerScoreSubstructs()) {
    terms.emplace_back(term);
  }
  return terms;
}

}  // namespace

void wrap_tautomer(nb::module_ &m) {
  nb::enum_<MolStandardize::TautomerEnumeratorStatus>(
      m, "TautomerEnumeratorStatus")
      .value("Completed", MolStandardize::TautomerEnumeratorStatus::Completed)
      .value("MaxTautomersReached",
             MolStandardize::TautomerEnumeratorStatus::MaxTautomersReached)
      .value("MaxTransformsReached",
             MolStandardize::TautomerEnumeratorStatus::MaxTransformsReached)
      .value("Canceled", MolStandardize::TautomerEnumeratorStatus::Canceled);

  nb::class_<MolStandardize::TautomerEnumeratorCallback,
             TautomerEnumeratorCallbackTrampoline>(
      m, "TautomerEnumeratorCallback",
      R"DOC(Create a derived class from this abstract base class and
implement the __call__() method.
The __call__() method is called in the innermost loop of the
algorithm, and provides a mechanism to monitor or stop
its progress.

To have your callback called, pass an instance of your
derived class to TautomerEnumerator.SetCallback())DOC")
      .def(nb::init<>())
      .def("__call__", &MolStandardize::TautomerEnumeratorCallback::operator(),
           "mol"_a, "res"_a,
           "This must be implemented in the derived class. "
           "Return True if the tautomer enumeration should continue; "
           "False if the tautomer enumeration should stop.\n");

  nb::class_<MolStandardize::Tautomer>(
      m, "Tautomer",
      "used to hold the aromatic and kekulized versions of each tautomer")
      .def_prop_ro(
          "tautomer",
          [](const MolStandardize::Tautomer &self) {
            return new ROMol(*self.tautomer);
          },
          nb::rv_policy::take_ownership, "aromatic version of the tautomer")
      .def_prop_ro(
          "kekulized",
          [](const MolStandardize::Tautomer &self) {
            return new ROMol(*self.getKekulized());
          },
          nb::rv_policy::take_ownership, "kekulized version of the tautomer");

  nb::class_<MolStandardize::SmilesTautomerMap>(
      m, "SmilesTautomerMap",
      "maps SMILES strings to the respective Tautomer objects")
      .def("keys", &smilesTautomerMapKeysHelper)
      .def("values", &smilesTautomerMapValuesHelper)
      .def("items", &smilesTautomerMapItemsHelper)
      .def("__len__", [](const MolStandardize::SmilesTautomerMap &self) {
        return self.size();
      });

  nb::class_<MolStandardize::TautomerEnumeratorResult>(
      m, "TautomerEnumeratorResult",
      "used to return tautomer enumeration results")
      .def_prop_ro("tautomers", &tERTautomersGetterHelper,
                   "tautomers generated by the enumerator")
      .def_prop_ro(
          "smiles",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return self.smiles();
          },
          "SMILES of tautomers generated by the enumerator")
      .def_prop_ro(
          "smilesTautomerMap",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return self.smilesTautomerMap();
          },
          nb::rv_policy::reference_internal,
          "dictionary mapping SMILES strings to the respective Tautomer objects")
      .def_prop_ro(
          "status",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return self.status();
          },
          "whether the enumeration completed or not; see "
          "TautomerEnumeratorStatus for possible values")
      .def_prop_ro(
          "modifiedAtoms",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return bitsetToTuple(self.modifiedAtoms());
          },
          "tuple of atom indices modified by the transforms")
      .def_prop_ro(
          "modifiedBonds",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return bitsetToTuple(self.modifiedBonds());
          },
          "tuple of bond indices modified by the transforms")
      .def("__call__", &tERTautomersGetterHelper,
           "tautomers generated by the enumerator")
      .def(
          "__iter__",
          [](const MolStandardize::TautomerEnumeratorResult &self) {
            return nb::make_iterator(
                nb::type<MolStandardize::TautomerEnumeratorResult>(),
                "TautomerEnumeratorResult iterator", self.begin(), self.end());
          },
          nb::keep_alive<0, 1>())
      .def(
          "__getitem__",
          [](const MolStandardize::TautomerEnumeratorResult &self, int pos) {
            auto sz = static_cast<int>(self.size());
            if (pos < 0) {
              pos += sz;
            }
            if (pos < 0 || pos >= sz) {
              throw std::out_of_range("index out of bounds");
            }
            return self[pos];
          },
          "pos"_a, nb::rv_policy::take_ownership)
      .def("__len__", [](const MolStandardize::TautomerEnumeratorResult &self) {
        return self.size();
      });

  nb::class_<MolStandardize::TautomerEnumerator>(m, "TautomerEnumerator")
      .def(nb::init<>())
      .def(nb::init<const MolStandardize::CleanupParameters &>(), "params"_a)
      .def(nb::init<const MolStandardize::TautomerEnumerator &>(), "other"_a)
      .def("Enumerate", &enumerateHelper, "mol"_a,
           R"DOC(Generates the tautomers for a molecule.

The enumeration rules are inspired by the publication:
M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
https://doi.org/10.1007/s10822-010-9346-4

Note: the definitions used here are that the atoms modified during
tautomerization are the atoms at the beginning and end of each tautomer
transform (the H "donor" and H "acceptor" in the transform) and the bonds
modified during transformation are any bonds whose order is changed during
the tautomer transform (these are the bonds between the "donor" and the
"acceptor").)DOC")
      .def("Canonicalize", &canonicalizeHelper, "mol"_a,
           nb::rv_policy::take_ownership,
           R"DOC(Returns the canonical tautomer for a molecule.

The default scoring scheme is inspired by the publication:
M. Sitzmann et al., "Tautomerism in Large Databases.", JCAMD 24:521 (2010)
https://doi.org/10.1007/s10822-010-9346-4

Note that the canonical tautomer is very likely not the most stable tautomer
for any given conditions. The default scoring rules are designed to produce
"reasonable" tautomers, but the primary concern is that the results are
canonical: you always get the same canonical tautomer for a molecule
regardless of what the input tautomer or atom ordering were.)DOC")
      .def("Canonicalize", &canonicalizeHelper2, "mol"_a, "scoreFunc"_a,
           nb::rv_policy::take_ownership,
           "picks the canonical tautomer from an iterable of molecules "
           "using a custom scoring function")
      .def("PickCanonical", &pickCanonicalHelper, "iterable"_a,
           nb::rv_policy::take_ownership,
           "picks the canonical tautomer from an iterable of molecules")
      .def("PickCanonical", &pickCanonicalHelper2, "iterable"_a, "scoreFunc"_a,
           nb::rv_policy::take_ownership,
           "returns the canonical tautomer for a molecule using a custom "
           "scoring function")
      .def_static(
          "ScoreTautomer",
          &MolStandardize::TautomerScoringFunctions::scoreTautomer, "mol"_a,
          "returns the score for a tautomer using the default scoring scheme.")
      .def("SetMaxTautomers",
           &MolStandardize::TautomerEnumerator::setMaxTautomers,
           "maxTautomers"_a,
           "set the maximum number of tautomers to be generated.")
      .def("GetMaxTautomers",
           &MolStandardize::TautomerEnumerator::getMaxTautomers,
           "returns the maximum number of tautomers to be generated.")
      .def("SetMaxTransforms",
           &MolStandardize::TautomerEnumerator::setMaxTransforms,
           "maxTransforms"_a,
           "set the maximum number of transformations to be applied. "
           "This limit is usually hit earlier than the maxTautomers limit "
           "and leads to a more linear scaling of CPU time with increasing "
           "number of tautomeric centers (see Sitzmann et al.).")
      .def("GetMaxTransforms",
           &MolStandardize::TautomerEnumerator::getMaxTransforms,
           "returns the maximum number of transformations to be applied.")
      .def("SetRemoveSp3Stereo",
           &MolStandardize::TautomerEnumerator::setRemoveSp3Stereo,
           "removeSp3Stereo"_a,
           "set to True if you wish stereochemistry information "
           "to be removed from sp3 atoms involved in tautomerism. "
           "This means that S-aminoacids will lose their stereochemistry "
           "after going through tautomer enumeration because of the "
           "amido-imidol tautomerism. This defaults to True in RDKit, "
           "and to False in the workflow described by Sitzmann et al.")
      .def("GetRemoveSp3Stereo",
           &MolStandardize::TautomerEnumerator::getRemoveSp3Stereo,
           "returns whether stereochemistry information will be removed from "
           "sp3 atoms involved in tautomerism.")
      .def("SetRemoveBondStereo",
           &MolStandardize::TautomerEnumerator::setRemoveBondStereo,
           "removeBondStereo"_a,
           "set to True if you wish stereochemistry information "
           "to be removed from double bonds involved in tautomerism. "
           "This means that enols will lose their E/Z stereochemistry "
           "after going through tautomer enumeration because of the "
           "keto-enolic tautomerism. This defaults to True in the "
           "RDKit and also in the workflow described by Sitzmann et al.")
      .def("GetRemoveBondStereo",
           &MolStandardize::TautomerEnumerator::getRemoveBondStereo,
           "returns whether stereochemistry information "
           "will be removed from double bonds involved in tautomerism.")
      .def("SetReassignStereo",
           &MolStandardize::TautomerEnumerator::setReassignStereo,
           "reassignStereo"_a,
           "set to True if you wish AssignStereochemistry to be called "
           "on each tautomer generated by the Enumerate() method. "
           "This defaults to True.")
      .def("GetReassignStereo",
           &MolStandardize::TautomerEnumerator::getReassignStereo,
           "returns whether AssignStereochemistry will be called "
           "on each tautomer generated by the Enumerate() method.")
      .def("SetCallback", &setCallbackHelper, "callback"_a,
           "Pass an instance of a class derived from\n"
           "TautomerEnumeratorCallback, which must implement the\n"
           "__call__() method.")
      .def("GetCallback", &getCallbackHelper,
           "Get the TautomerEnumeratorCallback subclass instance,\n"
           "or None if none was set.")
      .def_prop_ro_static("tautomerScoreVersion", [](nb::handle) {
        return MolStandardize::TautomerScoringFunctions::tautomerScoringVersion;
      });

  m.def("GetV1TautomerEnumerator", MolStandardize::getV1TautomerEnumerator,
        "return a TautomerEnumerator using v1 of the enumeration rules",
        nb::rv_policy::take_ownership);

  m.def("ScoreRings", MolStandardize::TautomerScoringFunctions::scoreRings,
        "mol"_a,
        "scores the ring system of the tautomer for canonicalization\n"
        "Aromatic rings score 100, all carbon aromatic rings score 250");

  m.def(
      "ScoreHeteroHs", MolStandardize::TautomerScoringFunctions::scoreHeteroHs,
      "mol"_a,
      "scores the number of heteroHs of the tautomer for canonicalization\n"
      "This gives a negative penalty to hydrogens attached to S,P, Se and Te");

  nb::class_<MolStandardize::TautomerScoringFunctions::SubstructTerm>(
      m, "SubstructTerm",
      "Sets the score of this particular tautomer substructure, higher scores "
      "are more preferable\n"
      "Aromatic rings score 100, all carbon aromatic rings score 250")
      .def(nb::init<std::string, std::string, int>(), "name"_a, "smarts"_a,
           "score"_a)
      .def_ro("name",
              &MolStandardize::TautomerScoringFunctions::SubstructTerm::name)
      .def_ro("smarts",
              &MolStandardize::TautomerScoringFunctions::SubstructTerm::smarts)
      .def_ro("score",
              &MolStandardize::TautomerScoringFunctions::SubstructTerm::score);

  nb::bind_vector<
      std::vector<MolStandardize::TautomerScoringFunctions::SubstructTerm>>(
      m, "SubstructTermVector");

  m.def(
      "ScoreSubstructs",
      [](const ROMol &mol,
         const std::optional<std::vector<
             MolStandardize::TautomerScoringFunctions::SubstructTerm>> &terms)
          -> int {
        if (!terms) {
          return MolStandardize::TautomerScoringFunctions::scoreSubstructs(mol);
        }
        return MolStandardize::TautomerScoringFunctions::scoreSubstructs(
            mol, *terms);
      },
      "mol"_a, "terms"_a = nb::none(), "scores the tautomer substructures");

  m.def("GetDefaultTautomerScoreSubstructs", GetDefaultTautomerSubstructsHelper,
        "Return the default tautomer substructure scoring terms");
}
