//
//  Copyright (c) 2015-2026 Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <RDBoost/Wrap_nb.h>
#include <nanobind/stl/shared_ptr.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSample.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSampleAllBBs.h>
#include <GraphMol/ChemReactions/Enumerate/EvenSamplePairs.h>
#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>
#include <cstdint>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

std::vector<RDKit::MOL_SPTR_VECT> ConvertToVect(nb::object bbs) {
  std::vector<RDKit::MOL_SPTR_VECT> vect;
  for (nb::handle row_handle : bbs) {
    RDKit::MOL_SPTR_VECT reacts;
    for (nb::handle mol_handle : row_handle) {
      reacts.push_back(RDKit::ROMOL_SPTR(nb::cast<RDKit::ROMol *>(mol_handle), [](RDKit::ROMol *) {}));
    }
    vect.push_back(std::move(reacts));
  }
  return vect;
}

nb::tuple EnumerateLibraryBase__next__(RDKit::EnumerateLibraryBase *base) {
  if (!static_cast<bool>(*base)) {
    PyErr_SetString(PyExc_StopIteration, "Enumerations exhausted");
    throw nb::python_error();
  }
  std::vector<RDKit::MOL_SPTR_VECT> mols;
  {
    NOGIL gil;
    mols = base->next();
  }
  nb::list res;
  for (const auto &mol_vec : mols) {
    nb::list inner;
    for (const auto &mol : mol_vec) {
      inner.append(std::shared_ptr<RDKit::ROMol>(mol.get(), [b = mol](RDKit::ROMol *) {}));
    }
    res.append(nb::tuple(inner));
  }
  return nb::tuple(res);
}

nb::bytes EnumerateLibraryBase_Serialize(const EnumerateLibraryBase &en) {
  const auto res = en.Serialize();
  return nb::bytes(res.c_str(), res.size());
}

void ToBBS(EnumerationStrategyBase &rgroup, ChemicalReaction &rxn,
           nb::object ob) {
  rgroup.initialize(rxn, ConvertToVect(ob));
}

}  // namespace RDKit

void wrap_enumeration(nb::module_ &m) {
  using namespace RDKit;

  nb::class_<RDKit::EnumerateLibraryBase>(m, "EnumerateLibraryBase")
      .def("__bool__",
           [](RDKit::EnumerateLibraryBase *base) {
             return static_cast<bool>(*base);
           })
      .def("__iter__", [](nb::object self) { return self; })
      .def("next", &EnumerateLibraryBase__next__,
           "Return the next molecule from the enumeration.")
      .def("__next__", &EnumerateLibraryBase__next__,
           "Return the next molecule from the enumeration.")
      .def("nextSmiles", &RDKit::EnumerateLibraryBase::nextSmiles,
           "Return the next smiles string from the enumeration.")
      .def("Serialize", &EnumerateLibraryBase_Serialize,
           R"DOC(Serialize the library to a binary string.
Note that the position in the library is serialized as well.  Care should
be taken when serializing.  See GetState/SetState for position manipulation.)DOC")
      .def("InitFromString", &RDKit::EnumerateLibraryBase::initFromString,
           "data"_a, "Initialize the library from a binary string")
      .def("InitFromString",
           [](RDKit::EnumerateLibraryBase &self, nb::bytes data) {
             self.initFromString(std::string(
                 static_cast<const char *>(data.data()), data.size()));
           },
           "data"_a, "Initialize the library from a binary string")
      .def(
          "GetPosition", &RDKit::EnumerateLibraryBase::getPosition,
          nb::rv_policy::reference_internal,
          R"DOC(Returns the current enumeration position into the reagent vectors, as
returned by GetReagents().  They do not necessarily refer to
the input reagent sets as they only refer to reagents compatible
with the reaction.)DOC")
      .def(
          "GetState", &RDKit::EnumerateLibraryBase::getState,
          R"DOC(Returns the current enumeration state (position) of the library.
This position can be used to restart the library from a known position)DOC")
      .def("SetState", &RDKit::EnumerateLibraryBase::setState, "state"_a,
           "Sets the enumeration state (position) of the library.")
      .def("ResetState", &RDKit::EnumerateLibraryBase::resetState,
           "Returns the current enumeration state (position) of the library "
           "to the start.")
      .def("GetReaction", &RDKit::EnumerateLibraryBase::getReaction,
           nb::rv_policy::reference_internal,
           "Returns the chemical reaction for this library")
      .def("GetEnumerator", &RDKit::EnumerateLibraryBase::getEnumerator,
           nb::rv_policy::reference_internal,
           "Returns the enumation strategy for the current library");

  nb::class_<RDKit::EnumerationParams>(m, "EnumerationParams",
                                       R"DOC(EnumerationParams
Controls some aspects of how the enumeration is performed.
Options:
  reagentMaxMatchCount [ default Infinite ]
    This specifies how many times the reactant template can match a reagent.

  sanePartialProducts [default false]
    If true, forces all products of the reagent plus the product templates
     pass chemical sanitization.  Note that if the product template itself
     does not pass sanitization, then none of the products will.
)DOC")
      .def(nb::init<>())
      .def_rw("reagentMaxMatchCount",
              &RDKit::EnumerationParams::reagentMaxMatchCount)
      .def_rw("sanePartialProducts",
              &RDKit::EnumerationParams::sanePartialProducts);

  nb::class_<RDKit::EnumerateLibrary, RDKit::EnumerateLibraryBase>(
      m, "EnumerateLibrary",
      R"DOC(EnumerateLibrary
This class allows easy enumeration of reactions.  Simply provide a reaction
and a set of reagents and you are off the races.

Note that this functionality should be considered beta and that the API may
change in a future release.

EnumerateLibrary follows the python enumerator protocol, for example:

library = EnumerateLibrary(rxn, bbs)
for products in library:
   ... do something with the product

It is useful to sanitize reactions before hand:

SanitizeRxn(rxn)
library = EnumerateLibrary(rxn, bbs)

If ChemDraw style reaction semantics are prefereed, you can apply
the ChemDraw parameters:

SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())

For one, this enforces only matching RGroups and assumes all atoms
have fully satisfied valences.

Each product has the same output as applying a set of reagents to
the libraries reaction.

This can be a bit confusing as each product can have multiple molecules
generated.  The returned data structure is as follows:

   [ [products1], [products2],... ]
Where products1 are the molecule products for the reactions first product
template and products2 are the molecule products for the second product
template.  Since each reactant can match more than once, there may be
multiple product molecules for each template.

for products in library:
    for results_for_product_template in products:
        for mol in results_for_product_template:
            Chem.MolToSmiles(mol) # finally have a molecule!

For sufficiently large libraries, using this iteration strategy is not
recommended as the library may contain more products than atoms in the
universe.  To help with this, you can supply an enumeration strategy.
The default strategy is a CartesianProductStrategy which enumerates
everything.  RandomSampleStrategy randomly samples the products but
this strategy never terminates, however, python supplies itertools:

import itertools
library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
for result in itertools.islice(library, 1000):
    # do something with the first 1000 samples

for result in itertools.islice(library, 1000):
    # do something with the next 1000 samples

Libraries are also serializable, including their current state:

s = library.Serialize()
library2 = EnumerateLibrary()
library2.InitFromString(s)
for result in itertools.islice(libary2, 1000):
    # do something with the next 1000 samples
)DOC")
      .def(nb::init<>())
      .def("__init__",
           [](RDKit::EnumerateLibrary *self,
              const RDKit::ChemicalReaction &rxn, nb::object bbs,
              const RDKit::EnumerationParams &params) {
             new (self)
                 RDKit::EnumerateLibrary(rxn, RDKit::ConvertToVect(bbs), params);
           },
           "rxn"_a, "reagents"_a, "params"_a = RDKit::EnumerationParams())
      .def("__init__",
           [](RDKit::EnumerateLibrary *self,
              const RDKit::ChemicalReaction &rxn, nb::object bbs,
              const RDKit::EnumerationStrategyBase &enumerator,
              const RDKit::EnumerationParams &params) {
             new (self) RDKit::EnumerateLibrary(
                 rxn, RDKit::ConvertToVect(bbs), enumerator, params);
           },
           "rxn"_a, "reagents"_a, "enumerator"_a,
           "params"_a = RDKit::EnumerationParams())
      .def(
          "GetReagents",
          [](const RDKit::EnumerateLibrary &lib) {
            nb::list outer;
            for (const auto &inner_vec : lib.getReagents()) {
              nb::list inner;
              for (const auto &mol : inner_vec) {
                inner.append(std::shared_ptr<RDKit::ROMol>(mol.get(), [b = mol](RDKit::ROMol *) {}));
              }
              outer.append(inner);
            }
            return outer;
          },
          R"DOC(Return the reagents used in this library.  These are the subset
of the input reagents that are compatible with the reaction so may
be smaller than the input reagent sets.)DOC");

  nb::class_<RDKit::EnumerationStrategyBase>(m, "EnumerationStrategyBase")
      .def("__bool__",
           [](RDKit::EnumerationStrategyBase *base) {
             return static_cast<bool>(*base);
           })
      .def("Type", &RDKit::EnumerationStrategyBase::type,
           "Returns the enumeration strategy type as a string.")
      .def("Skip", &RDKit::EnumerationStrategyBase::skip, "skipCount"_a,
           R"DOC(Skip the next Nth results. note: this may be an expensive operation
depending on the enumeration strategy used. It is recommended to use
the enumerator state to advance to a known position)DOC")
      .def("__copy__",
           [](const RDKit::EnumerationStrategyBase &self) {
             return self.copy();
           },
           nb::rv_policy::take_ownership)
      .def(
          "GetNumPermutations",
          &RDKit::EnumerationStrategyBase::getNumPermutations,
          R"DOC(Returns the total number of results for this enumeration strategy.
Note that some strategies are effectively infinite.)DOC")
      .def("GetPosition", &RDKit::EnumerationStrategyBase::getPosition,
           nb::rv_policy::reference_internal,
           R"DOC(Return the current indices into the arrays of reagents, as
returned by GetReagents().  They do not necessarily refer to
the input reagent sets as they only refer to reagents compatible
with the reaction.)DOC")
      .def("next", &RDKit::EnumerationStrategyBase::next,
           nb::rv_policy::reference_internal,
           "Return the next indices into the arrays of reagents")
      .def("__next__", &RDKit::EnumerationStrategyBase::next,
           nb::rv_policy::reference_internal,
           "Return the next indices into the arrays of reagents")
      .def("Initialize", &RDKit::ToBBS, "rxn"_a, "ob"_a);

  nb::class_<RDKit::CartesianProductStrategy,
             RDKit::EnumerationStrategyBase>(
      m, "CartesianProductStrategy",
      R"DOC(CartesianProductStrategy produces a standard walk through all possible
reagent combinations:

(0,0,0), (1,0,0), (2,0,0) ...)DOC")
      .def(nb::init<>())
      .def("__copy__",
           [](const RDKit::CartesianProductStrategy &self) {
             return self.copy();
           },
           nb::rv_policy::take_ownership);

  nb::class_<RDKit::RandomSampleStrategy, RDKit::EnumerationStrategyBase>(
      m, "RandomSampleStrategy",
      "RandomSampleStrategy simply randomly samples from the reagent sets.\n"
      "Note that this strategy never halts and can produce duplicates.")
      .def(nb::init<>())
      .def("__copy__",
           [](const RDKit::RandomSampleStrategy &self) {
             return self.copy();
           },
           nb::rv_policy::take_ownership);

  nb::class_<RDKit::RandomSampleAllBBsStrategy,
             RDKit::EnumerationStrategyBase>(
      m, "RandomSampleAllBBsStrategy",
      R"DOC(RandomSampleAllBBsStrategy randomly samples from the reagent sets
with the constraint that all building blocks are samples as early as possible.
Note that this strategy never halts and can produce duplicates.)DOC")
      .def(nb::init<>())
      .def("__copy__",
           [](const RDKit::RandomSampleAllBBsStrategy &self) {
             return self.copy();
           },
           nb::rv_policy::take_ownership);

  nb::class_<RDKit::EvenSamplePairsStrategy, RDKit::EnumerationStrategyBase>(
      m, "EvenSamplePairsStrategy",
      R"DOC(Randomly sample Pairs evenly from a collection of building blocks
This is a good strategy for choosing a relatively small selection
of building blocks from a larger set.  As the amount of work needed
to retrieve the next evenly sample building block grows with the
number of samples, this method performs progressively worse as the
number of samples gets larger.
See EnumerationStrategyBase for more details.)DOC")
      .def(nb::init<>())
      .def("__copy__",
           [](const RDKit::EvenSamplePairsStrategy &self) {
             return self.copy();
           },
           nb::rv_policy::take_ownership)
      .def("Stats", &RDKit::EvenSamplePairsStrategy::stats,
           "Return the statistics log of the pairs used in the current "
           "enumeration.");

  m.def("EnumerateLibraryCanSerialize", RDKit::EnumerateLibraryCanSerialize,
        "Returns True if the EnumerateLibrary is serializable "
        "(requires boost serialization)");
}
