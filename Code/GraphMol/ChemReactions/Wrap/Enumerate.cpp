//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutues for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSample.h>
#include <GraphMol/ChemReactions/Enumerate/RandomSampleAllBBs.h>
#include <GraphMol/ChemReactions/Enumerate/EvenSamplePairs.h>
#include <GraphMol/ChemReactions/Enumerate/Enumerate.h>
#include <boost/python/stl_iterator.hpp>
#include <cstdint>

namespace python = boost::python;

namespace RDKit {

template <class T>
std::vector<RDKit::MOL_SPTR_VECT> ConvertToVect(T bbs) {
  std::vector<RDKit::MOL_SPTR_VECT> vect;
  size_t num_bbs = python::extract<unsigned int>(bbs.attr("__len__")());
  vect.resize(num_bbs);
  for (size_t i = 0; i < num_bbs; ++i) {
    unsigned int len1 = python::extract<unsigned int>(bbs[i].attr("__len__")());
    RDKit::MOL_SPTR_VECT &reacts = vect[i];
    reacts.reserve(len1);
    for (unsigned int j = 0; j < len1; ++j) {
      RDKit::ROMOL_SPTR mol = python::extract<RDKit::ROMOL_SPTR>(bbs[i][j]);
      if (mol) {
        reacts.push_back(mol);
      } else {
        throw_value_error("reaction called with non molecule reactant");
      }
    }
  }
  return vect;
}

bool EnumerateLibraryBase__nonzero__(RDKit::EnumerateLibraryBase *base) {
  return static_cast<bool>(*base);
}
bool EnumerationStrategyBase__nonzero__(RDKit::EnumerationStrategyBase *base) {
  return static_cast<bool>(*base);
}

inline python::object pass_through(python::object const &o) { return o; }

PyObject *EnumerateLibraryBase__next__(RDKit::EnumerateLibraryBase *base) {
  if (!static_cast<bool>(*base)) {
    PyErr_SetString(PyExc_StopIteration, "Enumerations exhausted");
    boost::python::throw_error_already_set();
  }
  std::vector<RDKit::MOL_SPTR_VECT> mols;
  {
    NOGIL gil;
    mols = base->next();
  }
  PyObject *res = PyTuple_New(mols.size());

  for (unsigned int i = 0; i < mols.size(); ++i) {
    PyObject *lTpl = PyTuple_New(mols[i].size());
    for (unsigned int j = 0; j < mols[i].size(); ++j) {
      PyTuple_SetItem(lTpl, j,
                      python::converter::shared_ptr_to_python(mols[i][j]));
    }
    PyTuple_SetItem(res, i, lTpl);
  }
  return res;
}

python::object EnumerateLibraryBase_Serialize(const EnumerateLibraryBase &en) {
  std::string res = en.Serialize();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

class EnumerateLibraryWrap : public RDKit::EnumerateLibrary {
 public:
  ~EnumerateLibraryWrap() override {}
  EnumerateLibraryWrap() : RDKit::EnumerateLibrary() {}
  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn, python::list ob,
                       const EnumerationParams &params = EnumerationParams())
      : RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), params) {}

  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn, python::tuple ob,
                       const EnumerationParams &params = EnumerationParams())
      : RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), params) {}

  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn, python::list ob,
                       const EnumerationStrategyBase &enumerator,
                       const EnumerationParams &params = EnumerationParams())
      : RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), enumerator, params) {}

  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn, python::tuple ob,
                       const EnumerationStrategyBase &enumerator,
                       const EnumerationParams &params = EnumerationParams())
      : RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), enumerator, params) {}
};

namespace {
template <typename T>
inline std::vector<T> to_std_vector(const python::object &iterable) {
  return std::vector<T>(python::stl_input_iterator<T>(iterable),
                        python::stl_input_iterator<T>());
}
}  // namespace

void ToBBS(EnumerationStrategyBase &rgroup, ChemicalReaction &rxn,
           python::list ob) {
  rgroup.initialize(rxn, ConvertToVect(ob));
}

struct enumeration_wrapper {
  static void wrap() {
    std::string docString;

    RegisterVectorConverter<std::vector<std::string>>("VectorOfStringVectors");
    RegisterVectorConverter<boost::uint64_t>("VectSizeT");
    RegisterVectorConverter<MOL_SPTR_VECT>("VectMolVect");

    python::class_<RDKit::EnumerateLibraryBase, RDKit::EnumerateLibraryBase *,
                   RDKit::EnumerateLibraryBase &, boost::noncopyable>(
        "EnumerateLibraryBase", python::no_init)
        .def("__nonzero__", &EnumerateLibraryBase__nonzero__)
        .def("__bool__", &EnumerateLibraryBase__nonzero__)
        .def("__iter__", &pass_through)
        .def("next", &EnumerateLibraryBase__next__,
             "Return the next molecule from the enumeration.")
        .def("__next__", &EnumerateLibraryBase__next__,
             "Return the next molecule from the enumeration.")
        .def("nextSmiles", &RDKit::EnumerateLibraryBase::nextSmiles,
             "Return the next smiles string from the enumeration.")
        .def("Serialize", &EnumerateLibraryBase_Serialize,
             "Serialize the library to a binary string.\n"
             "Note that the position in the library is serialized as well.  "
             "Care should\n"
             "be taken when serializing.  See GetState/SetState for position "
             "manipulation.")
        .def("InitFromString", &RDKit::EnumerateLibraryBase::initFromString,
             python::arg("data"),
             "Inititialize the library from a binary string")
        .def(
            "GetPosition", &RDKit::EnumerateLibraryBase::getPosition,
            "Returns the current enumeration position into the reagent vectors",
            python::return_internal_reference<
                1, python::with_custodian_and_ward_postcall<0, 1>>())
        .def(
            "GetState", &RDKit::EnumerateLibraryBase::getState,
            "Returns the current enumeration state (position) of the library.\n"
            "This position can be used to restart the library from a known "
            "position")
        .def("SetState", &RDKit::EnumerateLibraryBase::setState,
             python::arg("state"),
             "Sets the enumeration state (position) of the library.")
        .def("ResetState", &RDKit::EnumerateLibraryBase::resetState,
             "Returns the current enumeration state (position) of the library "
             "to the start.")
        .def("GetReaction", &RDKit::EnumerateLibraryBase::getReaction,
             "Returns the chemical reaction for this library",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())
        .def("GetEnumerator", &RDKit::EnumerateLibraryBase::getEnumerator,
             "Returns the enumation strategy for the current library",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>());

    docString =
        "EnumerationParams\n\
Controls some aspects of how the enumeration is performed.\n\
Options:\n\
  reagentMaxMatchCount [ default Infinite ]\n\
    This specifies how many times the reactant template can match a reagent.\n\
\n\
  sanePartialProducts [default false]\n\
    If true, forces all products of the reagent plus the product templates\n\
     pass chemical sanitization.  Note that if the product template itself\n\
     does not pass sanitization, then none of the products will.\n\
";

    python::class_<RDKit::EnumerationParams, RDKit::EnumerationParams *,
                   RDKit::EnumerationParams &>(
        "EnumerationParams", docString.c_str(), python::init<>())
        .def_readwrite("reagentMaxMatchCount",
                       &RDKit::EnumerationParams::reagentMaxMatchCount)
        .def_readwrite("sanePartialProducts",
                       &RDKit::EnumerationParams::sanePartialProducts);

    docString =
        "EnumerateLibrary\n\
This class allows easy enumeration of reactions.  Simply provide a reaction\n\
and a set of reagents and you are off the races.\n\
\n\
Note that this functionality should be considered beta and that the API may\n\
change in a future release.\n\
\n\
EnumerateLibrary follows the python enumerator protocol, for example:\n\
\n\
library = EnumerateLibrary(rxn, bbs)\n\
for products in library:\n\
   ... do something with the product\n\
\n\
It is useful to sanitize reactions before hand:\n\
\n\
SanitizeRxn(rxn)\n\
library = EnumerateLibrary(rxn, bbs)\n\
\n\
If ChemDraw style reaction semantics are prefereed, you can apply\n\
the ChemDraw parameters:\n\
\n\
SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())\n\
\n\
For one, this enforces only matching RGroups and assumes all atoms\n\
have fully satisfied valences.\n\
\n\
Each product has the same output as applying a set of reagents to\n\
the libraries reaction.\n\
\n\
This can be a bit confusing as each product can have multiple molecules\n\
generated.  The returned data structure is as follows:\n\
\n\
   [ [products1], [products2],... ]\n\
Where products1 are the molecule products for the reactions first product\n\
template and products2 are the molecule products for the second product\n\
template.  Since each reactant can match more than once, there may be\n\
multiple product molecules for each template.\n\
\n\
for products in library:\n\
    for results_for_product_template in products:\n\
        for mol in results_for_product_template:\n\
            Chem.MolToSmiles(mol) # finally have a molecule!\n\
\n\
For sufficiently large libraries, using this iteration strategy is not\n\
recommended as the library may contain more products than atoms in the\n\
universe.  To help with this, you can supply an enumeration strategy.\n\
The default strategy is a CartesianProductStrategy which enumerates\n\
everything.  RandomSampleStrategy randomly samples the products but\n\
this strategy never terminates, however, python supplies itertools:\n\
\n\
import itertools\n\
library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())\n\
for result in itertools.islice(library, 1000):\n\
    # do something with the first 1000 samples\n\
\n\
for result in itertools.islice(library, 1000):\n\
    # do something with the next 1000 samples\n\
\n\
Libraries are also serializable, including their current state:\n\
\n\
s = library.Serialize()\n\
library2 = EnumerateLibrary()\n\
library2.InitFromString(s)\n\
for result in itertools.islice(libary2, 1000):\n\
    # do something with the next 1000 samples\n\
";
    python::class_<EnumerateLibraryWrap, boost::noncopyable,
                   python::bases<RDKit::EnumerateLibraryBase>>(
        "EnumerateLibrary", docString.c_str(), python::init<>())
        .def(python::init<const RDKit::ChemicalReaction &, python::list,
                          python::optional<const RDKit::EnumerationParams &>>(
            python::args("rxn", "reagents", "params")))
        .def(python::init<const RDKit::ChemicalReaction &, python::tuple,
                          python::optional<const RDKit::EnumerationParams &>>(
            python::args("rxn", "reagents", "params")))

        .def(python::init<const RDKit::ChemicalReaction &, python::list,
                          const RDKit::EnumerationStrategyBase &,
                          python::optional<const RDKit::EnumerationParams &>>(
            python::args("rxn", "reagents", "enumerator", "params")))
        .def(python::init<const RDKit::ChemicalReaction &, python::tuple,
                          const RDKit::EnumerationStrategyBase &,
                          python::optional<const RDKit::EnumerationParams &>>(
            python::args("rxn", "reagents", "enumerator", "params")))

        .def("GetReagents", &RDKit::EnumerateLibrary::getReagents,
             "Return the reagents used in this library.",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>());

    // iterator_wrappers<EnumerateLibrary>().wrap("EnumerateLibraryIterator");

    python::class_<RDKit::EnumerationStrategyBase,
                   RDKit::EnumerationStrategyBase *,
                   RDKit::EnumerationStrategyBase &, boost::noncopyable>(
        "EnumerationStrategyBase", python::no_init)
        .def("__nonzero__", &EnumerationStrategyBase__nonzero__)
        .def("__bool__", &EnumerationStrategyBase__nonzero__)
        .def("Type", &EnumerationStrategyBase::type,
             "Returns the enumeration strategy type as a string.")
        .def("Skip", &EnumerationStrategyBase::skip, python::args("skipCount"),
             "Skip the next Nth results. note: this may be an expensive "
             "operation\n"
             "depending on the enumeration strategy used. It is recommended to "
             "use\n"
             "the enumerator state to advance to a known position")
        .def("__copy__", python::pure_virtual(&EnumerationStrategyBase::copy),
             python::return_value_policy<python::manage_new_object>())
        .def("GetNumPermutations", &EnumerationStrategyBase::getNumPermutations,
             "Returns the total number of results for this enumeration "
             "strategy.\n"
             "Note that some strategies are effectively infinite.")
        .def("GetPosition", &EnumerationStrategyBase::getPosition,
             "Return the current indices into the arrays of reagents",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())
        .def("next", python::pure_virtual(&EnumerationStrategyBase::next),
             "Return the next indices into the arrays of reagents",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())
        .def("__next__", python::pure_virtual(&EnumerationStrategyBase::next),
             "Return the next indices into the arrays of reagents",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1>>())
        .def("Initialize", ToBBS);

    docString =
        "CartesianProductStrategy produces a standard walk through all "
        "possible\n"
        "reagent combinations:\n"
        "\n"
        "(0,0,0), (1,0,0), (2,0,0) ...\n";

    python::class_<RDKit::CartesianProductStrategy,
                   RDKit::CartesianProductStrategy *,
                   RDKit::CartesianProductStrategy &,
                   python::bases<EnumerationStrategyBase>>(
        "CartesianProductStrategy", docString.c_str(), python::init<>())
        .def("__copy__", &RDKit::CartesianProductStrategy::copy,
             python::return_value_policy<python::manage_new_object>());

    docString =
        "RandomSampleStrategy simply randomly samples from the reagent sets.\n"
        "Note that this strategy never halts and can produce duplicates.";
    python::class_<RDKit::RandomSampleStrategy, RDKit::RandomSampleStrategy *,
                   RDKit::RandomSampleStrategy &,
                   python::bases<EnumerationStrategyBase>>(
        "RandomSampleStrategy", docString.c_str(), python::init<>())
        .def("__copy__", &RDKit::RandomSampleStrategy::copy,
             python::return_value_policy<python::manage_new_object>());

    docString =
        "RandomSampleAllBBsStrategy randomly samples from the reagent sets\n"
        "with the constraint that all building blocks are samples as early as "
        "possible.\n"
        "Note that this strategy never halts and can produce duplicates.";
    python::class_<RDKit::RandomSampleAllBBsStrategy,
                   RDKit::RandomSampleAllBBsStrategy *,
                   RDKit::RandomSampleAllBBsStrategy &,
                   python::bases<EnumerationStrategyBase>>(
        "RandomSampleAllBBsStrategy", docString.c_str(), python::init<>())
        .def("__copy__", &RDKit::RandomSampleAllBBsStrategy::copy,
             python::return_value_policy<python::manage_new_object>());

    docString =
        "Randomly sample Pairs evenly from a collection of building blocks\n"
        "This is a good strategy for choosing a relatively small selection\n"
        "of building blocks from a larger set.  As the amount of work needed\n"
        "to retrieve the next evenly sample building block grows with the\n"
        "number of samples, this method performs progressively worse as the\n"
        "number of samples gets larger.\n"
        "See EnumerationStrategyBase for more details.\n";

    python::class_<RDKit::EvenSamplePairsStrategy,
                   RDKit::EvenSamplePairsStrategy *,
                   RDKit::EvenSamplePairsStrategy &,
                   python::bases<EnumerationStrategyBase>>(
        "EvenSamplePairsStrategy", docString.c_str(), python::init<>())
        .def("__copy__", &RDKit::EvenSamplePairsStrategy::copy,
             python::return_value_policy<python::manage_new_object>())
        .def("Stats", &RDKit::EvenSamplePairsStrategy::stats,
             "Return the statistics log of the pairs used in the current "
             "enumeration.");

    python::def("EnumerateLibraryCanSerialize", EnumerateLibraryCanSerialize,
                "Returns True if the EnumerateLibrary is serializable "
                "(requires boost serialization");
  }
};

}  // namespace RDKit

void wrap_enumeration() { RDKit::enumeration_wrapper::wrap(); }
