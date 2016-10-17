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
//       products derived from this software without specific prior written permission.
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
#include <GraphMol/Chemreactions/Enumerate/Enumerate.h>
#include <boost/python/stl_iterator.hpp>

namespace python = boost::python;


namespace RDKit {
  
template<class T>
std::vector<RDKit::MOL_SPTR_VECT> ConvertToVect(T bbs) {
  std::vector<RDKit::MOL_SPTR_VECT> vect;
  size_t num_bbs = python::extract<unsigned int>(bbs.attr("__len__")());
  vect.resize(num_bbs);
  for(size_t i=0; i<num_bbs; ++i) {
    unsigned int len1 = python::extract<unsigned int>(bbs[i].attr("__len__")());
    RDKit::MOL_SPTR_VECT &reacts = vect[i];
    reacts.reserve(len1);
    for(unsigned int j=0;j<len1;++j){
      RDKit::ROMOL_SPTR mol = python::extract<RDKit::ROMOL_SPTR>(bbs[i][j]);
      if(mol)
        reacts.push_back(mol);
      else {
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
RDKit::EnumerateLibraryBase &EnumerateLibraryBase__iter__(
           RDKit::EnumerateLibraryBase *base) {
  return *base;
}

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
  PyObject *res=PyTuple_New(mols.size());
  
  for(unsigned int i=0;i<mols.size();++i){
    PyObject *lTpl =PyTuple_New(mols[i].size());
    for(unsigned int j=0;j<mols[i].size();++j){
      PyTuple_SetItem(lTpl,j,
                      python::converter::shared_ptr_to_python(mols[i][j]));
    }
    PyTuple_SetItem(res,i,lTpl);
  }
  return res;
}

  
class EnumerateLibraryWrap : public RDKit::EnumerateLibrary {
public:
  EnumerateLibraryWrap() : RDKit::EnumerateLibrary() {}
  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn,
                       python::list ob) :
    RDKit::EnumerateLibrary(rxn, ConvertToVect(ob)) {
  }
  
  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn,
                       python::tuple ob) :
    RDKit::EnumerateLibrary(rxn, ConvertToVect(ob)) {
  }
  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn,
                       python::list ob,
                       const EnumerationStrategyBase &enumerator) :
    RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), enumerator) {
  }
  
  EnumerateLibraryWrap(const RDKit::ChemicalReaction &rxn,
                       python::tuple ob,
                       const EnumerationStrategyBase &enumerator) :
    RDKit::EnumerateLibrary(rxn, ConvertToVect(ob), enumerator) {
  }
};

namespace {
  template< typename T >
  inline
  std::vector< T > to_std_vector( const python::object& iterable )
  {
    return std::vector< T >( python::stl_input_iterator< T >( iterable ),
                             python::stl_input_iterator< T >( ) );
  }
}

void ToBBS(EnumerationStrategyBase &rgroup, ChemicalReaction &rxn, python::list ob) {
  rgroup.initialize(rxn, ConvertToVect(ob));
}
  
typedef std::vector<size_t> VectSizeT;
typedef std::vector<std::vector<std::string> > VectStringVect;

struct enumeration_wrapper {
  static void wrap() {
    std::string docString;
    python::class_<VectStringVect>("VectorOfStringVectors")
      .def(python::vector_indexing_suite<VectStringVect, false>() );

    python::class_<VectSizeT>("VectSizeT")
      .def(python::vector_indexing_suite<VectSizeT, false>() );

    python::class_<RDKit::RGroupPosition, RDKit::RGroupPosition *,
                   RDKit::RGroupPosition &>(
        "RGroupPosition",
        "class RGroupPosition.\n"
        " This class holds the current group indices (pos)"
        " and the enumeration state - which can be used to restart the "
        "                             enumeration from this position.",
        python::no_init)
        .def_readonly("pos", &RDKit::RGroupPosition::pos)
        .def_readonly("state", &RDKit::RGroupPosition::state);

    python::class_<RDKit::EnumerateLibraryBase, RDKit::EnumerateLibraryBase *,
                   RDKit::EnumerateLibraryBase &, boost::noncopyable>(
        "EnumerateLibraryBase", python::no_init)
        .def("__nonzero__", &EnumerateLibraryBase__nonzero__)
        .def("__iter__", &EnumerateLibraryBase__iter__,
             python::return_internal_reference<>())
        .def("next", &EnumerateLibraryBase__next__,
             "Return the next molecule from the enumeration.")
        .def("nextSmiles", &RDKit::EnumerateLibraryBase::nextSmiles,
             "Return the next smiles string from the enumeration.")
        .def("Serialize", &RDKit::EnumerateLibraryBase::Serialize,
             "Serialize the library to a binary string.\n"
             "Note that the position in the library is serialized as well.  Care should\n"
             "be taken when serializing.  See GetState/SetState for position manipulation.")
        .def("InitFromString", &RDKit::EnumerateLibraryBase::initFromString,
             python::arg("data"),
             "Inititialize the library from a binary string")
        .def("GetState", &RDKit::EnumerateLibraryBase::getState,
             "Returns the current enumeration state (position) of the library.\n"
             "This position can be used to restart the library from a known position")
        .def("SetState", &RDKit::EnumerateLibraryBase::setState,
             python::arg("state"),
             "Sets the enumeration state (position) of the library.")
        .def("GetEnumerator", &RDKit::EnumerateLibraryBase::getEnumerator,
             "Returns the enumation strategy for the current library",
             python::return_internal_reference<
             1, python::with_custodian_and_ward_postcall<0, 1> >());

    python::class_<EnumerateLibraryWrap,
                   EnumerateLibraryWrap*,EnumerateLibraryWrap&,
                   python::bases<RDKit::EnumerateLibraryBase> >("EnumerateLibrary", "foo",
                                                                python::init<>())
      .def(python::init<const RDKit::ChemicalReaction &,
           python::list>())
      .def(python::init<const RDKit::ChemicalReaction &,
           python::tuple>())
      .def(python::init<const RDKit::ChemicalReaction &,
           python::list,
           const RDKit::EnumerationStrategyBase &>())
      .def(python::init<const RDKit::ChemicalReaction &,
           python::tuple,
           const RDKit::EnumerationStrategyBase &>())
      ;

    python::class_<RDKit::EnumerationStrategyBase,
                   RDKit::EnumerationStrategyBase *,
                   RDKit::EnumerationStrategyBase &, boost::noncopyable>(
        "EnumerationStrategyBase", python::no_init)
        .def("__nonzero__", &EnumerationStrategyBase__nonzero__)
        .def("Type", &EnumerationStrategyBase::type,
             "Returns the enumeration strategy type as a string.")
        .def("Skip", &EnumerationStrategyBase::skip,
             python::args("skipCount"),
             "Skip the next Nth results. note: this may be an expensive "
             "operation\n"
             "depending on the enumeration strategy used. It is recommended to "
             "use\n"
             "the enumerator state to advance to a known position")
        .def("Clone", python::pure_virtual(&EnumerationStrategyBase::Clone),
             python::return_value_policy<python::manage_new_object>())
        .def("GetNumPermutations", &EnumerationStrategyBase::getNumPermutations,
             "Returns the total number of results for this enumeration strategy.\n"
             "Note that some strategies are effectively infinite.")
        .def("CurrentPosition", &EnumerationStrategyBase::currentPosition,
             "Return the current indices into the arrays of reagents",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1> >())
        .def("Next", python::pure_virtual(&EnumerationStrategyBase::next),
             "Return the next indices into the arrays of reagents",
             python::return_internal_reference<
                 1, python::with_custodian_and_ward_postcall<0, 1> >())
        .def("Initialize", ToBBS);

    docString = "CartesianProductStrategy produces a standard walk through all possible\n"
        "reagent combinations:\n"
        "\n"
        "(0,0,0), (1,0,0), (2,0,0) ...\n";
    
    python::class_<RDKit::CartesianProductStrategy,
                   RDKit::CartesianProductStrategy*,
                   RDKit::CartesianProductStrategy&,
                   python::bases<EnumerationStrategyBase> >("CartesianProductStrategy",
                                                            docString.c_str(),
                                                            python::init<>())
      ;

    docString = "RandomSampleStrategy simply randomly samples from the reagent sets.\n"
        "Note that this strategy never halts and can produce duplicates.";
    python::class_<RDKit::RandomSampleStrategy,
                   RDKit::RandomSampleStrategy*,
                   RDKit::RandomSampleStrategy&,
                   python::bases<EnumerationStrategyBase> >("RandomSampleStrategy",
                                                            docString.c_str(),
                                                            python::init<>())
      ;

    docString = "RandomSampleAllBBsStrategy randomly samples from the reagent sets\n"
        "with the constraint that all building blocks are samples as early as possible.\n"
        "Note that this strategy never halts and can produce duplicates.";
    python::class_<RDKit::RandomSampleAllBBsStrategy,
                   RDKit::RandomSampleAllBBsStrategy*,
                   RDKit::RandomSampleAllBBsStrategy&,
                   python::bases<EnumerationStrategyBase> >("RandomSampleAllBBsStrategy",
                                                            docString.c_str(),
                                                            python::init<>())
      ;

    docString = "Randomly sample Pairs evenly from a collection of building blocks\n"
        "This is a good strategy for choosing a relatively small selection\n"
        "of building blocks from a larger set.  As the amount of work needed\n"
        "to retrieve the next evenly sample building block grows with the\n"
        "number of samples, this method performs progressively worse as the\n"
        "number of samples gets larger.\n"
        "See EnumeartionStrategyBase for more details.\n";
    
    python::class_<RDKit::EvenSamplePairsStrategy,
                   RDKit::EvenSamplePairsStrategy*,
                   RDKit::EvenSamplePairsStrategy&,
                   python::bases<EnumerationStrategyBase> >("EvenSamplePairsStrategy",
                                                            docString.c_str(),
                                                            python::init<>())
      ;
    
  }
};

}// end of namespace

void wrap_enumeration() {
  RDKit::enumeration_wrapper::wrap();
}

