// $Id$
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>

namespace python = boost::python;


void rdChemicalReactionParserExceptionTranslator(RDKit::ChemicalReactionParserException const& x){
  std::ostringstream ss;
  ss << "ChemicalReactionParserException: " << x.message();
  PyErr_SetString(PyExc_ValueError,ss.str().c_str());
}
void rdChemicalReactionExceptionTranslator(RDKit::ChemicalReactionException const& x){
  std::ostringstream ss;
  ss << "ChemicalParserException: " << x.message();
  PyErr_SetString(PyExc_ValueError,ss.str().c_str());
}

namespace RDKit {
  template <typename T>
  PyObject* RunReactants(ChemicalReaction *self,T reactants){
    if(!self->isInitialized()){
      self->initReactantMatchers();
    }
    MOL_SPTR_VECT reacts;
    unsigned int len1 = python::extract<unsigned int>(reactants.attr("__len__")());
    reacts.resize(len1);
    for(unsigned int i=0;i<len1;++i){
      reacts[i] = python::extract<ROMOL_SPTR>(reactants[i]);
    }
    std::vector<MOL_SPTR_VECT> mols;
    mols = self->runReactants(reacts);
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

  python::tuple ValidateReaction(const ChemicalReaction *self,bool silent=false){
    unsigned int numWarn,numError;
    self->validate(numWarn,numError,silent);
    return python::make_tuple(numWarn,numError);
  }

}

BOOST_PYTHON_MODULE(rdChemReactions) {
  python::scope().attr("__doc__") =
    "Module containing classes and functions for working with chemical reactions."
    ;

  python::register_exception_translator<RDKit::ChemicalReactionParserException>(&rdChemicalReactionParserExceptionTranslator);
  python::register_exception_translator<RDKit::ChemicalReactionException>(&rdChemicalReactionExceptionTranslator);
    
  std::string docString = "A class for storing and applying chemical reactions.\n\
\n\
Sample Usage:\n\
>>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')\n\
>>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))\n\
>>> products = rxn.RunReactants(reacts)\n\
>>> len(products)\n\
1\n\
>>> len(products[0])\n\
1\n\
>>> Chem.MolToSmiles(products[0])\n\
'CN(C)C=O'\n\
\n\
";
  python::class_<RDKit::ChemicalReaction>("ChemicalReaction",docString.c_str())
    .def("GetNumReactantTemplates",&RDKit::ChemicalReaction::getNumReactantTemplates,
         "returns the number of reactants this reaction expects")
    .def("GetNumProductTemplates",&RDKit::ChemicalReaction::getNumProductTemplates,
         "returns the number of products this reaction generates")
    .def("AddReactantTemplate",&RDKit::ChemicalReaction::addReactantTemplate,
         "adds a reactant (a Molecule) to the reaction")
    .def("AddProductTemplate",&RDKit::ChemicalReaction::addProductTemplate,
         "adds a product (a Molecule)")
    .def("RunReactants",(PyObject *(*)(RDKit::ChemicalReaction *,python::tuple))RDKit::RunReactants,
         "apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples")
    .def("RunReactants",(PyObject *(*)(RDKit::ChemicalReaction *,python::list))RDKit::RunReactants,
         "apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples")
    .def("Initialize",&RDKit::ChemicalReaction::initReactantMatchers,
         "initializes the reaction so that it can be used")
    .def("IsInitialized",&RDKit::ChemicalReaction::isInitialized,
         "checks if the reaction is ready for use")
    .def("Validate",&RDKit::ValidateReaction,
         (python::arg("self"),python::arg("silent")=false),
         "checks the reaction for potential problems, returns (numWarnings,numErrors)")
  ;

  def("ReactionFromSmarts",RDKit::RxnSmartsToChemicalReaction,
      "construct a ChemicalReaction from a reaction SMARTS string",
      python::return_value_policy<python::manage_new_object>());
  def("ReactionFromRxnFile",RDKit::RxnFileToChemicalReaction,
      "construct a ChemicalReaction from an MDL rxn file",
      python::return_value_policy<python::manage_new_object>());
  def("ReactionFromRxnBlock",RDKit::RxnBlockToChemicalReaction,
      "construct a ChemicalReaction from an string in MDL rxn format",
      python::return_value_policy<python::manage_new_object>());
  
 

              
}
              
