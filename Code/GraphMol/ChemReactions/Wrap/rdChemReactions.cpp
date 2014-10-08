// $Id$
//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/Depictor/DepictUtils.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>

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
  python::object ReactionToBinary(const ChemicalReaction &self){
    std::string res;
    ReactionPickler::pickleReaction(self,res);
    python::object retval = python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length())));
    return retval;
  }
  //
  // allows reactions to be pickled.
  //
  struct reaction_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const ChemicalReaction & self)
    {
      return python::make_tuple(ReactionToBinary(self));
    };
  };


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
      if(!reacts[i]) throw_value_error("reaction called with None reactants");
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

  ROMol * GetProductTemplate(const ChemicalReaction *self,unsigned int which){
    if(which>=self->getNumProductTemplates()){
      throw_value_error("requested template index too high");
    }
    MOL_SPTR_VECT::const_iterator iter=self->beginProductTemplates();
    iter += which;
    ROMol *res = const_cast<ROMol *>(iter->get());
    return res;
  }
  ROMol * GetReactantTemplate(const ChemicalReaction *self,unsigned int which){
    if(which>=self->getNumReactantTemplates()){
      throw_value_error("requested template index too high");
    }
    MOL_SPTR_VECT::const_iterator iter=self->beginReactantTemplates();
    iter += which;
    ROMol *res = const_cast<ROMol *>(iter->get());
    return res;
  }
  ROMol * GetAgentTemplate(const ChemicalReaction *self,unsigned int which){
    if(which>=self->getNumAgentTemplates()){
      throw_value_error("requested template index too high");
    }
    MOL_SPTR_VECT::const_iterator iter=self->beginAgentTemplates();
    iter += which;
    ROMol *res = const_cast<ROMol *>(iter->get());
    return res;
  }

  void RemoveUnmappedReactantTemplates(ChemicalReaction *self, double thresholdUnmappedAtoms, bool moveToAgentTemplates,
		  python::object targetList){
    if(targetList==python::object()){
    	self->removeUnmappedReactantTemplates(thresholdUnmappedAtoms, moveToAgentTemplates);
	}
    else{
	  MOL_SPTR_VECT tmp;
	  self->removeUnmappedReactantTemplates(thresholdUnmappedAtoms, moveToAgentTemplates, &tmp);
	  python::list molList = python::extract<python::list>(targetList);
	  if(tmp.size() > 0){
	    for(unsigned int i = 0; i < tmp.size(); i++){
		  molList.append(tmp.at(i));
	    }
	  }
    }
  }

  void RemoveUnmappedProductTemplates(ChemicalReaction *self, double thresholdUnmappedAtoms, bool moveToAgentTemplates,
		  python::object targetList){
    if(targetList==python::object()){
    	self->removeUnmappedProductTemplates(thresholdUnmappedAtoms, moveToAgentTemplates);
	}
    else{
	  MOL_SPTR_VECT tmp;
	  self->removeUnmappedProductTemplates(thresholdUnmappedAtoms, moveToAgentTemplates, &tmp);
	  python::list molList = python::extract<python::list>(targetList);
	  if(tmp.size() > 0){
	    for(unsigned int i = 0; i < tmp.size(); i++){
		  molList.append(tmp.at(i));
	    }
	  }
    }
  }

  void Compute2DCoordsForReaction(RDKit::ChemicalReaction &rxn,
                                  double spacing=2.0,
                                  bool updateProps=true,
                                  bool canonOrient=false,
                                  unsigned int nFlipsPerSample=0,
                                  unsigned int nSamples=0,
                                  int sampleSeed=0,
                                  bool permuteDeg4Nodes=false,
                                  double bondLength=-1){
    double oBondLen=RDDepict::BOND_LEN;
    if(bondLength>0){
      RDDepict::BOND_LEN=bondLength;
    }
    RDDepict::compute2DCoordsForReaction(rxn,spacing,updateProps,canonOrient,
                                         nFlipsPerSample,nSamples,sampleSeed,
                                         permuteDeg4Nodes);
    if(bondLength>0){
      RDDepict::BOND_LEN=oBondLen;
    }
  }

  bool IsMoleculeReactantOfReaction(const ChemicalReaction &rxn,const ROMol &mol){
    unsigned int which;
    return isMoleculeReactantOfReaction(rxn,mol,which);
  }
  bool IsMoleculeProductOfReaction(const ChemicalReaction &rxn,const ROMol &mol){
    unsigned int which;
    return isMoleculeProductOfReaction(rxn,mol,which);
  }
  bool IsMoleculeAgentOfReaction(const ChemicalReaction &rxn,const ROMol &mol){
    unsigned int which;
    return isMoleculeAgentOfReaction(rxn,mol,which);
  }
  

  ChemicalReaction *ReactionFromSmarts(const char *smarts,
                                       python::dict replDict,
                                       bool useSmiles){
    PRECONDITION(smarts,"null SMARTS string");
    std::map<std::string,std::string> replacements;
    for(unsigned int i=0;i<python::extract<unsigned int>(replDict.keys().attr("__len__")());++i){
      replacements[python::extract<std::string>(replDict.keys()[i])]=python::extract<std::string>(replDict.values()[i]);
    }
    ChemicalReaction *res; 
    res = RxnSmartsToChemicalReaction(smarts,&replacements,useSmiles);
    return res;
  }
  
  python::object GetReactingAtoms(const ChemicalReaction &self,bool mappedAtomsOnly){
    python::list res;
    VECT_INT_VECT rAs=getReactingAtoms(self,mappedAtomsOnly);
    for(VECT_INT_VECT_I rIt=rAs.begin();rIt!=rAs.end();++rIt){
      res.append(python::tuple(*rIt));
    }
    return python::tuple(res);
  }

  python::object AddRecursiveQueriesToReaction(ChemicalReaction &self, python::dict queryDict,
                                     std::string propName, bool getLabels=false){
    // transform dictionary into map
    std::map<std::string, ROMOL_SPTR> queries;
    for(unsigned int i=0;i<python::extract<unsigned int>(queryDict.keys().attr("__len__")());++i){
      ROMol *m = python::extract<ROMol *>(queryDict.values()[i]);
      ROMOL_SPTR nm(new ROMol(*m));
      std::string k = python::extract<std::string>(queryDict.keys()[i]);
      queries[k]=nm;
    }

    if (getLabels) {
      std::vector<std::vector<std::pair<unsigned int, std::string> > > labels;
      addRecursiveQueriesToReaction(self, queries, propName, &labels);

      // transform labels into python::tuple(python::tuple(python::tuple))
      python::list reactantLabels;
      for (unsigned int i=0; i<labels.size(); ++i) {
        python::list tmpLabels;
        for (unsigned int j=0; j<labels[i].size(); ++j) {
          python::list tmpPair;
          tmpPair.append(labels[i][j].first);
          tmpPair.append(labels[i][j].second);
          tmpLabels.append(python::tuple(tmpPair));
        }
        reactantLabels.append(python::tuple(tmpLabels));
      }
      return python::tuple(reactantLabels);
    } else {
      addRecursiveQueriesToReaction(self, queries, propName);
      return python::object(); // this is None
    }
  }
}

BOOST_PYTHON_MODULE(rdChemReactions) {
  python::scope().attr("__doc__") =
    "Module containing classes and functions for working with chemical reactions."
    ;

  python::register_exception_translator<RDKit::ChemicalReactionParserException>(&rdChemicalReactionParserExceptionTranslator);
  python::register_exception_translator<RDKit::ChemicalReactionException>(&rdChemicalReactionExceptionTranslator);

  python::enum_<RDKit::FingerprintType>("FingerprintType")
          .value("AtomPairFP", RDKit::AtomPairFP)
          .value("TopologicalTorsion", RDKit::TopologicalTorsion)
          .value("MorganFP", RDKit::MorganFP)
          .value("RDKitFP", RDKit::RDKitFP)
          .value("PatternFP", RDKit::PatternFP)
          ;
  std::string docStringReactionFPParams =
		  "A class for storing parameters to manipulate the calculation of fingerprints of chemical reactions.";

  python::class_<RDKit::ReactionFingerprintParams>("ReactionFingerprintParams",docStringReactionFPParams.c_str(),
		  python::init<>("Constructor, takes no arguments"))
   .def(python::init<bool, double, unsigned int, int, unsigned int, RDKit::FingerprintType>())
   .def_readwrite("fpSize", &RDKit::ReactionFingerprintParams::fpSize)
   .def_readwrite("fpType", &RDKit::ReactionFingerprintParams::fpType)
   .def_readwrite("bitRatioAgents", &RDKit::ReactionFingerprintParams::bitRatioAgents)
   .def_readwrite("nonAgentWeight", &RDKit::ReactionFingerprintParams::nonAgentWeight)
   .def_readwrite("agentWeight", &RDKit::ReactionFingerprintParams::agentWeight)
   .def_readwrite("includeAgents", &RDKit::ReactionFingerprintParams::includeAgents)
  ;

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
  python::class_<RDKit::ChemicalReaction>("ChemicalReaction",docString.c_str(),
                                          python::init<>("Constructor, takes no arguments"))
    .def(python::init<const std::string &>())
    .def("GetNumReactantTemplates",&RDKit::ChemicalReaction::getNumReactantTemplates,
         "returns the number of reactants this reaction expects")
    .def("GetNumProductTemplates",&RDKit::ChemicalReaction::getNumProductTemplates,
         "returns the number of products this reaction generates")
    .def("GetNumAgentTemplates",&RDKit::ChemicalReaction::getNumAgentTemplates,
         "returns the number of agents this reaction expects")
    .def("AddReactantTemplate",&RDKit::ChemicalReaction::addReactantTemplate,
         "adds a reactant (a Molecule) to the reaction")
    .def("AddProductTemplate",&RDKit::ChemicalReaction::addProductTemplate,
         "adds a product (a Molecule)")
    .def("AddAgentTemplate",&RDKit::ChemicalReaction::addAgentTemplate,
         "adds a agent (a Molecule)")
    .def("RemoveUnmappedReactantTemplates",RDKit::RemoveUnmappedReactantTemplates,
        (python::arg("self"), python::arg("thresholdUnmappedAtoms")=0.2,python::arg("moveToAgentTemplates")=true,
         python::arg("targetList")=python::object()),
        "Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList")
    .def("RemoveUnmappedProductTemplates",RDKit::RemoveUnmappedProductTemplates,
        (python::arg("self"), python::arg("thresholdUnmappedAtoms")=0.2,python::arg("moveToAgentTemplates")=true,
         python::arg("targetList")=python::object()),
         "Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList")
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
    .def("GetProductTemplate",&RDKit::GetProductTemplate,
         (python::arg("self"),python::arg("which")),
         python::return_value_policy<python::reference_existing_object>(),
         "returns one of our product templates")
    .def("GetReactantTemplate",&RDKit::GetReactantTemplate,
         (python::arg("self"),python::arg("which")),
         python::return_value_policy<python::reference_existing_object>(),
         "returns one of our reactant templates")
    .def("GetAgentTemplate",&RDKit::GetAgentTemplate,
         (python::arg("self"),python::arg("which")),
         python::return_value_policy<python::reference_existing_object>(),
         "returns one of our agent templates")
    .def("_setImplicitPropertiesFlag",&RDKit::ChemicalReaction::setImplicitPropertiesFlag,
         (python::arg("self"),python::arg("val")),
         "EXPERT USER: indicates that the reaction can have implicit properties")
    .def("_getImplicitPropertiesFlag",&RDKit::ChemicalReaction::getImplicitPropertiesFlag,
         (python::arg("self")),
         "EXPERT USER: returns whether or not the reaction can have implicit properties")
    .def("ToBinary",RDKit::ReactionToBinary,
         "Returns a binary string representation of the reaction.")
    .def("IsMoleculeReactant",RDKit::IsMoleculeReactantOfReaction,
         "returns whether or not the molecule has a substructure match to one of the reactants.")
    .def("IsMoleculeProduct",RDKit::IsMoleculeProductOfReaction,
         "returns whether or not the molecule has a substructure match to one of the products.")
    .def("IsMoleculeAgent",RDKit::IsMoleculeAgentOfReaction,
         "returns whether or not the molecule has a substructure match to one of the agents.")
    .def("GetReactingAtoms",&RDKit::GetReactingAtoms,
         (python::arg("self"),python::arg("mappedAtomsOnly")=false),
         "returns a sequence of sequences with the atoms that change in the reaction")
    .def("AddRecursiveQueriesToReaction", RDKit::AddRecursiveQueriesToReaction,
         (python::arg("reaction"), python::arg("queries")=python::dict(),
          python::arg("propName")="molFileValue", python::arg("getLabels")=false),
         "adds recursive queries and returns reactant labels")
    // enable pickle support
    .def_pickle(RDKit::reaction_pickle_suite())
  ;


  python::def("ReactionFromSmarts",RDKit::ReactionFromSmarts,
              (python::arg("SMARTS"),
               python::arg("replacements")=python::dict(),
               python::arg("useSmiles")=false),
              "construct a ChemicalReaction from a reaction SMARTS string. \n\
see the documentation for rdkit.Chem.MolFromSmiles for an explanation\n\
of the replacements argument.",
      python::return_value_policy<python::manage_new_object>());
  python::def("ReactionFromRxnFile",RDKit::RxnFileToChemicalReaction,
      "construct a ChemicalReaction from an MDL rxn file",
      python::return_value_policy<python::manage_new_object>());
  python::def("ReactionFromRxnBlock",RDKit::RxnBlockToChemicalReaction,
      "construct a ChemicalReaction from an string in MDL rxn format",
      python::return_value_policy<python::manage_new_object>());
  python::def("ReactionFromMolecule",RDKit::RxnMolToChemicalReaction,
      "construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set",
      python::return_value_policy<python::manage_new_object>());

  python::def("ReactionToSmarts",RDKit::ChemicalReactionToRxnSmarts,
      (python::arg("reaction")),
      "construct a reaction SMARTS string for a ChemicalReaction");
  python::def("ReactionToSmiles",RDKit::ChemicalReactionToRxnSmiles,
      (python::arg("reaction")),
      "construct a reaction SMILES string for a ChemicalReaction");
  python::def("ReactionToRxnBlock",RDKit::ChemicalReactionToRxnBlock,
      (python::arg("reaction"), python::arg("separateAgents")=false),
      "construct a string in MDL rxn format for a ChemicalReaction");
  python::def("ReactionToMolecule",RDKit::ChemicalReactionToRxnMol,
      (python::arg("reaction")),
      "construct a molecule for a ChemicalReaction with RXN role property set",
      python::return_value_policy<python::manage_new_object>());

  docString = "Compute 2D coordinates for a reaction. \n\
  ARGUMENTS: \n\n\
     reaction - the reaction of interest\n\
     spacing - the amount of space left between components of the reaction\n\
     canonOrient - orient the reactants and products in a canonical way\n\
     updateProps - if set, properties such as conjugation and\n\
        hybridization will be calculated for the reactant and product\n\
        templates before generating coordinates. This should result in\n\
        better depictions, but can lead to errors in some cases.\n\
     nFlipsPerSample - number of rotatable bonds that are\n\
                flipped at random at a time.\n\
     nSample - Number of random samplings of rotatable bonds.\n\
     sampleSeed - seed for the random sampling process.\n\
     permuteDeg4Nodes - allow permutation of bonds at a degree 4\n\
                 node during the sampling process \n\
     bondLength - change the default bond length for depiction \n\n";
  python::def("Compute2DCoordsForReaction",
              RDKit::Compute2DCoordsForReaction,
	      (python::arg("reaction"),
	       python::arg("spacing")=2.0,
	       python::arg("updateProps")=true,
	       python::arg("canonOrient")=true,
               python::arg("nFlipsPerSample")=0,
               python::arg("nSample")=0,
               python::arg("sampleSeed")=0,
               python::arg("permuteDeg4Nodes")=false,
               python::arg("bondLength")=-1.0),
	      docString.c_str());

  python::def("CreateDifferenceFingerprintForReaction",RDKit::DifferenceFingerprintChemReaction,
              (python::arg("reaction"),
               python::arg("ReactionFingerPrintParams")=RDKit::DefaultDifferenceFPParams),
              "construct a difference fingerprint for a ChemicalReaction by subtracting the reactant "
              "fingerprint from the product fingerprint",
      python::return_value_policy<python::manage_new_object>());

  python::def("CreateStructuralFingerprintForReaction",RDKit::StructuralFingerprintChemReaction,
              (python::arg("reaction"),
               python::arg("ReactionFingerPrintParams")=RDKit::DefaultStructuralFPParams),
              "construct a structural fingerprint for a ChemicalReaction by concatenating the reactant "
              "fingerprint and the product fingerprint",
      python::return_value_policy<python::manage_new_object>());

  python::def("IsReactionTemplateMoleculeAgent",RDKit::isReactionTemplateMoleculeAgent,
		  (python::arg("molecule"), python::arg("agentThreshold")),
          "tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold");
  python::def("HasReactionAtomMapping",RDKit::hasReactionAtomMapping,
          "tests if a reaction obtains any atom mapping");
  python::def("HasReactionSubstructMatch",RDKit::hasReactionSubstructMatch,
		  (python::arg("reaction"), python::arg("queryReaction"), python::arg("includeAgents")=false),
          "tests if the queryReaction is a substructure of a reaction");
  python::def("HasAgentTemplateSubstructMatch",RDKit::hasAgentTemplateSubstructMatch,
		  (python::arg("reaction"), python::arg("queryReaction")),
          "tests if the agents of a queryReaction are the same as those of a reaction");
  python::def("HasProductTemplateSubstructMatch",RDKit::hasProductTemplateSubstructMatch,
		  (python::arg("reaction"), python::arg("queryReaction")),
          "tests if the products of a queryReaction are substructures of the products of a reaction");
  python::def("HasReactantTemplateSubstructMatch",RDKit::hasReactantTemplateSubstructMatch,
		  (python::arg("reaction"), python::arg("queryReaction")),
          "tests if the reactants of a queryReaction are substructures of the reactants of a reaction");

}
              
