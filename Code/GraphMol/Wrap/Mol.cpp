// $Id$
//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include "rdchem.h"
#include "seqs.hpp"
// ours
#include <RDBoost/pyint_api.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/python/iterator.hpp>
#include <boost/python/copy_non_const_reference.hpp>

namespace python = boost::python;


namespace RDKit {

  python::object MolToBinary(const ROMol &self){
    std::string res;
    MolPickler::pickleMol(self,res);
    python::object retval = python::object(python::handle<>(PyBytes_FromStringAndSize(res.c_str(),res.length())));
    return retval;
  }
  //
  // allows molecules to be pickled.
  //  since molecules have a constructor that takes a binary string
  //  we only need to provide getinitargs()
  //
  struct mol_pickle_suite : python::pickle_suite
  {
    static python::tuple
    getinitargs(const ROMol& self)
    {
      return python::make_tuple(MolToBinary(self));
    };
  };

  bool HasSubstructMatchStr(std::string pkl, const ROMol &query,
			    bool recursionPossible=true,bool useChirality=false,
                            bool useQueryQueryMatches=false){
    ROMol *mol;
    try {
      mol = new ROMol(pkl);
    } catch (...) {
      mol = NULL;
    }
    if(!mol){
      throw ValueErrorException("Null Molecule");
    }
    MatchVectType res;
    bool hasM=SubstructMatch(*mol,query,res,recursionPossible,useChirality,useQueryQueryMatches);
    delete mol;
    return hasM;
  }

  bool HasSubstructMatch(const ROMol &mol, const ROMol &query,
			 bool recursionPossible=true,bool useChirality=false,
                            bool useQueryQueryMatches=false){
    MatchVectType res;
    return SubstructMatch(mol,query,res,recursionPossible,useChirality,useQueryQueryMatches);
  }

  PyObject *convertMatches(MatchVectType &matches){
    PyObject *res = PyTuple_New(matches.size());
    MatchVectType::const_iterator i;
    for(i=matches.begin();i!=matches.end();i++){
      PyTuple_SetItem(res,i->first,PyInt_FromLong(i->second));
    }
    return res;
  }
  PyObject *GetSubstructMatch(const ROMol &mol, const ROMol &query,bool useChirality=false,
                            bool useQueryQueryMatches=false){
    MatchVectType matches;
    SubstructMatch(mol,query,matches,true,useChirality,useQueryQueryMatches);
    return convertMatches(matches);
  }

  PyObject *GetSubstructMatches(const ROMol &mol, const ROMol &query,bool uniquify=true,
                                bool useChirality=false,
                                bool useQueryQueryMatches=false){
    std::vector< MatchVectType >  matches;
    int matched = SubstructMatch(mol,query,matches,uniquify,true,useChirality,useQueryQueryMatches);
    PyObject *res = PyTuple_New(matched);
    for(int idx=0;idx<matched;idx++){
      PyTuple_SetItem(res,idx,convertMatches(matches[idx]));
    }
    return res;
  }

  unsigned int AddMolConformer(ROMol &mol, Conformer *conf, bool assignId=false) {
    Conformer *nconf = new Conformer(*conf);
    return mol.addConformer(nconf, assignId);
  }

  Conformer *GetMolConformer(ROMol &mol, int id=-1) {
    return &(mol.getConformer(id));
  }

  PyObject* GetMolConformers(ROMol &mol) {
    PyObject *res = PyTuple_New(mol.getNumConformers());
    ROMol::ConformerIterator ci;
    unsigned int i = 0;
    for (ci = mol.beginConformers(); ci != mol.endConformers(); ci++) {
      PyTuple_SetItem(res, i, python::converter::shared_ptr_to_python(*ci));
      i++;
    }
    return res;
  }
    
  std::string MolGetProp(const ROMol &mol,const char *key){
    if(!mol.hasProp(key)){
      PyErr_SetString(PyExc_KeyError,key);
      throw python::error_already_set();
    }
    std::string res;
    mol.getProp(key,res);
    return res;
  }

  int MolHasProp(const ROMol &mol,const char *key){
    int res = mol.hasProp(key);
    //std::cout << "key: "  << key << ": " << res << std::endl;
    return res;
  }
  void MolSetProp(const ROMol &mol,const char *key,std::string val,
		  bool computed=false){
    mol.setProp(key, val, computed);
  }

  void MolClearProp(const ROMol &mol,const char *key) {
    if(!mol.hasProp(key)){
      return;
    }
    mol.clearProp(key);
  }

  void MolClearComputedProps(const ROMol &mol) {
    mol.clearComputedProps();
  }

  void MolDebug(const ROMol &mol){
    mol.debugMol(std::cout);
  }

  // FIX: we should eventually figure out how to do iterators properly
  AtomIterSeq *MolGetAtoms(ROMol *mol){
    AtomIterSeq *res = new AtomIterSeq(mol->beginAtoms(),mol->endAtoms());
    return res;
  }
  QueryAtomIterSeq *MolGetAromaticAtoms(ROMol *mol){
    QueryAtom *qa=new QueryAtom();
    qa->setQuery(makeAtomAromaticQuery());
    QueryAtomIterSeq *res = new QueryAtomIterSeq(mol->beginQueryAtoms(qa),
                                                 mol->endQueryAtoms());
    return res;
  }
  QueryAtomIterSeq *MolGetQueryAtoms(ROMol *mol,QueryAtom *qa){
    QueryAtomIterSeq *res = new QueryAtomIterSeq(mol->beginQueryAtoms(qa),
                                                 mol->endQueryAtoms());
    return res;
  }

  //AtomIterSeq *MolGetHeteros(ROMol *mol){
  //  AtomIterSeq *res = new AtomIterSeq(mol->beginHeteros(),
  //                                     mol->endHeteros());
  //  return res;
  //}
  BondIterSeq *MolGetBonds(ROMol *mol){
    BondIterSeq *res = new BondIterSeq(mol->beginBonds(),mol->endBonds());
    return res;
  }

  int getMolNumAtoms(const ROMol &mol, int onlyHeavy, bool onlyExplicit){
    if(onlyHeavy>-1){
      BOOST_LOG(rdWarningLog)<<"WARNING: the onlyHeavy argument to mol.GetNumAtoms() has been deprecated. Please use the onlyExplicit argument instead or mol.GetNumHeavyAtoms() if you want the heavy atom count."<<std::endl;
      return mol.getNumAtoms(onlyHeavy);
    }
    return mol.getNumAtoms(onlyExplicit);
  }



  std::string molClassDoc = "The Molecule class.\n\n\
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
struct mol_wrapper {
  static void wrap(){
    python::register_exception_translator<ConformerException>(&rdExceptionTranslator);

    python::class_<ROMol,ROMOL_SPTR,boost::noncopyable>("Mol",
			  molClassDoc.c_str(),
			  python::init<>("Constructor, takes no arguments"))
      .def(python::init<const std::string &>())
      .def(python::init<const ROMol &>())
      .def("GetNumAtoms",getMolNumAtoms,
	   (python::arg("onlyHeavy")=-1,
            python::arg("onlyExplicit")=true),
	   "Returns the number of atoms in the molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - onlyExplicit: (optional) include only explicit atoms (atoms in the molecular graph)\n"
	   "                    defaults to 1.\n"
	   "  NOTE: the onlyHeavy argument is deprecated\n"

)
      .def("GetNumHeavyAtoms",&ROMol::getNumHeavyAtoms,
	   "Returns the number of heavy atoms (atomic number >1) in the molecule.\n\n"
           )
      .def("GetAtomWithIdx",(Atom * (ROMol::*)(unsigned int))&ROMol::getAtomWithIdx,
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a particular Atom.\n\n"
	   "  ARGUMENTS:\n"
	   "    - idx: which Atom to return\n\n"
	   "  NOTE: atom indices start at 0\n")

      .def("GetNumBonds",&ROMol::getNumBonds,
	   (python::arg("onlyHeavy")=true),
	   "Returns the number of Bonds in the molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)\n"
	   "                  defaults to 1.\n")

      .def("GetBondWithIdx",(Bond * (ROMol::*)(unsigned int))&ROMol::getBondWithIdx,
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a particular Bond.\n\n"
	   "  ARGUMENTS:\n"
	   "    - idx: which Bond to return\n\n"
	   "  NOTE: bond indices start at 0\n")

      .def("GetNumConformers", &ROMol::getNumConformers,
           "Return the number of conformations on the molecule")

      .def("AddConformer", AddMolConformer,
	   (python::arg("self"),python::arg("conf"),
	    python::arg("assignId")=false),
	   "Add a conformer to the molecule and return the conformer ID")

      .def("GetConformer", GetMolConformer,
	   (python::arg("self"),python::arg("id")=-1),
	   "Get the conformer with a specified ID",
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >())

      .def("GetConformers", GetMolConformers,
           "Get all the conformers as a tuple")

      .def("RemoveAllConformers", &ROMol::clearConformers,
           "Remove all the conformations on the molecule")

      .def("RemoveConformer", &ROMol::removeConformer,
           "Remove the conformer with the specified ID")
      .def("GetBondBetweenAtoms", 
	   (Bond *(ROMol::*)(unsigned int,unsigned int))&ROMol::getBondBetweenAtoms,
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns the bond between two atoms, if there is one.\n\n"
	   "  ARGUMENTS:\n"
	   "    - idx1,idx2: the Atom indices\n\n"
	   "  Returns:\n"
	   "    The Bond between the two atoms, if such a bond exists.\n"
	   "    If there is no Bond between the atoms, None is returned instead.\n\n"
	   "  NOTE: bond indices start at 0\n"
	   )

      // substructures
      .def("HasSubstructMatch",HasSubstructMatch,
	   (python::arg("self"),python::arg("query"),
	    python::arg("recursionPossible")=true,
	    python::arg("useChirality")=false,
            python::arg("useQueryQueryMatches")=false),
	   "Queries whether or not the molecule contains a particular substructure.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule\n\n"
	   "    - recursionPossible: (optional)\n\n"
	   "    - useChirality: enables the use of stereochemistry in the matching\n\n"
	   "    - useQueryQueryMatches: use query-query matching logic\n\n"
	   "  RETURNS: True or False\n")
      .def("GetSubstructMatch",GetSubstructMatch,
	   (python::arg("self"),python::arg("query"),
	    python::arg("useChirality")=false,
            python::arg("useQueryQueryMatches")=false),
	   "Returns the indices of the molecule's atoms that match a substructure query.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule\n\n"
	   "    - useChirality: enables the use of stereochemistry in the matching\n\n"
	   "    - useQueryQueryMatches: use query-query matching logic\n\n"
	   "  RETURNS: a tuple of integers\n\n"
	   "  NOTES:\n"
	   "     - only a single match is returned\n"
	   "     - the ordering of the indices corresponds to the atom ordering\n"
	   "         in the query. For example, the first index is for the atom in\n"
	   "         this molecule that matches the first atom in the query.\n" 
	   )

      .def("GetSubstructMatches",
	   GetSubstructMatches,
	   (python::arg("self"),python::arg("query"),
	    python::arg("uniquify")=true,
	    python::arg("useChirality")=false,
            python::arg("useQueryQueryMatches")=false),
	   "Returns tuples of the indices of the molecule's atoms that match a substructure query.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule.\n"
	   "    - uniquify: (optional) determines whether or not the matches are uniquified.\n"
	   "                Defaults to 1.\n\n"
	   "    - useChirality: enables the use of stereochemistry in the matching\n\n"
	   "    - useQueryQueryMatches: use query-query matching logic\n\n"
	   "  RETURNS: a tuple of tuples of integers\n\n"
	   "  NOTE:\n"
	   "     - the ordering of the indices corresponds to the atom ordering\n"
	   "         in the query. For example, the first index is for the atom in\n"
	   "         this molecule that matches the first atom in the query.\n")

      // properties
      .def("SetProp",MolSetProp,
	   (python::arg("self"), python::arg("key"),
	    python::arg("val"), python::arg("computed")=false),
	   "Sets a molecular property\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to be set (a string).\n"
	   "    - value: the property value (a string).\n"
	   "    - computed: (optional) marks the property as being computed.\n"
	   "                Defaults to 0.\n\n")
      .def("HasProp",MolHasProp,
	   "Queries a molecule to see if a particular property has been assigned.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to check for (a string).\n")

      .def("GetProp",MolGetProp,
	   "Returns the value of the property.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to return (a string).\n\n"
	   "  RETURNS: a string\n\n"
	   "  NOTE:\n"
	   "    - If the property has not been set, a KeyError exception will be raised.\n")

      .def("ClearProp", MolClearProp,
	   "Removes a property from the molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to clear (a string).\n")

      .def("ClearComputedProps", MolClearComputedProps,
	   "Removes all computed properties from the molecule.\n\n")

      .def("UpdatePropertyCache", &ROMol::updatePropertyCache,
	   (python::arg("self"),python::arg("strict")=true),
	    "Regenerates computed properties like implicit valence and ring information.\n\n")


      .def("GetPropNames",&ROMol::getPropList,
	   (python::arg("self"),python::arg("includePrivate")=false,
	    python::arg("includeComputed")=false),
	   "Returns a tuple with all property names for this molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - includePrivate: (optional) toggles inclusion of private properties in the result set.\n"
	   "                      Defaults to 0.\n"
	   "    - includeComputed: (optional) toggles inclusion of computed properties in the result set.\n"
	   "                      Defaults to 0.\n\n"
	   "  RETURNS: a tuple of strings\n")

      .def("GetAtoms",MolGetAtoms,
           python::return_value_policy<python::manage_new_object,
           python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a read-only sequence containing all of the molecule's Atoms.\n")
      .def("GetAromaticAtoms",MolGetAromaticAtoms,
           python::return_value_policy<python::manage_new_object,
           python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a read-only sequence containing all of the molecule's aromatic Atoms.\n")
      .def("GetAtomsMatchingQuery",MolGetQueryAtoms,
           python::return_value_policy<python::manage_new_object,
           python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a read-only sequence containing all of the atoms in a molecule that match the query atom.\n")

      .def("GetBonds",MolGetBonds,
           python::return_value_policy<python::manage_new_object,
           python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns a read-only sequence containing all of the molecule's Bonds.\n")

      // enable pickle support
      .def_pickle(mol_pickle_suite())

      .def("Debug",MolDebug,
	   "Prints debugging information about the molecule.\n")

      .def("ToBinary",MolToBinary,
	   "Returns a binary string representation of the molecule.\n")

      .def("GetRingInfo",&ROMol::getRingInfo,
      python::return_value_policy<python::reference_existing_object>(),
       "Returns the number of molecule's RingInfo object.\n\n")
      ;
        
    // ---------------------------------------------------------------------------------------------
    python::def("_HasSubstructMatchStr",
                HasSubstructMatchStr,
                (python::arg("pkl"),python::arg("query"),
		 python::arg("recursionPossible")=true,
		 python::arg("useChirality")=false,
		 python::arg("useQueryQueryMatches")=false
                 ),
		"This function is included to speed substructure queries from databases, \n"
		"it's probably not of\n"
		"general interest.\n\n"
		"  ARGUMENTS:\n"
		"    - pkl: a Molecule pickle\n\n"
		"    - query: a Molecule\n\n"
		"    - recursionPossible: (optional)\n\n"
		"    - useChirality: (optional)\n\n"
                "    - useQueryQueryMatches: use query-query matching logic\n\n"
		"  RETURNS: True or False\n");

    
  };
};
}// end of namespace
void wrap_mol() {
  RDKit::mol_wrapper::wrap();
}
