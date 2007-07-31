// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

#include "rdchem.h"
#include "seqs.hpp"
// ours
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/python/iterator.hpp>
#include <boost/python/copy_non_const_reference.hpp>

namespace python = boost::python;


namespace RDKit {

  std::string MolToBinary(const ROMol &self){
    std::string res;
    MolPickler::pickleMol(self,res);
    return res;
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
			    bool recursionPossible=true,bool useChirality=false){
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
    bool hasM=SubstructMatch(*mol,query,res,recursionPossible,useChirality);
    delete mol;
    return hasM;
  }

  bool HasSubstructMatch(const ROMol &mol, const ROMol &query,
			 bool recursionPossible=true,bool useChirality=false){
    MatchVectType res;
    return SubstructMatch(mol,query,res,recursionPossible,useChirality);
  }

  PyObject *convertMatches(MatchVectType &matches){
    PyObject *res = PyTuple_New(matches.size());
    MatchVectType::const_iterator i;
    for(i=matches.begin();i!=matches.end();i++){
      PyTuple_SetItem(res,i->first,PyInt_FromLong(i->second));
    }
    return res;
  }
  PyObject *GetSubstructMatch(const ROMol &mol, const ROMol &query,bool useChirality=false){
    MatchVectType matches;
    SubstructMatch(mol,query,matches,true,useChirality);
    return convertMatches(matches);
  }

  PyObject *GetSubstructMatches(const ROMol &mol, const ROMol &query,bool uniquify=true,bool useChirality=false){
    std::vector< MatchVectType >  matches;
    int matched = SubstructMatch(mol,query,matches,uniquify,true,useChirality);
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

#if 1
  std::vector<std::string> MolGetPropNames(const ROMol &mol,
					   bool includePrivate=false,
					   bool includeComputed=false){
    std::vector<std::string> res,computed,tmp;
    tmp=mol.getPropList();
    if(!includeComputed) mol.getProp("computedProps",computed);
    // we'll never return this:
    computed.push_back("computedProps");
    std::vector<std::string>::iterator pos = tmp.begin();
    while(pos!=tmp.end()){
      if((includePrivate || (*pos)[0]!='_') &&
	 std::find(computed.begin(),computed.end(),*pos)==computed.end()){
	std::string holder;
	try {
	  mol.getProp(*pos,holder);
	  res.push_back(*pos);
	} catch (const boost::bad_any_cast &) {
	  ;
	}
	
      }
      pos++;
    }
    return res;
  }
#else
  python::list MolGetPropNames(const ROMol &mol,
					   bool includePrivate=false,
					   bool includeComputed=false){
    std::vector<std::string> computed,tmp;
    python::list res;
    tmp=mol.getPropList();
    if(!includeComputed) mol.getProp("computedProps",computed);
    // we'll never return this:
    computed.push_back("computedProps");
    std::vector<std::string>::iterator pos = tmp.begin();
    while(pos!=tmp.end()){
      if((includePrivate || (*pos)[0]!='_') &&
	 std::find(computed.begin(),computed.end(),*pos)==computed.end()){
	std::string holder;
	try {
	  mol.getProp(*pos,holder);
	  res.append(*pos);
	} catch (const boost::bad_any_cast &) {
	  ;
	}
	
      }
      pos++;
    }
    return res;
  }
#endif
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
    mol.clearProp(key);
  }

  void MolClearComputedProps(const ROMol &mol) {
    mol.clearComputedProps();
  }

  void MolDebug(const ROMol &mol){
    mol.debugMol(std::cout);
  }

#if 0
  // FIX: we should eventually figure out how to do iterators properly
  //  so that these tuples don't have to be built
  PyObject *MolGetAtoms(ROMol *mol){
    python::list res;
    for(ROMol::AtomIterator i=mol->beginAtoms();i!=mol->endAtoms();i++){
      res.append(*i);
    }
    //return python::incref(python::tuple(res).ptr());
    return python::incref(res.ptr());
  }
#else
  // FIX: we should eventually figure out how to do iterators properly
  AtomIterSeq *MolGetAtoms(ROMol *mol){
    AtomIterSeq *res = new AtomIterSeq(mol->beginAtoms(),mol->endAtoms());
    return res;
  }
  AromaticAtomIterSeq *MolGetAromaticAtoms(ROMol *mol){
    AromaticAtomIterSeq *res = new AromaticAtomIterSeq(mol->beginAromaticAtoms(),
						       mol->endAromaticAtoms());
    return res;
  }
  HeteroatomIterSeq *MolGetHeteros(ROMol *mol){
    HeteroatomIterSeq *res = new HeteroatomIterSeq(mol->beginHeteros(),
						   mol->endHeteros());
    return res;
  }
  BondIterSeq *MolGetBonds(ROMol *mol){
    BondIterSeq *res = new BondIterSeq(mol->beginBonds(),mol->endBonds());
    return res;
  }
#endif

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
      .def("GetNumAtoms",&ROMol::getNumAtoms,
	   (python::arg("onlyHeavy")=true),
	   "Returns the number of Atoms in the molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - onlyHeavy: (optional) include only heavy atoms (not Hs)\n"
	   "                 defaults to 1.\n")
      .def("GetAtomWithIdx",(ROMol::GRAPH_NODE_TYPE (ROMol::*)(unsigned int))&ROMol::getAtomWithIdx,
	   python::return_value_policy<python::reference_existing_object>(),
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

      .def("GetBondWithIdx",(ROMol::GRAPH_EDGE_TYPE (ROMol::*)(unsigned int))&ROMol::getBondWithIdx,
	   python::return_value_policy<python::reference_existing_object>(),
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
	   python::return_value_policy<python::reference_existing_object>())

      .def("GetConformers", GetMolConformers,
           "Get all the conformers as a tuple")

      .def("RemoveAllConformers", &ROMol::clearConformers,
           "Remove all the conformations on the molecule")

      .def("RemoveConformer", &ROMol::removeConformer,
           "Remove the conformer with the specified ID")
      .def("GetBondBetweenAtoms", 
	   (ROMol::GRAPH_EDGE_TYPE (ROMol::*)(unsigned int,unsigned int))&ROMol::getBondBetweenAtoms,
	   python::return_value_policy<python::reference_existing_object>(),
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
	    python::arg("useChirality")=false),
	   "Queries whether or not the molecule contains a particular substructure.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule\n\n"
	   "    - recursionPossible: (optional)\n\n"
	   "    - useChirality: (optional)\n\n"
	   "  RETURNS: 1 or 0\n")
      .def("GetSubstructMatch",GetSubstructMatch,
	   (python::arg("self"),python::arg("query"),
	    python::arg("useChirality")=false),
	   "Returns the indices of the molecule's atoms that match a substructure query.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule\n\n"
	   "    - useChirality: (optional)\n\n"
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
	    python::arg("useChirality")=false),
	   "Returns tuples of the indices of the molecule's atoms that match a substructure query.\n\n"
	   "  ARGUMENTS:\n"
	   "    - query: a Molecule.\n"
	   "    - uniquify: (optional) determines whether or not the matches are uniquified.\n"
	   "                Defaults to 1.\n\n"
	   "    - useChirality: (optional)\n\n"
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
	   "    - key: the name of the property to clear (a string).\n\n"
	   "  NOTE:\n"
	   "    - If the property has not been set, a KeyError exception will be raised.\n")

      .def("ClearComputedProps", MolClearComputedProps,
	   "Removes all computed properties from the molecule.\n\n")

      .def("UpdatePropertyCache", &ROMol::updatePropertyCache,
	   (python::arg("self"),python::arg("strict")=true),
	    "Regenerates computed properties like implicit valence and ring information.\n\n")


      .def("GetPropNames",MolGetPropNames,
	   (python::arg("self"),python::arg("includePrivate")=false,
	    python::arg("includeComputed")=false),
	   "Returns a tuple with all property names for this molecule.\n\n"
	   "  ARGUMENTS:\n"
	   "    - includePrivate: (optional) toggles inclusion of private properties in the result set.\n"
	   "                      Defaults to 0.\n"
	   "    - includeComputed: (optional) toggles inclusion of computed properties in the result set.\n"
	   "                      Defaults to 0.\n\n"
	   "  RETURNS: a tuple of strings\n")

#if 0
      .def("SetAtomBookmark",
	   (void (ROMol::*)(Atom *,int))&ROMol::setAtomBookmark,
	   "Sets an atom bookmark.\n\n"
	   "  ARGUMENTS:\n"
	   "    - atom: the Atom to bookmark\n"
	   "    - mark: an integer bookmark\n")

      .def("GetAtomWithBookmark",&ROMol::getAtomWithBookmark,
	   "Returns the Atom with a particular bookmark (if there is one).\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to be retrieved\n\n"
	   "  RETURNS:\n"
	   "    an Atom, if the bookmark exists. None otherwise.\n",
	   python::return_value_policy<python::reference_existing_object>())
      // This one makes g++ barf (internal compiler error), so skip it for now
      //.def("GetAllAtomsWithBookmark",&(ROMol::getAllAtomsWithBookmark))
      .def("ClearAtomBookmark",(void (ROMol::*)(int))&ROMol::clearAtomBookmark)
      .def("ClearAtomBookmark",(void (ROMol::*)(int,const Atom *))&ROMol::clearAtomBookmark,
	   "Clears an atom bookmark.\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to clear\n"
	   "    - atom: (optional) an Atom. If this argument is provided, the\n"
	   "            bookmark is only cleared on that particular atom.\n")
      .def("ClearAllAtomBookmarks",&ROMol::clearAllAtomBookmarks,
	   "Removes all atom bookmarks.\n\n")
      .def("HasAtomBookmark",&ROMol::hasAtomBookmark,
	   "Queries whether or not a particular atom bookmark is set.\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to look for\n\n"
	   "  RETURNS: 1 or 0\n")

	   
      .def("SetBondBookmark",(void (ROMol::*)(Bond *,int))&ROMol::setBondBookmark,
	   "Sets a bond bookmark.\n\n"
	   "  ARGUMENTS:\n"
	   "    - bond: the Bond to bookmark\n"
	   "    - mark: an integer bookmark\n")

      .def("GetBondWithBookmark",&ROMol::getBondWithBookmark,
	   "Returns the Bond with a particular bookmark (if there is one).\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to be retrieved\n\n"
	   "  RETURNS:\n"
	   "    a Bond, if the bookmark exists. None otherwise.\n",
	   python::return_value_policy<python::reference_existing_object>())
      .def("ClearBondBookmark",(void (ROMol::*)(int))&ROMol::clearBondBookmark)
      .def("ClearBondBookmark",(void (ROMol::*)(int, const Bond *))&ROMol::clearBondBookmark,
	   "Clears a bond bookmark.\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to clear\n"
	   "    - bond: (optional) an Bond. If this argument is provided, the\n"
	   "            bookmark is only cleared on that particular bond.\n")
      .def("ClearAllBondBookmarks",&ROMol::clearAllBondBookmarks,
	   "Removes all bond bookmarks.\n\n")

      .def("HasBondBookmark",&ROMol::hasBondBookmark,
	   "Queries whether or not a particular bond bookmark is set.\n\n"
	   "  ARGUMENTS:\n"
	   "    - mark: an integer with the bookmark to look for\n\n"
	   "  RETURNS: 1 or 0\n")

      // Iterators
      .def("GetAromaticAtoms",MolGetAromaticAtoms,
	   "Returns a read-only sequence containing all of the molecule's aromatic atoms.\n",
      	   python::return_value_policy<python::manage_new_object>())
      .def("GetHeteros",MolGetHeteros,
	   "Returns a read-only sequence containing all of the molecule's heteroatoms.\n",
      	   python::return_value_policy<python::manage_new_object>())
#endif
      .def("GetAtoms",MolGetAtoms,
	   "Returns a read-only sequence containing all of the molecule's Atoms.\n",
	   python::return_value_policy<python::manage_new_object>())
      .def("GetBonds",MolGetBonds,
	   "Returns a read-only sequence containing all of the molecule's Bonds.\n",
	   python::return_value_policy<python::manage_new_object>())


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
		 python::arg("useChirality")=false),
		"This function is included to speed substructure queries from databases, \n"
		"it's probably not of\n"
		"general interest.\n\n"
		"  ARGUMENTS:\n"
		"    - pkl: a Molecule pickle\n\n"
		"    - query: a Molecule\n\n"
		"    - recursionPossible: (optional)\n\n"
		"    - useChirality: (optional)\n\n"
		"  RETURNS: 1 or 0\n");

    
  };
};
}// end of namespace
void wrap_mol() {
  RDKit::mol_wrapper::wrap();
}
