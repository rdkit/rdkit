// $Id$
//
//  Copyright (C) 2003-2014 Rational Discovery LLC
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

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/types.h>

namespace python = boost::python;
namespace RDKit{

  namespace {
    std::string qhelper(Bond::QUERYBOND_QUERY *q,unsigned int depth){
      std::ostringstream res;

      if(q){
	for (unsigned int i=0;i<depth;++i) res<<"  ";
	res << q->getFullDescription();
        if(q->getNegation()) res <<" ! ";
        
        res << "\n";
	for(Bond::QUERYBOND_QUERY::CHILD_VECT_CI ci=q->beginChildren();
	    ci!=q->endChildren();++ci){
	  res << qhelper((*ci).get(),depth+1);
	}
      }
      return res.str();
    }
  } // end of local namespace
  std::string describeQuery(const Bond *bond){
    std::string res="";
    if(bond->hasQuery()){
      res=qhelper(bond->getQuery(),0);
    }
    return res;
  }



  int BondHasProp(const Bond *bond, const char *key) {
    int res = bond->hasProp(key);
    return res;
  }

  void BondSetProp(const Bond *bond, const char *key,std::string val) {
    bond->setProp(key, val);
  }
  void BondClearProp(const Bond *bond, const char *key){
    if (!bond->hasProp(key)) {
      return;
    }
    bond->clearProp(key);
  }


  std::string BondGetProp(const Bond *bond, const char *key) {
    if (!bond->hasProp(key)) {
      PyErr_SetString(PyExc_KeyError,key);
      throw python::error_already_set();
    }
    std::string res;
    bond->getProp(key, res);
    return res;
  }

  bool BondIsInRing(const Bond *bond){
    if(!bond->getOwningMol().getRingInfo()->isInitialized()){
      MolOps::findSSSR(bond->getOwningMol());
    }
    return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx())!=0;

  }

  bool BondIsInRingSize(const Bond *bond,int size){
    if(!bond->getOwningMol().getRingInfo()->isInitialized()){
      MolOps::findSSSR(bond->getOwningMol());
    }
    return bond->getOwningMol().getRingInfo()->isBondInRingOfSize(bond->getIdx(),size);
    return false;
  }

  INT_VECT getBondStereoAtoms(const Bond *bond){
    return bond->getStereoAtoms();
  }

  std::string BondGetSmarts(const Bond *bond,bool allBondsExplicit){
    std::string res;
    if(bond->hasQuery()){      
      res=SmartsWrite::GetBondSmarts(static_cast<const QueryBond *>(bond));
    } else {
      res=SmilesWrite::GetBondSmiles(bond,-1,false,allBondsExplicit);
    }
    return res;
  }

  
  std::string bondClassDoc="The class to store Bonds.\n\
Note: unlike Atoms, is it currently impossible to construct Bonds from\n\
Python.\n";
struct bond_wrapper {
  static void wrap(){
    python::class_<Bond>("Bond",bondClassDoc.c_str(),python::no_init)

      .def("GetOwningMol",&Bond::getOwningMol,
	   "Returns the Mol that owns this bond.\n",
	   python::return_value_policy<python::reference_existing_object>())

      .def("GetBondType",&Bond::getBondType,
	   "Returns the type of the bond as a BondType\n")
      .def("SetBondType",&Bond::setBondType,
	   "Set the type of the bond as a BondType\n")
      .def("GetBondTypeAsDouble",&Bond::getBondTypeAsDouble,
	   "Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)\n")

      .def("GetBondDir",&Bond::getBondDir,
	   "Returns the type of the bond as a BondDir\n")
      .def("SetBondDir",&Bond::setBondDir,
	   "Set the type of the bond as a BondDir\n")

      .def("GetStereo",&Bond::getStereo,
	   "Returns the CIP-classification of the bond as a BondStereo\n")
      // this is no longer exposed because it requires that stereo atoms
      // be set. This is a task that is tricky and "dangerous".
      //.def("SetStereo",&Bond::setStereo,
      //	   "Set the CIP-classification of the bond as a BondStereo\n")
      .def("GetStereoAtoms",getBondStereoAtoms,
	   "Returns the indices of the atoms setting this bond's stereochemistry.")
      .def("GetValenceContrib",
	   (double (Bond::*)(const Atom *) const)&Bond::getValenceContrib,
	   "Returns the contribution of the bond to the valence of an Atom.\n\n"
	   "  ARGUMENTS:\n\n"
	   "    - atom: the Atom to consider.\n")

      .def("GetIsAromatic",&Bond::getIsAromatic)
      .def("SetIsAromatic",&Bond::setIsAromatic)

      .def("GetIsConjugated",&Bond::getIsConjugated,
	   "Returns whether or not the bond is considered to be conjugated.")
      .def("SetIsConjugated",&Bond::setIsConjugated)

      .def("GetIdx",&Bond::getIdx,
	   "Returns the bond's index (ordering in the molecule)\n")

      .def("GetBeginAtomIdx",&Bond::getBeginAtomIdx,
	   "Returns the index of the bond's first atom.\n")
      .def("GetEndAtomIdx",&Bond::getEndAtomIdx,
	   "Returns the index of the bond's first atom.\n")
      .def("GetOtherAtomIdx",&Bond::getOtherAtomIdx,
	   "Given the index of one of the bond's atoms, returns the\n"
	   "index of the other.\n")

      .def("GetBeginAtom",&Bond::getBeginAtom,
	   python::return_value_policy<python::reference_existing_object>(),
	   "Returns the bond's first atom.\n")
      .def("GetEndAtom",&Bond::getEndAtom,
	   python::return_value_policy<python::reference_existing_object>(),
	   "Returns the bond's second atom.\n")
      .def("GetOtherAtom",&Bond::getOtherAtom,
	   python::return_value_policy<python::reference_existing_object>(),
	   "Given one of the bond's atoms, returns the other one.\n")

      // FIX: query stuff
      .def("Match",(bool (Bond::*)(const Bond *) const)&Bond::Match,
	   "Returns whether or not this bond matches another Bond.\n\n"
	   "  Each Bond (or query Bond) has a query function which is\n"
	   "  used for this type of matching.\n\n"
	   "  ARGUMENTS:\n"
	   "    - other: the other Bond to which to compare\n")

      .def("IsInRingSize",BondIsInRingSize,
	   "Returns whether or not the bond is in a ring of a particular size.\n\n"
	   "  ARGUMENTS:\n"
	   "    - size: the ring size to look for\n")
      .def("IsInRing",BondIsInRing,
	   "Returns whether or not the bond is in a ring of any size.\n\n")

      .def("HasQuery",&Bond::hasQuery,
           "Returns whether or not the bond has an associated query\n\n")

      .def("DescribeQuery",describeQuery,
	   "returns a text description of the query. Primarily intended for debugging purposes.\n\n")

      .def("GetSmarts",BondGetSmarts,
           (python::arg("bond"),
            python::arg("allBondsExplicit")=false),
              "returns the SMARTS (or SMILES) string for a Bond")

      // properties
      .def("SetProp",BondSetProp,
	   (python::arg("self"), python::arg("key"),
	    python::arg("val")),
	   "Sets a bond property\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to be set (a string).\n"
	   "    - value: the property value (a string).\n\n"
           )

      .def("GetProp", BondGetProp,
           "Returns the value of the property.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to return (a string).\n\n"
	   "  RETURNS: a string\n\n"
	   "  NOTE:\n"
	   "    - If the property has not been set, a KeyError exception will be raised.\n")

      .def("HasProp", BondHasProp,
           "Queries a Bond to see if a particular property has been assigned.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to check for (a string).\n")

      .def("ClearProp", BondClearProp,
           "Removes a particular property from an Bond (does nothing if not already set).\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to be removed.\n")

      .def("GetPropNames",&Bond::getPropList,
	   (python::arg("self")),
           "Returns a list of the properties set on the Bond.\n\n"
           )
      
      ;

    python::enum_<Bond::BondType>("BondType")
      .value("UNSPECIFIED",Bond::UNSPECIFIED)
      .value("SINGLE",Bond::SINGLE)
      .value("DOUBLE",Bond::DOUBLE)
      .value("TRIPLE",Bond::TRIPLE)
      .value("QUADRUPLE",Bond::QUADRUPLE)
      .value("QUINTUPLE",Bond::QUINTUPLE)
      .value("HEXTUPLE",Bond::HEXTUPLE)
      .value("ONEANDAHALF",Bond::ONEANDAHALF)
      .value("TWOANDAHALF",Bond::TWOANDAHALF)
      .value("THREEANDAHALF",Bond::THREEANDAHALF)
      .value("FOURANDAHALF",Bond::FOURANDAHALF)
      .value("FIVEANDAHALF",Bond::FIVEANDAHALF)
      .value("AROMATIC",Bond::AROMATIC)
      .value("IONIC",Bond::IONIC)
      .value("HYDROGEN",Bond::HYDROGEN)
      .value("THREECENTER",Bond::THREECENTER)
      .value("DATIVEONE",Bond::DATIVEONE)
      .value("DATIVE",Bond::DATIVE)
      .value("DATIVEL",Bond::DATIVEL)
      .value("DATIVER",Bond::DATIVER)
      .value("OTHER",Bond::OTHER)
      .value("ZERO",Bond::ZERO)
      ;
    python::enum_<Bond::BondDir>("BondDir")
      .value("NONE",Bond::NONE)
      .value("BEGINWEDGE",Bond::BEGINWEDGE)
      .value("BEGINDASH",Bond::BEGINDASH)
      .value("ENDDOWNRIGHT",Bond::ENDDOWNRIGHT)
      .value("ENDUPRIGHT",Bond::ENDUPRIGHT)
      .value("UNKNOWN",Bond::UNKNOWN)
      ;
    python::enum_<Bond::BondStereo>("BondStereo")
      .value("STEREONONE",Bond::STEREONONE)
      .value("STEREOANY",Bond::STEREOANY)
      .value("STEREOZ",Bond::STEREOZ)
      .value("STEREOE",Bond::STEREOE)
      ;
  };
  
};
}// end of namespace
void wrap_bond() {
  RDKit::bond_wrapper::wrap();
}
