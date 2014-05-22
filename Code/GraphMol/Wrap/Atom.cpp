// $Id$
//
//  Copyright (C) 2003-2013 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/QueryAtom.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/types.h>
#include <Geometry/point.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDBoost/Wrap.h>

#include "seqs.hpp"
#include <algorithm>


namespace python = boost::python;
namespace RDKit{
  namespace {
    std::string qhelper(Atom::QUERYATOM_QUERY *q,unsigned int depth){
      std::string res="";
      if(q){
	for (unsigned int i=0;i<depth;++i) res+="  ";
	res += q->getFullDescription()+"\n";
	for(Atom::QUERYATOM_QUERY::CHILD_VECT_CI ci=q->beginChildren();
	    ci!=q->endChildren();++ci){
	  res +=  qhelper((*ci).get(),depth+1);
	}
      }
      return res;
    }
  } // end of local namespace
  std::string describeQuery(const Atom*atom){
    std::string res="";
    if(atom->hasQuery()){
      res=qhelper(atom->getQuery(),0);
    }
    return res;
  }
  void expandQuery(QueryAtom *self,const QueryAtom *other,
                          Queries::CompositeQueryType how,
                          bool maintainOrder){
    if(other->hasQuery()){
      const QueryAtom::QUERYATOM_QUERY *qry=other->getQuery();
      self->expandQuery(qry->copy(),how,maintainOrder);
    }
  }

  void AtomSetProp(const Atom *atom, const char *key,std::string val) {
    //std::cerr<<"asp: "<<atom<<" " << key<<" - " << val << std::endl;
    atom->setProp(key, val);
  }
  
  int AtomHasProp(const Atom *atom, const char *key) {
    //std::cerr<<"ahp: "<<atom<<" " << key<< std::endl;
    int res = atom->hasProp(key);
    return res;
  }

  std::string AtomGetProp(const Atom *atom, const char *key) {
    if (!atom->hasProp(key)) {
      PyErr_SetString(PyExc_KeyError,key);
      throw python::error_already_set();
    }
    std::string res;
    atom->getProp(key, res);
    return res;
  }

  void AtomClearProp(const Atom *atom, const char *key) {
    if (!atom->hasProp(key)) {
      return;
    }
    atom->clearProp(key);
  }

  python::tuple AtomGetNeighbors(Atom *atom){
    python::list res;
    const ROMol *parent = &atom->getOwningMol();
    ROMol::ADJ_ITER begin,end;
    boost::tie(begin,end) = parent->getAtomNeighbors(atom);
    while(begin!=end){
      res.append(python::ptr(parent->getAtomWithIdx(*begin)));
      begin++;
    }
    return python::tuple(res);
  }

  python::tuple AtomGetBonds(Atom *atom){
    python::list res;
    const ROMol *parent = &atom->getOwningMol();
    ROMol::OEDGE_ITER begin,end;
    boost::tie(begin,end) = parent->getAtomBonds(atom);
    while(begin!=end){
      Bond *tmpB = (*parent)[*begin].get();
      res.append(python::ptr(tmpB));
      begin++;
    }
    return python::tuple(res);
  }

  bool AtomIsInRing(const Atom *atom){
    if(!atom->getOwningMol().getRingInfo()->isInitialized()){
      MolOps::findSSSR(atom->getOwningMol());
    }
    return atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())!=0;
  }
  bool AtomIsInRingSize(const Atom *atom,int size){
    if(!atom->getOwningMol().getRingInfo()->isInitialized()){
      MolOps::findSSSR(atom->getOwningMol());
    }
    return atom->getOwningMol().getRingInfo()->isAtomInRingOfSize(atom->getIdx(),size);
  }

  std::string AtomGetSmarts(const Atom *atom){
    std::string res;      
    if(atom->hasQuery()){
      res=SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
    } else {
      res = SmilesWrite::GetAtomSmiles(atom);
    }
    return res;
  }

  void SetAtomMonomerInfo(Atom *atom,const AtomMonomerInfo *info){
    atom->setMonomerInfo(info->copy());
  }

  AtomMonomerInfo *AtomGetMonomerInfo(Atom *atom){
    return atom->getMonomerInfo();
  }
  AtomPDBResidueInfo *AtomGetPDBResidueInfo(Atom *atom){
    AtomMonomerInfo *res=atom->getMonomerInfo();
    if(!res) return NULL;
    if(res->getMonomerType()!=AtomMonomerInfo::PDBRESIDUE){
      throw_value_error("MonomerInfo is not a PDB Residue");
    }
    return (AtomPDBResidueInfo *)res;
  }



  // FIX: is there any reason at all to not just prevent the construction of Atoms?
  std::string atomClassDoc="The class to store Atoms.\n\
Note that, though it is possible to create one, having an Atom on its own\n\
(i.e not associated with a molecule) is not particularly useful.\n";
struct atom_wrapper {
  static void wrap(){
    python::class_<Atom>("Atom",atomClassDoc.c_str(),python::init<std::string>())

      .def(python::init<unsigned int>("Constructor, takes either an int (atomic number) or a string (atomic symbol).\n"))

      .def("GetAtomicNum",&Atom::getAtomicNum,
	   "Returns the atomic number.")

      .def("SetAtomicNum",&Atom::setAtomicNum,
	   "Sets the atomic number, takes an integer value as an argument")

      .def("GetSymbol",&Atom::getSymbol,
	   "Returns the atomic symbol (a string)\n")

      .def("GetIdx",&Atom::getIdx,
	   "Returns the atom's index (ordering in the molecule)\n")

      .def("GetDegree",&Atom::getDegree,
	   "Returns the degree of the atom in the molecule.\n\n"
	   "  The degree of an atom is defined to be its number of\n"
	   "  directly-bonded neighbors.\n"
	   "  The degree is independent of bond orders, but is dependent\n"
	   "    on whether or not Hs are explicit in the graph.\n"
           )
      .def("GetTotalDegree",&Atom::getTotalDegree,
	   "Returns the degree of the atom in the molecule including Hs.\n\n"
	   "  The degree of an atom is defined to be its number of\n"
	   "  directly-bonded neighbors.\n"
	   "  The degree is independent of bond orders.\n")

      .def("GetTotalNumHs",&Atom::getTotalNumHs,
           (python::arg("self"),python::arg("includeNeighbors")=false),
           "Returns the total number of Hs (explicit and implicit) on the atom.\n\n"
           "  ARGUMENTS:\n\n"
           "    - includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.\n"
           "      Defaults to 0.\n")
      .def("GetNumImplicitHs",&Atom::getNumImplicitHs,
	   "Returns the total number of implicit Hs on the atom.\n")

      .def("GetExplicitValence",&Atom::getExplicitValence,
           "Returns the number of explicit Hs on the atom.\n")
      .def("GetImplicitValence",&Atom::getImplicitValence,
           "Returns the number of implicit Hs on the atom.\n")
      .def("GetTotalValence",&Atom::getTotalValence,
	   "Returns the total valence (explicit + implicit) of the atom.\n\n")

      .def("GetFormalCharge",&Atom::getFormalCharge)
      .def("SetFormalCharge",&Atom::setFormalCharge)


      .def("SetNoImplicit",&Atom::setNoImplicit,
	   "Sets a marker on the atom that *disallows* implicit Hs.\n"
	   "  This holds even if the atom would otherwise have implicit Hs added.\n")
      .def("GetNoImplicit",&Atom::getNoImplicit,
	   "Returns whether or not the atom is *allowed* to have implicit Hs.\n")

      .def("SetNumExplicitHs",&Atom::setNumExplicitHs)
      .def("GetNumExplicitHs",&Atom::getNumExplicitHs)
      .def("SetIsAromatic",&Atom::setIsAromatic)
      .def("GetIsAromatic",&Atom::getIsAromatic)
      .def("SetMass",&Atom::setMass)
      .def("GetMass",&Atom::getMass)
      .def("SetIsotope",&Atom::setIsotope)
      .def("GetIsotope",&Atom::getIsotope)
      .def("SetNumRadicalElectrons",&Atom::setNumRadicalElectrons)
      .def("GetNumRadicalElectrons",&Atom::getNumRadicalElectrons)

      // NOTE: these may be used at some point in the future, but they
      //  aren't now, so there's no point in confusing things.
      //.def("SetDativeFlag",&Atom::setDativeFlag)
      //.def("GetDativeFlag",&Atom::getDativeFlag)
      //.def("ClearDativeFlag",&Atom::clearDativeFlag)

      .def("SetChiralTag",&Atom::setChiralTag)
      .def("InvertChirality",&Atom::invertChirality)
      .def("GetChiralTag",&Atom::getChiralTag)

      .def("SetHybridization",&Atom::setHybridization,
	   "Sets the hybridization of the atom.\n"
	   "  The argument should be a HybridizationType\n")
      .def("GetHybridization",&Atom::getHybridization,
	   "Returns the atom's hybridization.\n")

      .def("GetOwningMol",&Atom::getOwningMol,
	   "Returns the Mol that owns this atom.\n",
	   python::return_value_policy<python::reference_existing_object>())

      .def("GetNeighbors",AtomGetNeighbors,
	   "Returns a read-only sequence of the atom's neighbors\n")

      .def("GetBonds",AtomGetBonds,
	   "Returns a read-only sequence of the atom's bonds\n")

      .def("Match",(bool (Atom::*)(const Atom *) const)&Atom::Match,
	   "Returns whether or not this atom matches another Atom.\n\n"
	   "  Each Atom (or query Atom) has a query function which is\n"
	   "  used for this type of matching.\n\n"
	   "  ARGUMENTS:\n"
	   "    - other: the other Atom to which to compare\n")

      .def("IsInRingSize",AtomIsInRingSize,
	   "Returns whether or not the atom is in a ring of a particular size.\n\n"
	   "  ARGUMENTS:\n"
	   "    - size: the ring size to look for\n") 

      .def("IsInRing",AtomIsInRing,
	   "Returns whether or not the atom is in a ring\n\n")

      .def("HasQuery",&Atom::hasQuery,
     "Returns whether or not the atom has an associated query\n\n")

      .def("DescribeQuery",describeQuery,
	   "returns a text description of the query. Primarily intended for debugging purposes.\n\n")

      .def("GetSmarts",AtomGetSmarts,
              "returns the SMARTS (or SMILES) string for an Atom\n\n")

      // properties
      .def("SetProp",AtomSetProp,
	   (python::arg("self"), python::arg("key"),
	    python::arg("val")),
	   "Sets an atomic property\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to be set (a string).\n"
	   "    - value: the property value (a string).\n\n"
           )

      .def("GetProp", AtomGetProp,
           "Returns the value of the property.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to return (a string).\n\n"
	   "  RETURNS: a string\n\n"
	   "  NOTE:\n"
	   "    - If the property has not been set, a KeyError exception will be raised.\n")

      .def("HasProp", AtomHasProp,
           "Queries a Atom to see if a particular property has been assigned.\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to check for (a string).\n")

      .def("ClearProp", AtomClearProp,
           "Removes a particular property from an Atom (does nothing if not already set).\n\n"
	   "  ARGUMENTS:\n"
	   "    - key: the name of the property to be removed.\n")

      .def("GetPropNames",&Atom::getPropList,
	   (python::arg("self")),
           "Returns a list of the properties set on the Atom.\n\n"
           )

      .def("GetMonomerInfo",
	   AtomGetMonomerInfo,
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns the atom's MonomerInfo object, if there is one.\n\n"
	   )
      .def("GetPDBResidueInfo",
	   AtomGetPDBResidueInfo,
	   python::return_internal_reference<1,
	   python::with_custodian_and_ward_postcall<0,1> >(),
	   "Returns the atom's MonomerInfo object, if there is one.\n\n"
	   )
      .def("SetMonomerInfo", 
	   SetAtomMonomerInfo,
	   "Sets the atom's MonomerInfo object.\n\n"
	   )
      ;

    python::enum_<Atom::HybridizationType>("HybridizationType")
      .value("UNSPECIFIED",Atom::UNSPECIFIED)
      .value("SP",Atom::SP)
      .value("SP2",Atom::SP2)
      .value("SP3",Atom::SP3)
      .value("SP3D",Atom::SP3D)
      .value("SP3D2",Atom::SP3D2)
      .value("OTHER",Atom::OTHER)
      ;
    python::enum_<Atom::ChiralType>("ChiralType")
      .value("CHI_UNSPECIFIED",Atom::CHI_UNSPECIFIED)
      .value("CHI_TETRAHEDRAL_CW",Atom::CHI_TETRAHEDRAL_CW)
      .value("CHI_TETRAHEDRAL_CCW",Atom::CHI_TETRAHEDRAL_CCW)
      .value("CHI_OTHER",Atom::CHI_OTHER)
      ;


    python::enum_<Queries::CompositeQueryType>("CompositeQueryType")
      .value("COMPOSITE_AND",Queries::COMPOSITE_AND)
      .value("COMPOSITE_OR",Queries::COMPOSITE_OR)
      .value("COMPOSITE_XOR",Queries::COMPOSITE_XOR)
      ;
      

    atomClassDoc="The class to store QueryAtoms.\n\
These cannot currently be constructed directly from Python\n";
    python::class_<QueryAtom,python::bases<Atom> >("QueryAtom",atomClassDoc.c_str(),python::no_init)
      .def("ExpandQuery",expandQuery,
           (python::arg("self"),python::arg("other"),python::arg("how")=Queries::COMPOSITE_AND,
            python::arg("maintainOrder")=true),
              "combines the query from other with ours")
      ;

  };
};
}// end of namespace
void wrap_atom() {
  RDKit::atom_wrapper::wrap();
}
