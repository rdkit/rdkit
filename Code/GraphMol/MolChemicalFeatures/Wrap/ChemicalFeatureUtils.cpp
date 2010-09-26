// $Id$
//
//  Copyright (C) 2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <boost/dynamic_bitset.hpp>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>

namespace python = boost::python;
namespace RDKit {
  python::object GetAtomMatch(python::object featMatch,int maxAts=1024){
    python::list res;
    unsigned int nEntries=python::extract<unsigned int>(featMatch.attr("__len__")());
    
    boost::dynamic_bitset<> indices(maxAts);
    for(unsigned int i=0;i<nEntries;++i){
      MolChemicalFeature *feat=python::extract<MolChemicalFeature *>(featMatch[i]);
      const MolChemicalFeature::AtomPtrContainer &atoms = feat->getAtoms();
      MolChemicalFeature::AtomPtrContainer_CI aci;
      python::list local;
      for(aci=atoms.begin();aci!=atoms.end();++aci){
	unsigned int idx=(*aci)->getIdx();
	if(indices[idx]){
	  return python::list();
	} else {
	  indices[idx]=1;
	}
	local.append(idx);
      }
      res.append(local);
    }
    return res;
  }

#if 0
  // defined in MolChemicalFeatureFactory.cpp
  python::tuple getFeatsForMol(const MolChemicalFeatureFactory &factory, 
			   const ROMol &mol,
			   std::string includeOnly);

  // Moving this into C++ didn't provide any speed improvement in our profiling,
  // so, rather than adding additional complications to the wrapper interface,
  // we'll just leave it out
  //
  // In the event it should be re-integrated, here's the def:
  //python::def("GetFeatureMatchDict",GetFeatMatchDict,(python::arg("mol"),
  //						      python::arg("factory"),
  //						      python::arg("feats")),
  //           );

  python::dict GetFeatMatchDict(const ROMol *mol,const MolChemicalFeatureFactory *factory,
				python::object feats){
    python::dict res;
    unsigned int nEntries=python::extract<unsigned int>(feats.attr("__len__")());
    
    for(unsigned int i=0;i<nEntries;++i){
      const char *family;
      python::extract<MolChemicalFeature *> cfX(feats[i]);
      if(cfX.check()){
	family=cfX()->getFamily().c_str();
      } else {
	python::extract<ChemicalFeatures::FreeChemicalFeature *> fcfX(feats[i]);
	family=fcfX()->getFamily().c_str();
      }
      if(!res.has_key(family)){
	res[family] = getFeatsForMol(factory,mol,family);
      }
    }
    return res;
  }


#endif  
  struct ChemicalFeatureUtils_wrapper {
    static void wrap() {
      python::def("GetAtomMatch",GetAtomMatch,(python::arg("featMatch"),python::arg("maxAts")=1024),
		  "Returns an empty list if any of the features passed in share an atom.\n\
 Otherwise a list of lists of atom indices is returned.\n");
    }
  };
} // end of namespace RDKit

void wrap_ChemicalFeatureUtils() {
  RDKit::ChemicalFeatureUtils_wrapper::wrap();
}

