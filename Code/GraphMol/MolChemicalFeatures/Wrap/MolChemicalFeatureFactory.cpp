// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include <boost/python.hpp>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDBoost/Exceptions.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <boost/python/converter/shared_ptr_to_python.hpp>

namespace python = boost::python;

namespace RDKit {
  // ----------------------------------------------------------------------------
  // NOTE: the reason we don't provide an interface to get all the
  // features from a molecule at the same time directly from C++ is a
  // potential problem with the lifetime of molecules and features:
  // Each Feature carries around a pointer to its owning molecule,
  //  and that pointer must remain valid. So we need to use the BPL's
  //  with_custodian_and_ward() functionality to tie the lifetime
  //  of the molecule to that of the feature. This is not currently
  //  possible in BPL v1.32 if we return a tuple of objects; so if 
  //  we return a set of features, there's no way to use
  //  with_custodian_and_ward() to ensure that the molecule doesn't vanish
  //  underneath them.  Access to all of a molecule's features is provided
  //  from within python (looping over getMolFeature() and gaining efficiency
  //  by not recomputing after the first call).
  // ----------------------------------------------------------------------------
  
  int getNumMolFeatures(const MolChemicalFeatureFactory &factory, 
			const ROMol &mol,
			std::string includeOnly="") {
    FeatSPtrList feats = factory.getFeaturesForMol(mol,includeOnly.c_str());
    return feats.size();
  }

  FeatSPtr getMolFeature(const MolChemicalFeatureFactory &factory, 
			 const ROMol &mol,
			 int idx,
			 std::string includeOnly="",
			 bool recompute=true) {
    static FeatSPtrList feats;
    if(recompute){
      feats = factory.getFeaturesForMol(mol,includeOnly.c_str());
    }
    if(idx<0 || idx>=static_cast<int>(feats.size())){
      throw IndexErrorException(idx);
    }

    FeatSPtrList_I fi=feats.begin();
    for(int i=0;i<idx;++i){
      ++fi;
    }
    return (*fi);
  }

  python::tuple getFeatureFamilies(const MolChemicalFeatureFactory &factory){
    python::list res;

    MolChemicalFeatureDef::CollectionType::const_iterator iter;
    for(iter=factory.beginFeatureDefs();iter!=factory.endFeatureDefs();++iter){
      std::string fam= (*iter)->getFamily();
      if(res.count(fam)==0) res.append(fam);
    }
    return python::tuple(res);
  }

  python::dict getFeatureDefs(const MolChemicalFeatureFactory &factory){
    python::dict res;
    MolChemicalFeatureDef::CollectionType::const_iterator iter;
    for(iter=factory.beginFeatureDefs();iter!=factory.endFeatureDefs();++iter){
      std::string key= (*iter)->getFamily()+"."+(*iter)->getType();
      res[key]=(*iter)->getSmarts();
    }
    return res;
  }
  

  struct featfactory_wrapper {
    static void wrap() {
      std::string factoryClassDoc="Class to featurize a molecule\n";
      // no direct constructor - use buildFactory function
      python::class_<MolChemicalFeatureFactory>("MolChemicalFeatureFactory",
						factoryClassDoc.c_str(), python::no_init)
        .def("GetNumFeatureDefs", &MolChemicalFeatureFactory::getNumFeatureDefs,
             "Get the number of feature definitions")
        .def("GetFeatureFamilies", getFeatureFamilies,
             "Get a tuple of feature types")
        .def("GetFeatureDefs", getFeatureDefs,
             "Get a dictionary with SMARTS definitions for each feature type")
        .def("GetNumMolFeatures", getNumMolFeatures,
	     (python::arg("mol"),python::arg("includeOnly")=std::string("")),
             "Get the number of features the molecule has")
        .def("GetMolFeature", getMolFeature,
	     (python::arg("mol"),python::arg("idx"),
	      python::arg("includeOnly")=std::string(""),
	      python::arg("recompute")=true),
	     python::with_custodian_and_ward_postcall<0,2>(),
             "returns a particular feature (by index)")
        ;
    };
  };
}

void wrap_factory() {
  RDKit::featfactory_wrapper::wrap();
}

