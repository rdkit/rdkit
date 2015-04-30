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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "MolChemicalFeature.h"
#include "MolChemicalFeatureDef.h"

namespace RDKit {

  const std::string &MolChemicalFeature::getFamily() const {
    return dp_def->getFamily(); //return d_family;
  };

  const std::string &MolChemicalFeature::getType() const {
    return dp_def->getType(); //return d_type;
  };

  void MolChemicalFeature::setActiveConformer(int confId) {
    PRECONDITION(dp_mol,"bad molecule");
    d_activeConf = confId;
  }

  RDGeom::Point3D MolChemicalFeature::getPos() const {
    return this->getPos(d_activeConf);
  }

  RDGeom::Point3D MolChemicalFeature::getPos(int confId) const {
    PRECONDITION(dp_mol,"bad molecule");
    PRECONDITION(dp_mol->getNumConformers(),"molecule has no conformers");
    if(confId==-1) confId = (*dp_mol->beginConformers())->getId();

    // ------------- 
    // Check to see if we've got the value cached:
    PointCacheType::const_iterator cacheIt=d_locs.find(confId);
    if(cacheIt != d_locs.end()) {
      return cacheIt->second;
    }

    // --------------
    // Nope, we have to figure it out on our own:
    RDGeom::Point3D res(0,0,0);
    bool setNeg1=false;
    if(confId==-1){
      setNeg1=true;
    }
    
    if(d_atoms.size()==1){
      res=dp_mol->getConformer(confId).getAtomPos((*d_atoms.begin())->getIdx());
    }else{
      PRECONDITION(dp_def,"bad definition");
      PRECONDITION(dp_def->getNumWeights()==this->getNumAtoms(),"weight/atom mismatch");
      std::vector<double>::const_iterator weightIt=dp_def->beginWeights();
      const Conformer &conf=dp_mol->getConformer(confId);
      for(AtomPtrContainer::const_iterator atomIt=d_atoms.begin();
	  atomIt!=d_atoms.end();atomIt++,weightIt++){
	const Atom *atom=*atomIt;
	RDGeom::Point3D p=conf.getAtomPos(atom->getIdx());
	p *= *weightIt;
	res += p;
      }
    }
    d_locs[confId]=res;

    return res;
  }
}
