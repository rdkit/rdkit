// $Id$
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "MolChemicalFeature.h"
#include "MolChemicalFeatureDef.h"

namespace RDKit {

  MolChemicalFeatureDef::MolChemicalFeatureDef(const std::string &smarts,const std::string &family,
					 const std::string &type) : d_family(family),
								    d_type(type),
								    d_smarts(smarts) {
    ROMol *mol=static_cast<ROMol *>(SmartsToMol(smarts));
    dp_pattern.reset(mol);
  }
  void MolChemicalFeatureDef::normalizeWeights(){
    double accum=0.0;
    for(std::vector<double>::iterator i=this->beginWeights();
	i!=this->endWeights();i++){
      accum += *i;
    }
    for(std::vector<double>::iterator i=this->beginWeights();
	i!=this->endWeights();i++){
      *i /= accum;
    }
  }
  
}
