// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
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

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <boost/foreach.hpp>

namespace RDDepict {
  void compute2DCoordsForReaction(RDKit::ChemicalReaction &rxn,
                                  double spacing,
                                  bool updateProps,
                                  bool canonOrient,
                                  unsigned int nFlipsPerSample,
                                  unsigned int nSamples,
                                  int sampleSeed,
                                  bool permuteDeg4Nodes){
    double xOffset=0.0;
    for(RDKit::MOL_SPTR_VECT::iterator templIt=rxn.beginReactantTemplates();
        templIt!=rxn.endReactantTemplates();++templIt){
      if(updateProps){
        (*templIt)->updatePropertyCache(false);
        RDKit::MolOps::setConjugation(**templIt);
        RDKit::MolOps::setHybridization(**templIt);
      }
      compute2DCoords(**templIt,0,canonOrient,true,nFlipsPerSample,
                      nSamples,sampleSeed,permuteDeg4Nodes);
      double minX=100.,maxX=-100.;
      BOOST_FOREACH(RDGeom::Point3D &pt,(*templIt)->getConformer().getPositions()){
        minX=std::min(pt.x,minX);
      }
      xOffset += minX;
      BOOST_FOREACH(RDGeom::Point3D &pt,(*templIt)->getConformer().getPositions()){
        pt.x += xOffset;
        maxX=std::max(pt.x,maxX);
      }
      xOffset=maxX+spacing;
    }
    for(RDKit::MOL_SPTR_VECT::iterator templIt=rxn.beginProductTemplates();
        templIt!=rxn.endProductTemplates();++templIt){
      if(updateProps){
        (*templIt)->updatePropertyCache(false);
        RDKit::MolOps::setConjugation(**templIt);
        RDKit::MolOps::setHybridization(**templIt);
      }
      compute2DCoords(**templIt,0,canonOrient,true,nFlipsPerSample,
                      nSamples,sampleSeed,permuteDeg4Nodes);
      double minX=100.,maxX=-100.;
      BOOST_FOREACH(RDGeom::Point3D &pt,(*templIt)->getConformer().getPositions()){
        minX=std::min(pt.x,minX);
      }
      xOffset += minX;
      BOOST_FOREACH(RDGeom::Point3D &pt,(*templIt)->getConformer().getPositions()){
        pt.x += xOffset;
        maxX=std::max(pt.x,maxX);
      }
      xOffset=maxX+spacing;
    }

  }
} // end of namespace RDKit
    
    
