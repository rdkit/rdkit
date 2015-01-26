//
//  Copyright (C) 2015 Greg Landrum 
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <cmath>
#define M_PI 3.14159265358979323846

// ****************************************************************************

namespace RDKit {
  namespace MolDraw2D_detail {
    // implementation from $RDBASE/rdkit/sping/pid.py
    void arcPoints(const std::pair<float,float> &cds1,const std::pair<float,float> &cds2,
                   std::vector< std::pair<float,float> > &res,
                   float startAng,float extent){
      // Note: this implementation is simple and not particularly efficient.
      float xScale = (cds2.first-cds1.first)/2.0;
      float yScale = (cds2.second-cds1.second)/2.0;
      if(xScale<0) xScale *= -1;
      if(yScale<0) yScale *= -1;
        
      float x = std::min(cds1.first,cds2.first)+xScale;
      float y = std::min(cds1.second,cds2.second)+yScale;
        
      int steps = std::max(static_cast<int>(extent*2),5);
      float step = M_PI*extent/(180*steps);
      float angle = M_PI*startAng/180;
      for(unsigned int i=0;i<=steps;++i){
        std::pair<float,float> point = std::make_pair(x+xScale*cos(angle),
                                                      y-yScale*sin(angle));
        res.push_back(point);
        angle += step;
      }
    }
  }
}

