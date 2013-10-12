/*******************************************************************************
shapeAlignment.h - Shape-it
 
Copyright 2012 by Silicos-it, a division of Imacosi BVBA
 
This file is part of Shape-it.

	Shape-it is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published 
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Shape-it is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with Shape-it.  If not, see <http://www.gnu.org/licenses/>.

Shape-it is linked against OpenBabel version 2.

	OpenBabel is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

***********************************************************************/



#ifndef __SILICOSIT_SHAPEIT_SHAPEALIGNMENT_H__
#define __SILICOSIT_SHAPEIT_SHAPEALIGNMENT_H__



// General
#include <queue>
#include <vector>

// OpenBabel

// Shape-it
#include <Shape/gaussianVolume.h>
#include <Shape/atomGaussian.h>
#include <Shape/alignmentInfo.h>
#include <Shape/siMath.h>



typedef std::map<unsigned int, double *> MatrixMap;
typedef std::map<unsigned int, double *>::iterator MatIter;
	
   
   
class ShapeAlignment 
{
   private:
      
      GaussianVolume * _gRef;
      GaussianVolume * _gDb;
				
      unsigned int _rAtoms;
      unsigned int _rGauss;
      unsigned int _dAtoms;
      unsigned int _dGauss;
      unsigned int _maxSize;
      unsigned int _maxIter;
			
      MatrixMap _matrixMap;
		
      double* _updateMatrixMap(AtomGaussian&, AtomGaussian&);
			
   public:
      
      ShapeAlignment(GaussianVolume&, GaussianVolume&);
      ~ShapeAlignment(void);
			
      AlignmentInfo gradientAscent(SiMath::Vector);
      AlignmentInfo simulatedAnnealing(SiMath::Vector);
			
      void setMaxIterations(unsigned int);
};


#endif