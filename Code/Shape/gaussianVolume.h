/*******************************************************************************
gaussianVolume.h - Shape-it
 
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



#ifndef __SILICOSIT_SHAPEIT_GAUSSIANVOLUME_H__
#define __SILICOSIT_SHAPEIT_GAUSSIANVOLUME_H__



// General
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>
#include <queue>

// RDKit
#include <GraphMol/ROMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/AtomIterators.h>
#include <Geometry/point.h>
#include <Numerics/Matrix.h>
#include <Numerics/SymmMatrix.h>
#include <Numerics/Alignment/AlignPoints.h>


// Shape-it
#include <Shape/siMath.h>
#include <Shape/atomGaussian.h>
#include <Shape/alignmentInfo.h>




const double GCI = 2.828427125;
const double GLI = 1.480960979;
const double VCUTOFF = 2.0;
const unsigned int LEVEL = 6;
const double EPS = 0.03;
const double GRADSCALE = 0.9;
const double PENALTY = 5.00;



class GaussianVolume {
public:
  double volume;		///< Molecular volume
  double overlap;		///< Self-overlap of the molecule
  RDGeom::Point3D centroid;	///< center of the gaussian volume

  std::vector < AtomGaussian > gaussians;	///< vector of all atom gaussians and their overlaps
  std::vector < std::vector < unsigned int >*>childOverlaps;	///< vector to keep track of which overlaps are formed with one gaussian
  std::vector < unsigned int >levels;	///< indicates where in the vector the level of overlaps changes
  double rotation[3][3];        ///< rotation matrix to align molecule to principal axes

  GaussianVolume(void);
  ~GaussianVolume(void);
};



void listAtomVolumes(RDKit::ROMol & mol, GaussianVolume & gv);
void initOrientation(GaussianVolume &);
double atomOverlap(GaussianVolume &, GaussianVolume &);
double GAlpha(unsigned int);
double getScore(std::string &, double, double, double);
void checkVolumes(GaussianVolume &, GaussianVolume &, AlignmentInfo &);



#endif
