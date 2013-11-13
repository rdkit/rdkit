/*******************************************************************************
solutionInfo.cpp - Shape-it
 
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



#include <Shape/solutionInfo.h>



SolutionInfo::SolutionInfo(void)
:refName("")
, refAtomVolume(0.0)
, refCenter(0, 0, 0)
, refRotation(3, 3, 0)
, dbName("")
, dbAtomVolume(0.0)
, dbMol()
, dbCenter(0, 0, 0)
, dbRotation(3, 3, 0)
, atomOverlap(0.0)
, score(0.0)
, rotor(4, 0.0)
{
    rotor[0] = 1.0;
}



SolutionInfo::~SolutionInfo(void)
{
}


void
updateSolutionInfo(SolutionInfo & s, AlignmentInfo & res, double score,
		   GaussianVolume & gv)
{
    s.dbAtomVolume = gv.overlap;
    s.dbCenter = gv.centroid;
    s.dbRotation = gv.rotation;
    s.atomOverlap = res.overlap;
    s.score = score;
    s.rotor = res.rotor;
    return;
}



void setAllScores(SolutionInfo & res)
{
    RDKit::ROMol dbMol = *res.dbMol;

    dbMol.setProp(tanimoto,
		  res.atomOverlap / (res.refAtomVolume + res.dbAtomVolume -
				     res.atomOverlap));
    dbMol.setProp(tversky_ref,
		  res.atomOverlap / (0.95 * res.refAtomVolume +
				     0.05 * res.dbAtomVolume));
    dbMol.setProp(tversky_db,
		  res.atomOverlap / (0.05 * res.refAtomVolume +
				     0.95 * res.dbAtomVolume));

    return;
}
