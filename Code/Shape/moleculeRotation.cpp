/*******************************************************************************
moleculeRotation.cpp - Shape-it
 
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



#include <Shape/moleculeRotation.h>



void
positionMolecule(RDKit::ROMol& m, Coordinate& centroid, SiMath::Matrix& rotation)
{
   RDKit::ROMol::AtomIterator ai;
   RDKit::Conformer conf = m.getConformer(0);
   for (ai=m.beginAtoms(); ai!=m.endAtoms(); ai++)
   {
        RDKit::Atom *a = *ai;
        int idx = a->getIdx();
        
        // Get coordinates
        double x = conf.getAtomPos(idx).x - centroid.x;
        double y = conf.getAtomPos(idx).y - centroid.y;
        double z = conf.getAtomPos(idx).z - centroid.z;
		// Rotate according to eigenvectors SVD
        RDGeom::Point3D point;
        point.x = rotation[0][0] * x + rotation[1][0] * y + rotation[2][0] * z;
		point.y = rotation[0][1] * x + rotation[1][1] * y + rotation[2][1] * z;
		point.z = rotation[0][2] * x + rotation[1][2] * y + rotation[2][2] * z;
		
		m.getConformer(0).setAtomPos(idx, point);
	}
	return;
}


void
repositionMolecule(RDKit::ROMol& m, SiMath::Matrix& rotation, Coordinate& centroid)
{
   RDKit::ROMol::AtomIterator ai;
   RDKit::Conformer conf = m.getConformer(0);
   for (ai=m.beginAtoms(); ai!=m.endAtoms(); ai++)
   {
        RDKit::Atom *a = *ai;
        int idx = a->getIdx();
        
        // Get coordinates
        double x = conf.getAtomPos(idx).x;
        double y = conf.getAtomPos(idx).y;
        double z = conf.getAtomPos(idx).z;
		
		// Rotate according to eigenvectors SVD
		double xx = rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z;
		double yy = rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z;
		double zz = rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z;
		
		RDGeom::Point3D point;
        point.x = xx + centroid.x;
        point.y = yy + centroid.y;
        point.z = zz + centroid.z;
		m.getConformer(0).setAtomPos(idx, point);
	}
	return;
}



void 
rotateMolecule(RDKit::ROMol& m, SiMath::Vector& rotor)
{
	// Build rotation matrix
	SiMath::Matrix rot(3,3,0.0);
	double r1 = rotor[1] * rotor[1];
	double r2 = rotor[2] * rotor[2];
	double r3 = rotor[3] * rotor[3];
	
	rot[0][0] = 1.0 - 2.0*r2 - 2.0*r3;
	rot[0][1] = 2.0 * (rotor[1]*rotor[2] - rotor[0]*rotor[3]);
	rot[0][2] = 2.0 * (rotor[1]*rotor[3] + rotor[0]*rotor[2]);
	rot[1][0] = 2.0 * (rotor[1]*rotor[2] + rotor[0]*rotor[3]);
	rot[1][1] = 1.0 - 2*r3 - 2*r1;
	rot[1][2] = 2.0 * (rotor[2]*rotor[3] - rotor[0]*rotor[1]);
	rot[2][0] = 2.0 * (rotor[1]*rotor[3] - rotor[0]*rotor[2]);
	rot[2][1] = 2.0 * (rotor[2]*rotor[3] + rotor[0]*rotor[1]);
	rot[2][2] = 1.0 - 2*r2 - 2*r1;
   

   RDKit::ROMol::AtomIterator ai;
   RDKit::Conformer conf = m.getConformer(0);
   for (ai=m.beginAtoms(); ai!=m.endAtoms(); ai++)
   {
       
        RDKit::Atom *a = *ai;
        int idx = a->getIdx();
        
        if (idx == 0) {
            //std::cerr << conf.getAtomPos(0).x  << " " << m.getConformer(-1).getAtomPos(0).x << std::endl;
        }
        
        // Get coordinates
        double x = conf.getAtomPos(idx).x;
        double y = conf.getAtomPos(idx).y;
        double z = conf.getAtomPos(idx).z;
		
		// rotate according to eigenvectors SVD     
        RDGeom::Point3D point;
        point.x = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z;
        point.y = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z;
        point.z = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z;
        m.getConformer(0).setAtomPos(idx, point);

    }
	return;
}
