/*******************************************************************************
gaussianVolume.cpp - Shape-it
 
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



#include <Shape/gaussianVolume.h>
#include <Shape/solutionInfo.h>
#include <eigen3/Eigen/Dense>


GaussianVolume::GaussianVolume(void)
:volume(0.0)
, overlap(0.0)
, centroid(0.0, 0.0, 0.0)
, gaussians()
, childOverlaps()
, levels()
{
}



GaussianVolume::~GaussianVolume(void)
{
}



double GAlpha(unsigned int an)
{
    switch (an) {
    case 1:			///< H
	return 1.679158285;
	break;
    case 3:			///< Li
	return 0.729980658;
	break;
    case 5:			///< B
	return 0.604496983;
	break;
    case 6:			///< C
	return 0.836674025;
	break;
    case 7:			///< N
	return 1.006446589;
	break;
    case 8:			///< O
	return 1.046566798;
	break;
    case 9:			///< F
	return 1.118972618;
	break;
    case 11:			///< Na
	return 0.469247983;
	break;
    case 12:			///< Mg
	return 0.807908026;
	break;
    case 14:			///< Si
	return 0.548296583;
	break;
    case 15:			///< P
	return 0.746292571;
	break;
    case 16:			///< S
	return 0.746292571;
	break;
    case 17:			///< Cl
	return 0.789547080;
	break;
    case 19:			///< K
	return 0.319733941;
	break;
    case 20:			///< Ca
	return 0.604496983;
	break;
    case 26:			///< Fe
	return 1.998337133;
	break;
    case 29:			///< Cu
	return 1.233667312;
	break;
    case 30:			///< Zn
	return 1.251481772;
	break;
    case 35:			///< Br
	return 0.706497569;
	break;
    case 53:			///< I
	return 0.616770720;
	break;
    default:			///< *
	return 1.074661303;
    }
    return 1.074661303;
};



void listAtomVolumes(RDKit::ROMol & mol, GaussianVolume & gv)
{
    // Prepare the vector to store the atom and overlap volumes;
    unsigned int N(0);

    RDKit::ROMol::AtomIterator ai0;
    for (ai0 = mol.beginAtoms(); ai0 != mol.endAtoms(); ai0++) {
	RDKit::Atom * a = *ai0;
	int idx = a->getIdx();

	if (a->getAtomicNum() == 1) {
	    continue;
	} else {
	    ++N;
	}
    }
    gv.gaussians.resize(N);
    gv.childOverlaps.resize(N);
    gv.levels.push_back(N);	// first level
    gv.volume = 0.0;
    gv.centroid.x = 0;
    gv.centroid.y = 0;
    gv.centroid.z = 0;

    // Vector to keep track of parents of an overlap
    std::vector < std::pair < unsigned int, unsigned int > >parents(N);

    // Create a vector to keep track of the overlaps 
    // Overlaps are stored as sets
    std::vector < std::set < unsigned int >*>overlaps(N);
    std::set < unsigned int >::iterator setIter;

    // Start by iterating over the single atoms and build map of overlaps
    const RDKit::PeriodicTable * tbl = RDKit::PeriodicTable::getTable();

    int atomIndex = 0;		// keeps track of the atoms processed so far
    int vecIndex = N;		// keeps track of the last element added to the vectors
    RDKit::ROMol::AtomIterator ai1;
    RDKit::Conformer conf = mol.getConformer(0);
    for (ai1 = mol.beginAtoms(); ai1 != mol.endAtoms(); ai1++) {
	RDKit::Atom * a = *ai1;
	int idx = a->getIdx();
	// Skip hydrogens
	if (a->getAtomicNum() == 1) {
	    continue;
	}
	// First atom append self to the list
	// Store it at [index] 
	gv.gaussians[atomIndex].center.x = conf.getAtomPos(idx).x;
	gv.gaussians[atomIndex].center.y = conf.getAtomPos(idx).y;
	gv.gaussians[atomIndex].center.z = conf.getAtomPos(idx).z;
	gv.gaussians[atomIndex].alpha = GAlpha(a->getAtomicNum());
	gv.gaussians[atomIndex].C = GCI;
	double radius = tbl->getRvdw(a->getAtomicNum());
	gv.gaussians[atomIndex].volume =
	    (4.0 * PI / 3.0) * radius * radius * radius;
	gv.gaussians[atomIndex].nbr = 1;

	// Add empty child overlaps
	std::vector < unsigned int >*vec =
	    new std::vector < unsigned int >();
	gv.childOverlaps[atomIndex] = vec;

	// Update volume and centroid
	gv.volume += gv.gaussians[atomIndex].volume;
	gv.centroid.x +=
	    gv.gaussians[atomIndex].volume *
	    gv.gaussians[atomIndex].center.x;
	gv.centroid.y +=
	    gv.gaussians[atomIndex].volume *
	    gv.gaussians[atomIndex].center.y;
	gv.centroid.z +=
	    gv.gaussians[atomIndex].volume *
	    gv.gaussians[atomIndex].center.z;

	// Add new empty set of possible overlaps
	std::set < unsigned int >*tmp = new std::set < unsigned int >();
	overlaps[atomIndex] = tmp;

	// Loop over the current list of processed atoms and add overlaps
	for (int i = 0; i < atomIndex; ++i) {
	    // Create overlap gaussian
	    AtomGaussian ga =
		atomIntersection(gv.gaussians[i], gv.gaussians[atomIndex]);

	    // Check if the atom-atom overlap volume is large enough
	    if (ga.volume /
		(gv.gaussians[i].volume + gv.gaussians[atomIndex].volume -
		 ga.volume) < EPS) {
		continue;
	    }
	    // Add gaussian volume, and empty overlap set
	    gv.gaussians.push_back(ga);
	    std::vector < unsigned int >*vec =
		new std::vector < unsigned int >();
	    gv.childOverlaps.push_back(vec);

	    // Update local variables of parents and possible overlaps
	    parents.push_back(std::make_pair < unsigned int,
			      unsigned int >(i, atomIndex));
	    std::set < unsigned int >*dummy =
		new std::set < unsigned int >();
	    overlaps.push_back(dummy);

	    // Update volume and centroid (negative contribution of atom-atom overlap)
	    gv.volume -= ga.volume;
	    gv.centroid.x -= ga.volume * ga.center.x;
	    gv.centroid.y -= ga.volume * ga.center.y;
	    gv.centroid.z -= ga.volume * ga.center.z;

	    // Update overlap information of the parent atom
	    overlaps[i]->insert(atomIndex);
	    gv.childOverlaps[i]->push_back(vecIndex);

	    // Move to next index in vector
	    ++vecIndex;
	}

	// Update atom index
	++atomIndex;
    }

    // Position in list of gaussians where atom gaussians end
    unsigned int startLevel = atomIndex;
    unsigned int nextLevel = gv.gaussians.size();

    // Update level information
    gv.levels.push_back(nextLevel);

    // Loop overall possible levels of overlaps from 2 to 6
    for (unsigned int l = 2; l < LEVEL; ++l) {
	// List of atom-atom overlaps is made => gv.gaussians[startLevel .. nextLevel-1];
	// Now update the overlap lists for each overlap in this level
	// Create the next overlap Gaussian 
	// And add it to the vector of overlaps
	for (unsigned int i = startLevel; i < nextLevel; ++i) {
	    // Parent indices
	    unsigned int a1 = parents[i].first;
	    unsigned int a2 = parents[i].second;

	    // Append volume to end of overlap vector
	    // Add new empty set
	    std::set < unsigned int >*tmp = overlaps[i];
	    std::set_intersection(overlaps[a1]->begin(),
				  overlaps[a1]->end(),
				  overlaps[a2]->begin(),
				  overlaps[a2]->end(),
				  std::insert_iterator < std::set <
				  unsigned int > >(*tmp, tmp->begin()));

	    // Check if the overlap list is empty
	    if (overlaps[i]->empty()) {
		continue;
	    }
	    // Get the possible overlaps from the parent gaussians
	    // and create the new overlap volume
	    for (setIter = overlaps[i]->begin();
		 setIter != overlaps[i]->end(); ++setIter) {
		if (*setIter <= a2) {
		    continue;
		}
		// Create a new overlap gaussian
		AtomGaussian ga = atomIntersection(gv.gaussians[i],
						   gv.gaussians[*setIter]);

		// Check if the volume is large enough
		if (ga.volume /
		    (gv.gaussians[i].volume +
		     gv.gaussians[*setIter].volume - ga.volume) < EPS) {
		    continue;
		}

		gv.gaussians.push_back(ga);
		std::vector < unsigned int >*vec =
		    new std::vector < unsigned int >();
		gv.childOverlaps.push_back(vec);

		// Update local variables
		parents.push_back(std::make_pair < unsigned int,
				  unsigned int >(i, *setIter));
		std::set < unsigned int >*tmp =
		    new std::set < unsigned int >();
		overlaps.push_back(tmp);

		// Update volume, centroid and moments
		// Overlaps consisting of an even number of atoms have a negative contribution
		if ((ga.nbr % 2) == 0) {
		    // Update volume and centroid
		    gv.volume -= ga.volume;
		    gv.centroid.x -= ga.volume * ga.center.x;
		    gv.centroid.y -= ga.volume * ga.center.y;
		    gv.centroid.z -= ga.volume * ga.center.z;
		} else {
		    // Update volume and centroid
		    gv.volume += ga.volume;
		    gv.centroid.x += ga.volume * ga.center.x;
		    gv.centroid.y += ga.volume * ga.center.y;
		    gv.centroid.z += ga.volume * ga.center.z;
		}

		// Update child list of the first
		gv.childOverlaps[i]->push_back(vecIndex);

		// Move to next index in vector
		++vecIndex;
	    }
	}

	// Update levels
	startLevel = nextLevel;
	nextLevel = gv.gaussians.size();

	// Update level information
	gv.levels.push_back(nextLevel);
    }

    // cleanup current set of computed overlaps
    for (std::vector < std::set < unsigned int >*>::iterator si =
	 overlaps.begin(); si != overlaps.end(); ++si) {
	if (*si != NULL) {
	    delete *si;
	    *si = NULL;
	}
    }

    parents.clear();

    // Update self-overlap
    gv.overlap = atomOverlap(gv, gv);

    return;
}



void initOrientation(GaussianVolume & gv)
{
    double x(0.0), y(0.0), z(0.0);

    // Scale centroid and moments with self volume  
    gv.centroid.x /= gv.volume;
    gv.centroid.y /= gv.volume;
    gv.centroid.z /= gv.volume;

    //#define USE_EIGEN2
    // Compute moments of inertia from mass matrix
#ifdef USE_EIGEN2
    double sumXX=0,sumXY=0,sumXZ=0,sumYY=0,sumYZ=0,sumZZ=0;
#else
    RDNumeric::DoubleSymmMatrix mass(3,0.0);
#endif
    // Loop over all gaussians
    for (std::vector < AtomGaussian >::iterator i = gv.gaussians.begin();
	 i != gv.gaussians.end(); ++i) {
	// Translate to center
	i->center.x -= gv.centroid.x;
	i->center.y -= gv.centroid.y;
	i->center.z -= gv.centroid.z;

	x = i->center.x;
	y = i->center.y;
	z = i->center.z;

        int mult=(i->nbr % 2)?1:-1;
#ifdef USE_EIGEN2
        sumXX += mult*i->volume * x * x;
        sumXY += mult*i->volume * x * y;
        sumXZ += mult*i->volume * x * z;
        sumYY += mult*i->volume * y * y;
        sumYZ += mult*i->volume * y * z;
        sumZZ += mult*i->volume * z * z;
#else
        mass.setVal(0,0,mass.getVal(0,0) + mult*i->volume * x * x);
        mass.setVal(0,1,mass.getVal(0,1) + mult*i->volume * x * y);
        mass.setVal(0,2,mass.getVal(0,2) + mult*i->volume * x * z);
        mass.setVal(1,1,mass.getVal(1,1) + mult*i->volume * y * y);
        mass.setVal(1,2,mass.getVal(1,2) + mult*i->volume * y * z);
        mass.setVal(2,2,mass.getVal(2,2) + mult*i->volume * z * z);
#endif
    }
    // Normalize mass matrix
#ifdef USE_EIGEN2
    sumXX/=gv.volume;
    sumXY/=gv.volume;
    sumXZ/=gv.volume;
    sumYY/=gv.volume;
    sumYZ/=gv.volume;
    sumZZ/=gv.volume;
    Eigen::Matrix3d mat;
    mat << sumXX, sumXY, sumXZ,
      sumXY, sumYY, sumYZ,
      sumXZ, sumYZ, sumZZ;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
    if(eigensolver.info()!=Eigen::Success){
      BOOST_LOG(rdErrorLog)<<"eigenvalue calculation did not converge"<<std::endl;
      return;
    }
    //    std::cerr << " eigenvalues: " << eigensolver.eigenvalues() << std::endl;
    for(unsigned int i=0;i<3;++i){
      for(unsigned int j=0;j<3;++j){
        // the eigenvectors come back in increasing order, so reverse it here:
        gv.rotation[i][j] = eigensolver.eigenvectors()(i,2-j);
      }
    }
#else
    mass /= gv.volume;
    int nPts = gv.gaussians.size();
    RDGeom::Transform3D *tf=RDNumeric::computeCanonicalTransformFromCovMat(&mass,nPts);
    for(unsigned int i=0;i<3;++i){
      for(unsigned int j=0;j<3;++j){
        gv.rotation[i][j] = tf->getVal(j,i);
      }
    }
    delete tf;
#endif
#if 0
    std::cerr<<std::endl;
    for(unsigned int i=0;i<3;++i){
      for(unsigned int j=0;j<3;++j){
        std::cerr<<std::setprecision(5)<<gv.rotation[i][j]<<" ";
      }
      std::cerr<<std::endl;
    }
#endif    
    double det = gv.rotation[0][0] * gv.rotation[1][1] * gv.rotation[2][2]
	+ gv.rotation[2][1] * gv.rotation[1][0] * gv.rotation[0][2]
	+ gv.rotation[0][1] * gv.rotation[1][2] * gv.rotation[2][0]
	- gv.rotation[0][0] * gv.rotation[1][2] * gv.rotation[2][1]
	- gv.rotation[1][1] * gv.rotation[2][0] * gv.rotation[0][2]
	- gv.rotation[2][2] * gv.rotation[0][1] * gv.rotation[1][0];
    //std::cerr<<"  DET: "<<std::setprecision(3)<<det<<std::endl;
    // Check if it is a rotation matrix and not a mirroring
    if (det < 0) {
      //std::cerr<<"FLIP"<<std::endl;
#if 1
	// Switch sign of third column
	gv.rotation[0][2] = -gv.rotation[0][2];
	gv.rotation[1][2] = -gv.rotation[1][2];
	gv.rotation[2][2] = -gv.rotation[2][2];
#else
        for(unsigned int i=0;i<3;++i){
          for(unsigned int j=0;j<3;++j){
            gv.rotation[i][j]*=-1.;
          }
        }
#endif
    }

    // Rotate all gaussians
    for (std::vector < AtomGaussian >::iterator i = gv.gaussians.begin();
	 i != gv.gaussians.end(); ++i) {
	x = i->center.x;
	y = i->center.y;
	z = i->center.z;
	i->center.x =
	    gv.rotation[0][0] * x + gv.rotation[1][0] * y +
	    gv.rotation[2][0] * z;
	i->center.y =
	    gv.rotation[0][1] * x + gv.rotation[1][1] * y +
	    gv.rotation[2][1] * z;
	i->center.z =
	    gv.rotation[0][2] * x + gv.rotation[1][2] * y +
	    gv.rotation[2][2] * z;
    }

    return;
}



double atomOverlap(GaussianVolume & gRef, GaussianVolume & gDb)
{
    // Create a queue to hold the pairs to process
    std::queue < std::pair < unsigned int, unsigned int > >processQueue;

    // loop over the single atom volumes of both molecules and make the combinations
    unsigned int N1(gRef.levels[0]);
    unsigned int N2(gDb.levels[0]);

    double Cij(0.0), Vij(0.0);

    double dx(0.0), dy(0.0), dz(0.0);

    std::vector < unsigned int >*d1(NULL), *d2(NULL);
    std::vector < unsigned int >::iterator it1;

    // Overlap volume
    double overlapVol(0.0);

    // First compute atom-atom overlaps 
    for (unsigned int i(0); i < N1; ++i) {
	for (unsigned int j(0); j < N2; ++j) {
	    // Scaling constant
	    Cij =
		gRef.gaussians[i].alpha * gDb.gaussians[j].alpha /
		(gRef.gaussians[i].alpha + gDb.gaussians[j].alpha);

	    // Variables to store sum and difference of components
	    dx = (gRef.gaussians[i].center.x - gDb.gaussians[j].center.x);
	    dx *= dx;
	    dy = (gRef.gaussians[i].center.y - gDb.gaussians[j].center.y);
	    dy *= dy;
	    dz = (gRef.gaussians[i].center.z - gDb.gaussians[j].center.z);
	    dz *= dz;

	    // Compute overlap volume
	    Vij =
		gRef.gaussians[i].C * gDb.gaussians[j].C * pow(PI /
							       (gRef.gaussians
								[i].alpha +
								gDb.gaussians
								[j].alpha),
							       1.5) *
		exp(-Cij * (dx + dy + dz));

	    // Check if overlap is sufficient enough
	    if (Vij /
		(gRef.gaussians[i].volume + gDb.gaussians[j].volume -
		 Vij) < EPS) {
		continue;
	    }
	    // Add to overlap volume
	    overlapVol += Vij;

	    // Loop over child nodes and add to queue
	    d1 = gRef.childOverlaps[i];
	    d2 = gDb.childOverlaps[j];

	    // First add (i,child(j))
	    if (d2 != NULL) {
		for (it1 = d2->begin(); it1 != d2->end(); ++it1) {
		    processQueue.push(std::make_pair < unsigned int,
				      unsigned int >(i, *it1));
		}
	    }
	    // Second add (child(i,j))
	    if (d1 != NULL) {
		for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
		    // add (child(i),j)
		    processQueue.push(std::make_pair < unsigned int,
				      unsigned int >(*it1, j));
		}
	    }
	}
    }

    while (!processQueue.empty()) {
	// Get next element from queue
	std::pair < unsigned int, unsigned int >nextPair =
	    processQueue.front();
	processQueue.pop();

	unsigned int i = nextPair.first;
	unsigned int j = nextPair.second;

	// Scaling constant
	Cij =
	    gRef.gaussians[i].alpha * gDb.gaussians[j].alpha /
	    (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha);

	// Variables to store sum and difference of components
	dx = (gRef.gaussians[i].center.x - gDb.gaussians[j].center.x);
	dx *= dx;
	dy = (gRef.gaussians[i].center.y - gDb.gaussians[j].center.y);
	dy *= dy;
	dz = (gRef.gaussians[i].center.z - gDb.gaussians[j].center.z);
	dz *= dz;

	// Compute overlap volume
	Vij = gRef.gaussians[i].C * gDb.gaussians[j].C *
	    pow(PI / (gRef.gaussians[i].alpha + gDb.gaussians[j].alpha),
		1.5) * exp(-Cij * (dx + dy + dz));

	// Check if overlap is sufficient enough
	if (Vij /
	    (gRef.gaussians[i].volume + gDb.gaussians[j].volume - Vij) <
	    EPS) {
	    continue;
	}
	// Even number of overlap atoms => addition to volume
	// Odd number => substraction
	if ((gRef.gaussians[i].nbr + gDb.gaussians[j].nbr) % 2 == 0) {
	    overlapVol += Vij;
	} else {
	    overlapVol -= Vij;
	}

	// Loop over child nodes and add to queue
	d1 = gRef.childOverlaps[i];
	d2 = gDb.childOverlaps[j];
	if (d1 != NULL && gRef.gaussians[i].nbr > gDb.gaussians[j].nbr) {
	    for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
		// Add (child(i),j)
		processQueue.push(std::make_pair < unsigned int,
				  unsigned int >(*it1, j));
	    }
	} else {
	    // First add (i,child(j))
	    if (d2 != NULL) {
		for (it1 = d2->begin(); it1 != d2->end(); ++it1) {
		    processQueue.push(std::make_pair < unsigned int,
				      unsigned int >(i, *it1));
		}
	    }
	    if (d1 != NULL
		&& gDb.gaussians[j].nbr - gRef.gaussians[i].nbr < 2) {
		for (it1 = d1->begin(); it1 != d1->end(); ++it1) {
		    // add (child(i),j)
		    processQueue.push(std::make_pair < unsigned int,
				      unsigned int >(*it1, j));
		}
	    }
	}
    }

    return overlapVol;
}



double getScore(std::string & id, double Voa, double Vra, double Vda)
{
    // set the score by which molecules are being compared
    if (id == tanimoto) {
	return Voa / (Vra + Vda - Voa);
    } else if (id == tversky_ref) {
	return Voa / (0.95 * Vra + 0.05 * Vda);
    } else if (id == tversky_db) {
	return Voa / (0.05 * Vra + 0.95 * Vda);
    }

    return 0.0;
}




void
checkVolumes(GaussianVolume & gRef, GaussianVolume & gDb,
	     AlignmentInfo & res)
{
    if (res.overlap > gRef.overlap) {
	res.overlap = gRef.overlap;
    }
    if (res.overlap > gDb.overlap) {
	res.overlap = gDb.overlap;
    }
    return;
}
