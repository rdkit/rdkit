

#include <Shape/shape.h>

SolutionInfo shape(RDKit::ROMol & refMol, RDKit::ROMol & dbMol,
		   unsigned int nBestHits, bool scoreOnly,
		   unsigned int maxIter, std::string whichScore)
{

    //std::string whichScore = "Shape-it::Tanimoto";

    // Create a list to store the best results
    BestResults *bestHits = NULL;
    if (bestHits != 0) {
	bestHits = new BestResults(nBestHits);
    }

    std::string refName = "";
    if (refMol.hasProp("Name")) {
	refMol.getProp("Name", refName);
    } else {
	refName = "Unnamed_ref";
    }

    // Create the refence set of Gaussians
    GaussianVolume refVolume;

    // List all Gaussians and their respective intersections
    listAtomVolumes(refMol, refVolume);

    // Move the Gaussian towards its center of geometry and align with principal axes
    if (!scoreOnly) {
	initOrientation(refVolume);
    }

    // Create a class to hold the best solution of an iteration
    SolutionInfo bestSolution;
    bestSolution.refName = refName;
    bestSolution.refAtomVolume = refVolume.overlap;
    bestSolution.refCenter = refVolume.centroid;
    bestSolution.refRotation = refVolume.rotation;



    // Open database stream
    unsigned molCount(0);

    std::ostringstream ss;


    ++molCount;

    // Keep track of the number of molecules processed so far
    if ((molCount % 10) == 0) {
	std::cerr << ".";
	if ((molCount % 500) == 0) {
	    std::cerr << " " << molCount << " molecules" << std::endl;
	}
    }

    std::string dbName = "";

    if (dbMol.hasProp("Name")) {
	dbMol.getProp("Name", dbName);
    }
    if (dbName == "") {
	ss.str("");
	ss << "MOL_" << molCount;
	dbName = ss.str();
	dbMol.setProp("Name", dbName);
    }
    bestSolution.dbName = dbName;

    // Create the set of Gaussians of database molecule
    GaussianVolume dbVolume;
    listAtomVolumes(dbMol, dbVolume);

    // Overlap with reference
    AlignmentInfo res;
    double bestScore(0.0);

    if (scoreOnly) {
	res.overlap = atomOverlap(refVolume, dbVolume);
	res.rotor[0] = 1.0;
	bestScore =
	    getScore(whichScore, res.overlap, refVolume.overlap,
		     dbVolume.overlap);
    } else {
	initOrientation(dbVolume);
	ShapeAlignment aligner(refVolume, dbVolume);
	aligner.setMaxIterations(maxIter);

	for (unsigned int l(0); l < 4; ++l) {
	    SiMath::Vector quat(4, 0.0);
	    quat[l] = 1.0;
	    AlignmentInfo nextRes = aligner.gradientAscent(quat);
	    checkVolumes(refVolume, dbVolume, nextRes);
	    double ss =
		getScore(whichScore, nextRes.overlap, refVolume.overlap,
			 dbVolume.overlap);
	    if (ss > bestScore) {
		res = nextRes;
		bestScore = ss;
	    }

	    if (bestScore > 0.98) {
		break;
	    }
	}

	// Check if additional simulated annealing steps are requested and start from the current best solution
	if (maxIter > 0) {
	    AlignmentInfo nextRes = aligner.simulatedAnnealing(res.rotor);
	    checkVolumes(refVolume, dbVolume, nextRes);
	    double ss =
		getScore(whichScore, nextRes.overlap, refVolume.overlap,
			 dbVolume.overlap);
	    if (ss > bestScore) {
		bestScore = ss;
		res = nextRes;
	    }
	}
    }



    // Cleanup local pointers to atom-gaussians 
    dbVolume.gaussians.clear();
    dbVolume.levels.clear();
    for (std::vector < std::vector < unsigned int >*>::iterator si =
	 dbVolume.childOverlaps.begin();
	 si != dbVolume.childOverlaps.end(); ++si) {
	if (*si != NULL) {
	    delete *si;
	    *si = NULL;
	}
    }
    dbVolume.childOverlaps.clear();

    // Optimal alignment information is stored in res and bestScore
    // => result reporting and post-processing
    updateSolutionInfo(bestSolution, res, bestScore, dbVolume);

    bestSolution.dbMol = &dbMol;
    bestSolution.dbName = dbName;

    // Post-process molecules 
    if (nBestHits) {

	// Add the score properties
	setAllScores(bestSolution);

	if (!scoreOnly) {
	    // Translate and rotate the molecule towards its centroid and inertia axes
	    positionMolecule(*bestSolution.dbMol, bestSolution.dbCenter,
			     bestSolution.dbRotation);

	    // Rotate molecule with the optimal
	    rotateMolecule(*bestSolution.dbMol, bestSolution.rotor);

	    // Rotate and translate the molecule with the inverse rotation and translation of the reference molecule
	    repositionMolecule(*bestSolution.dbMol, refVolume.rotation,
			       refVolume.centroid);

	}

	if (bestHits) {
	    bestHits->add(bestSolution);
	}
    }

    return bestSolution;
}
