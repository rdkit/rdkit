/*******************************************************************************
main.cpp - Shape-it
 
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



// General
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

// OpenBabel

// rdkit
#include <GraphMol/ROMol.h>
#include <RDGeneral/versions.h>

// Pharao
#include <Shape/options.h>
#include <Shape/gaussianVolume.h>
#include <Shape/atomGaussian.h>
#include <Shape/alignmentInfo.h>
#include <Shape/mainErr.h>
#include <Shape/bestResults.h>
#include <Shape/shapeAlignment.h>
#include <Shape/moleculeRotation.h>
#include <Shape/parseCommandLine.h>
#include <Shape/printHeader.h>
#include <Shape/printUsage.h>



//*--------------------------------------------------------------------------*//
//* MAIN                                                                MAIN *//
//*--------------------------------------------------------------------------*//
int main(int argc, char* argv[])
{  
	// Initialise random number generator
	srandom(time(NULL));
	clock_t t0 = clock();
   
   	// Print header
   	printHeader();
    std::cerr << "  -> RDKit: " << RDKit::rdkitVersion << std::endl;     
	// Read options
   	Options uo = parseCommandLine(argc,argv);
   	if (uo.version)
   	{
		printHeader();
	    std::cerr << "  -> RDKit: " << RDKit::rdkitVersion << std::endl;  
		exit(0);
	}
   
	if (uo.help)
   	{
		printUsage();
		exit(0);
	}
   	std::cerr << uo.print();
	
	// Files
	if (uo.dbInpFile.empty())
   	{
		mainErr("Missing database file. This is a required option (-d).");
	}
	if (uo.refInpFile.empty())
   	{
		mainErr("Missing ref file. This is a required option (-r).");
	}
   
	if (uo.molOutFile.empty() && uo.scoreOutFile.empty())
   	{
		mainErr("At least one of the -o or -s option should be used.");
	}	
	
	// Create a list to store the best results
	BestResults* bestHits = NULL;
	if (uo.bestHits != 0) { bestHits = new BestResults(uo.bestHits); }
   
	// Print header line to score output file
	if (!uo.scoreOutFile.empty())
   	{
		*(uo.scoreOutStream) << "dbName"
			<< "\t" << "refName"
			<< "\t" << tanimoto
			<< "\t" << tversky_ref
			<< "\t" << tversky_db
			<< "\t" << "overlap"
			<< "\t" << "refVolume"
			<< "\t" << "dbVolume"
			<< std::endl;
	}
	
	// Create reference molecule
	RDKit::ROMol refMol = *(uo.refInpReader)->next();
	std::string refName = "";
 	if (refMol.hasProp("Name"))  { refMol.getProp("Name", refName); } 
 	else { refName = "Unnamed_ref"; }
   
	// Create the refence set of Gaussians
	GaussianVolume refVolume;
   	
	// List all Gaussians and their respective intersections
	listAtomVolumes(refMol, refVolume);

	// Move the Gaussian towards its center of geometry and align with principal axes
	if ( !uo.scoreOnly ) { initOrientation(refVolume); }
	
	// Write reference molecule to output
	if (uo.showRef && !uo.molOutFile.empty())
   	{
      uo.molOutWriter->write(refMol);
	}
		
	// Create a class to hold the best solution of an iteration
	SolutionInfo bestSolution;
	bestSolution.refName = refName;
	bestSolution.refAtomVolume = refVolume.overlap;
	bestSolution.refCenter = refVolume.centroid;
	bestSolution.refRotation = refVolume.rotation;
	
   	// Open database stream
   	unsigned molCount(0);
 	RDKit::ROMol* dbMol_ptr;
   	std::ostringstream ss;
   	while (! uo.dbInpReader->atEnd() )
   	{
   	    dbMol_ptr = uo.dbInpReader->next();
   	    RDKit::ROMol dbMol = *dbMol_ptr;
      	++molCount;

      	// Keep track of the number of molecules processed so far
      	if ((molCount % 10) == 0)
      	{
         	std::cerr << ".";
         	if ((molCount % 500) == 0)
         	{
            	std::cerr << " " << molCount << " molecules" << std::endl;
         	}
      	}
      
      	std::string dbName = "";
      	
      	if (dbMol.hasProp("Name")) { dbMol.getProp("Name", dbName); }
      	if (dbName == "")
      	{
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
			
      	if (uo.scoreOnly)
      	{
         	res.overlap = atomOverlap(refVolume, dbVolume);
         	res.rotor[0] = 1.0;
         	bestScore = getScore(uo.whichScore, res.overlap, refVolume.overlap, dbVolume.overlap);
      	}
      	else
      	{
         	initOrientation(dbVolume);
         	ShapeAlignment aligner(refVolume, dbVolume);
         	aligner.setMaxIterations(uo.maxIter);
				
         	for (unsigned int l(0); l < 4; ++l)
         	{
            	SiMath::Vector quat(4,0.0);
            	quat[l] = 1.0;
            	AlignmentInfo nextRes = aligner.gradientAscent(quat);
            	checkVolumes(refVolume, dbVolume, nextRes);
            	double ss = getScore(uo.whichScore, nextRes.overlap, refVolume.overlap, dbVolume.overlap);
            	if (ss > bestScore)
            	{
               		res = nextRes;
               		bestScore = ss;
            	}
					
            	if (bestScore > 0.98)
            	{
               		break;
            	}
         	}
         
         	// Check if additional simulated annealing steps are requested and start from the current best solution
         	if (uo.maxIter > 0)
         	{
            	AlignmentInfo nextRes = aligner.simulatedAnnealing(res.rotor);
            	checkVolumes(refVolume, dbVolume, nextRes);
            	double ss = getScore(uo.whichScore, nextRes.overlap, refVolume.overlap, dbVolume.overlap);
            	if (ss > bestScore)
            	{
               		bestScore = ss;
               		res = nextRes;
            	}
         	}
      	}
			
      	// Cleanup local pointers to atom-gaussians 
      	dbVolume.gaussians.clear();
      	dbVolume.levels.clear();
      	for (std::vector<std::vector<unsigned int> *>::iterator si = dbVolume.childOverlaps.begin(); 
           si != dbVolume.childOverlaps.end();
           ++si)
      	{
         	if (*si != NULL)
         	{
            	delete *si;
				*si = NULL;
         	}
      	}
      	dbVolume.childOverlaps.clear();

      	// Optimal alignment information is stored in res and bestScore
      	// => result reporting and post-processing
      	updateSolutionInfo(bestSolution, res, bestScore, dbVolume);
        bestSolution.dbMol = dbMol_ptr;
        bestSolution.dbName = dbName;
        
      	// At this point the information of the solution is stored in bestSolution
      	// Check if the result is better than the cutoff
      	if (bestSolution.score < uo.cutOff)
      	{
         	continue;
      	}
			
      	// Post-process molecules 
      	if (uo.bestHits || !uo.molOutFile.empty())
      	{
         	// Add the score properties
         	setAllScores(bestSolution);

         	if (!uo.scoreOnly)
         	{

            	// Translate and rotate the molecule towards its centroid and inertia axes
            	positionMolecule(*bestSolution.dbMol, bestSolution.dbCenter, bestSolution.dbRotation);

            	// Rotate molecule with the optimal
            	rotateMolecule(*bestSolution.dbMol, bestSolution.rotor);

            	// Rotate and translate the molecule with the inverse rotation and translation of the reference molecule
            	repositionMolecule(*bestSolution.dbMol, refVolume.rotation, refVolume.centroid);

         	}
				
         	if (uo.bestHits)
         	{
            	bestHits->add(bestSolution);
         	}
         	else if (!uo.molOutFile.empty())
         	{
            	(uo.molOutWriter)->write(*bestSolution.dbMol);
            	(uo.molOutWriter)->flush();
         	}
      	}	
			
      	if ((uo.bestHits == 0) && !uo.scoreOutFile.empty())
      	{
			bestSolution.printScores(uo);
      	}
			
      	// reset best solution
      	bestSolution.score = 0.0;

	}

	if (uo.bestHits)
   	{
		if (!uo.molOutFile.empty())
      	{
			bestHits->writeMolecules(&uo);
			uo.molOutStream->flush();
			uo.molOutStream->close();
         	delete uo.molOutStream;
         	uo.molOutStream = NULL;
         	delete uo.molOutWriter;
         	uo.molOutWriter = NULL;
      	}
		if (!uo.scoreOutFile.empty())
      	{
			bestHits->writeScores(&uo);
         	delete uo.scoreOutStream;
         	uo.scoreOutStream = NULL;
      	}
	}

   	// Clear current streams
   	if (0) {
       	if (uo.dbInpStream != NULL)
       	{
          	delete uo.dbInpReader;
          	uo.dbInpReader = NULL;
          	delete uo.dbInpStream;
          	uo.dbInpStream = NULL;
       	}
       	if (uo.refInpStream != NULL)
       	{
          	delete uo.refInpReader;
          	uo.refInpReader = NULL;
          	delete uo.refInpStream;
          	uo.refInpStream = NULL;
       	}
   	}
	
	// Done processing database
   	std::cerr << std::endl;
   	std::cerr << "Processed " << molCount << " molecules" << std::endl;
   	double tt = (double)(clock() - t0 )/CLOCKS_PER_SEC;
   	std::cerr << molCount << " molecules in " << tt << " seconds (";
   	std::cerr << molCount/tt << " molecules per second)" << std::endl;
	
	// Cleanup local db volume
	refVolume.gaussians.clear();
	for (std::vector<std::vector<unsigned int> *>::iterator si = refVolume.childOverlaps.begin();
        si != refVolume.childOverlaps.end();
        ++si)
   	{
		if (*si != NULL)
      	{
			delete *si;
         	*si = NULL;
		}
	}
    
	exit(0);
}

