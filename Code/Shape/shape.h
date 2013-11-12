// General
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>

// OpenBabel

// rdkit
#include <GraphMol/ROMol.h>
//#include <RDGeneral/versions.h>

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

SolutionInfo shape(RDKit::ROMol& refMol, RDKit::ROMol& dbMol, unsigned int nBestHits, bool scoreOnly, unsigned int maxIter, std::string whichScore);