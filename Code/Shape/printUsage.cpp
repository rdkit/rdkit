/*******************************************************************************
printUsage.cpp - Shape-it
 
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



#include <Shape/printUsage.h>


void printUsage(void)
{
	std::cerr << std::endl;
	std::cerr << "TASK:" << std::endl << std::endl;
	std::cerr << "  Shape-it is a tool to align pairs of molecules based on their maximal" << std::endl;
	std::cerr << "  volume overlap." << std::endl;
	std::cerr << std::endl;
	std::cerr << std::endl;
	std::cerr << "REQUIRED OPTIONS: " << std::endl;
	std::cerr << "  -r, --reference <file>" << std::endl;
	std::cerr << "           File of the reference molecule with 3D coordinates." << std::endl;
	std::cerr << "           Only the first molecule in the reference file will be used." << std::endl;
	std::cerr << "           Shape-it can also handle a gzipped files if the extension is '.gz'" << std::endl;
	std::cerr << "           All input formats which are recognized by OpenBabel are allowed." << std::endl;
	std::cerr << "  -d, --dbase <file>" << std::endl;
	std::cerr << "           File of the database molecules with 3D coordinates." << std::endl;
	std::cerr << "           Shape-it can also handle gzipped files if the extension is '.gz'" << std::endl;
	std::cerr << "           All input formats which are recognized by OpenBabel are allowed." << std::endl;
	std::cerr << std::endl;
	std::cerr << std::endl;
	std::cerr << "OUTPUT OPTIONS: " << std::endl;
	std::cerr << "One of these two output options is required:" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -o, --out <file>" << std::endl;
	std::cerr << "           File to write all database or the N best molecules such that their" << std::endl;
	std::cerr << "           coordinates correspond to the best alignment with the reference molecule." << std::endl;
	std::cerr << "           The first molecule in the file is the reference molecule. When this file" << std::endl;
	std::cerr << "           if of type 'sdf', then each molecule contains a set of properties in which" << std::endl;
	std::cerr << "           the respective scores are reported. These fields are labeled with an" << std::endl;
	std::cerr << "           identifier starting with the tag Shape-it::" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -s, --scores <file>" << std::endl;
	std::cerr << "           Tab-delimited text file with the scores of molecules." << std::endl;
	std::cerr << "           When the N best scoring molecules are reported the molecules are ranked" << std::endl;
	std::cerr << "           with the descending scores." << std::endl;
	std::cerr << std::endl;
	std::cerr << std::endl;
	std::cerr << "OPTIONAL OPTIONS: " << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -f, --format <format>" << std::endl;
	std::cerr << "          Specifies the format of the reference, database and output files. If not" << std::endl;
	std::cerr << "          provided, then the formats are determined from the respective file extensions." << std::endl;
	std::cerr << "          The specified format string should be one of the formats recognised" << std::endl;
	std::cerr << "          by OpenBabel." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --best <N> " << std::endl;
	std::cerr << "           When this option is used, only the N best scoring alignments will be" <<  std::endl;
	std::cerr << "           reported. The scoring function is defined by the --rankBy option." << std::endl;
	std::cerr << "           By default all molecules in the database are reported with their" << std::endl;
	std::cerr << "           respective scores without any ordering." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --scoreOnly" << std::endl;
	std::cerr << "           When this option is used the molecules are not aligned, only the volume" << std::endl;
	std::cerr << "           overlap between the reference and the given pose is computed." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --addIterations <nbr>" << std::endl;
	std::cerr << "           Sets the number of additional iterations in the simulated annealing" << std::endl;
	std::cerr << "           optimization step. The default value is set to 0, which refers to only" << std::endl;
	std::cerr << "           a local gradient ascent. Increasing the number of iterations will add" << std::endl;
	std::cerr << "           additional steps, and might give better alignments but it also takes" << std::endl;
	std::cerr << "           more time." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --rankBy <code>" << std::endl;
	std::cerr << "           This option can be used in combination with --best of --cutoff to rank" << std::endl;
	std::cerr << "           the molecules according to a given scoring function. The type of scoring" << std::endl;
	std::cerr << "           function is indicated with a code:" << std::endl;
	std::cerr << "             - TANIMOTO = Taninoto" << std::endl;
	std::cerr << "             - TVERSKY_REF = reference Tversky" << std::endl;
	std::cerr << "             - TVERSKY_DB = database Tversky" << std::endl;
	std::cerr << "           By default TANIMOTO is used." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --cutoff <value>" << std::endl;
	std::cerr << "           Defines a cutoff value. Only molecules with a score higher than the" << std::endl;
	std::cerr << "           cutoff are reported in the results files. Default value is set to" << std::endl;
	std::cerr << "           0.0. The scoring function is defined by the --rankBy option." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  --noRef" << std::endl;
	std::cerr << "           By default the reference molecule is written in the output files." << std::endl;
	std::cerr << "           Use this option to switch off this behavior." << std::endl;
	std::cerr << std::endl;
	std::cerr << std::endl;
	std::cerr << "HELP: " << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -h, --help" << std::endl;
	std::cerr << "           Prints this help overview." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -v, --version" << std::endl;
	std::cerr << "           Prints the version of the program." << std::endl;	
	std::cerr << std::endl;
   	return;
}
