/*******************************************************************************
parseCommandLine.cpp - Shape-it
 
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



#include <Shape/parseCommandLine.h>



Options
parseCommandLine(int argc, char* argv[])
{
	static struct option Arguments[] = {
		{ "version",           no_argument,         NULL,   'v' },
		{ "reference",         required_argument,   NULL,   'r' },
		{ "dbase",             required_argument,   NULL,   'd' },
		{ "scores",            required_argument,   NULL,   's' },
		{ "out",               required_argument,   NULL,   'o' },
		{ "format",            required_argument,	NULL,   'f' },
		{ "scoreOnly",         no_argument,         NULL,    1  },
		{ "rankBy",            required_argument,   NULL,    2  }, 
		{ "best",              required_argument,   NULL,    4  },
		{ "addIterations",     required_argument,   NULL,    5  },
		{ "cutoff",            required_argument,   NULL,    6  },
		{ "noRef",             no_argument,         NULL,   11  },
		{ "help",              no_argument,         NULL,   'h' }
	};
	
	Options o;
	
	int choice;
	opterr = 0;
	int optionIndex = 0;
	std::string s;
	std::string format;
	format.clear();
	
	while((choice = getopt_long(argc, argv,"vhpr:d:s:o:f:", Arguments, &optionIndex )) != -1)
	{
		switch (choice)
		{
			case 'v': //....................................................version 
			o.version = true;
            break;
            
			case 'r': //..................................................reference 
            o.refInpFile = optarg;
            o.refInpStream = new std::ifstream(optarg);
            if (!o.refInpStream->good()) { mainErr("Error opening input file for reference (-r)"); }
            o.refInpReader = new RDKit::SDMolSupplier(o.refInpStream);
            break;
            
			case 'd': //......................................................dbase
            o.dbInpFile = optarg;
            o.dbInpStream = new std::ifstream(optarg);
            if (!o.dbInpStream->good()) { mainErr("Error opening input file for database (-d)"); }
            o.dbInpReader = new RDKit::SDMolSupplier(o.dbInpStream);
            break;
            
			case 's': //.....................................................scores
            o.scoreOutFile = optarg;
            o.scoreOutStream = new std::ofstream(optarg);
            if (!o.scoreOutStream->good()) { mainErr("Error opening output file for scores (-s)"); }
            break;
            
			case 'o': //........................................................out
            o.molOutFile = optarg;
            o.molOutStream = new std::ofstream(optarg);
            if (!o.molOutStream->good()) { mainErr("Error opening output file for molecules (-o)"); }
            o.molOutWriter = new RDKit::SDWriter(o.molOutStream);
            break;

			case 'f': //.....................................................format
			format = optarg;
			break;
         
			case 1: //....................................................scoreOnly
			o.scoreOnly = true;
            break;
         
			case 2: //.......................................................rankBy
			s = optarg;
            transform(s.begin(), s.end(), s.begin(), toupper);
            if      (s == "TANIMOTO") { o.whichScore = tanimoto; }
            else if (s == "TVERSKY_DB") { o.whichScore = tversky_db; }
            else if (s == "TVERSKY_REF") { o.whichScore = tversky_ref; }
			break;
            
			case 4: //.........................................................best 
            o.bestHits = strtol(optarg, NULL, 10);
            break;
            
         	case 5: //................................................addIterations
            o.maxIter = strtol(optarg, NULL, 10);
			break;
            
         	case 6: //.......................................................cutoff
			o.cutOff = strtod(optarg, NULL);
            if      (o.cutOff > 1) { o.cutOff = 1.0; }
            else if (o.cutOff < 0) { o.cutOff = 0.0; }
            break;
         				
         	case 11: //.......................................................noRef
			o.showRef = false;
            break;
            
         	case 'h': //.......................................................help
            o.help = true;
            break;
				
         	default:
			mainErr("Unknown command line option");
		}
	}
	
	// If no options are given print the help
	if (optind == 1) { o.help = true; }
	
	argc -= optind;
	argv += optind;
	return o;
}
