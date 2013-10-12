/*******************************************************************************
options.h - Shape-it
 
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



#ifndef __SILICOSIT_SHAPEIT_OPTIONS_H__
#define __SILICOSIT_SHAPEIT_OPTIONS_H__



// General
#include <string>
#include <iostream>
#include <fstream>

// OpenBabel

// RDKit
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

// Shape-it
const std::string tanimoto = "Shape-it::Tanimoto"; 
const std::string tversky_ref = "Shape-it::Tversky_Ref"; 
const std::string tversky_db = "Shape-it::Tversky_Db"; 



class Options
{
   public:
   
      std::string             refInpFile;       //  -r  --reference
      std::ifstream*          refInpStream;
      
      RDKit::SDMolSupplier* refInpReader;
      
      std::string             dbInpFile;        //  -d  --dbase
      std::ifstream*          dbInpStream;
      RDKit::SDMolSupplier*   dbInpReader;
      
      std::string             molOutFile;       //  -o  --out
      std::ofstream*          molOutStream;
      RDKit::SDWriter*        molOutWriter;
      
      std::string             scoreOutFile;     //  -s  --score
      std::ofstream*          scoreOutStream;
      
      unsigned int            bestHits;         //      --best
      double                  cutOff;           //      --cutOff
      double                  maxIter;          //      --maxIterations
  
      std::string             whichScore;       //      --rankBy
      
      bool                    scoreOnly;        //      --scoreOnly
      bool                    showRef;          //      --noRef
  
      bool                    version;          //  -v  --version
      bool                    help;             //  -h  --help
   
      Options(void);
      ~Options(void);
      
      std::string print(void) const;
};



#endif
