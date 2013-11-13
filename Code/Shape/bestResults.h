/*******************************************************************************
bestResults.h - Shape-it
 
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



#ifndef __SILICOSIT_SHAPEIT_BESTRESULTS_H__
#define __SILICOSIT_SHAPEIT_BESTRESULTS_H__




// General
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

// OpenBabel

// Shape-it
#include <Shape/solutionInfo.h>



class BestResults {
  private:

    std::vector < SolutionInfo * >_bestList;	///< Local list to best N solutions

    double _lowest;		///< lowest score in the list
    unsigned int _size;		///< total number of elements to be stored in the list
    unsigned int _filled;	///< number of elements stored in the list sofar

    class _compInfo {
      public:

	bool operator() (const SolutionInfo * a, const SolutionInfo * b) {
	    return a->score > b->score;
	};
    };

  public:

    BestResults(unsigned int n = 100);
    ~BestResults(void);

    bool add(SolutionInfo & res);
};



#endif
