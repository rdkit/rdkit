/*******************************************************************************
bestResults.cpp - Shape-it
 
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



#include <Shape/bestResults.h>



BestResults::BestResults(unsigned int n)
{
   _bestList.clear();
   _size = n;
   _lowest = 0.0;
   _filled = 0;
}



BestResults::~BestResults(void)
{
	std::vector<SolutionInfo*>::iterator it;
	for (it = _bestList.begin(); it != _bestList.end(); ++it)
	{
		if (*it != NULL)
      {
			delete *it;
         *it = NULL;
		}
	}
}



bool
BestResults::add(SolutionInfo& res)
{
	std::vector<SolutionInfo* >::reverse_iterator it;
	if (_filled < _size)
	{
		SolutionInfo* i = new SolutionInfo(res);
		_bestList.push_back(i);
		++_filled;
	}
	else if (res.score < _lowest)
	{
		return false;
	}
	else
	{
		// delete last element
		it = _bestList.rbegin();
		if (*it != NULL)
      {
			delete *it;
         *it = NULL;
      }
		
		// make new info element in the list
		*it = new SolutionInfo(res);
	}
		
	std::sort(_bestList.begin(), _bestList.end(), BestResults::_compInfo());
	it = _bestList.rbegin();
	_lowest = (*it)->score;
	
	return true;
}


