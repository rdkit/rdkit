/*******************************************************************************
printHeader.cpp - Shape-it
 
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



#include <Shape/printHeader.h>



void
printHeader(void)
{
	std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cerr << "  Shape-it v" << SHAPEIT_VERSION << "." << SHAPEIT_RELEASE << "." << SHAPEIT_SUBRELEASE << " | ";
    std::cerr << __DATE__ " " << __TIME__ << std::endl;
	std::cerr << std::endl;
	std::cerr << "  -> GCC:       " << __VERSION__ << std::endl;
	std::cerr << "  -> OpenBabel: " << BABEL_VERSION << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Copyright 2012 by Silicos-it, a division of Imacosi BVBA" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Shape-it is free software: you can redistribute it and/or modify" << std::endl;
	std::cerr << "  it under the terms of the GNU Lesser General Public License as published" << std::endl;
	std::cerr << "  by the Free Software Foundation, either version 3 of the License, or" << std::endl;
	std::cerr << "  (at your option) any later version." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Shape-it is distributed in the hope that it will be useful," << std::endl;
	std::cerr << "  but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
	std::cerr << "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
	std::cerr << "  GNU Lesser General Public License for more details." << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Shape-it is linked against OpenBabel version 2." << std::endl;
	std::cerr << "  OpenBabel is free software; you can redistribute it and/or modify" << std::endl;
	std::cerr << "  it under the terms of the GNU General Public License as published by" << std::endl;
	std::cerr << "  the Free Software Foundation version 2 of the License." << std::endl;
	std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cerr << std::endl;
   return;
}
