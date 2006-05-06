
////////////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or              //
// modify it under the terms of the GNU General Public License                //
// version 2 as published by the Free Software Foundation.                    //
//                                                                            //
// This program is distributed in the hope that it will be useful, but        //
// WITHOUT ANY WARRANTY; without even the implied warranty of                 //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          //
// General Public License for more details.                                   //
//                                                                            //
// Written by François Fleuret                                                //
// Contact <francois.fleuret@epfl.ch> for comments & bug reports              //
// Copyright (C) 2004 EPFL                                                    //
////////////////////////////////////////////////////////////////////////////////

// $Id: misc.h,v 1.1 2005/03/03 15:52:35 fleuret Exp $

#ifndef MISC_H
#define MISC_H

template<class T> T sq(T x) { return x*x; }

struct Couple {
  int index;
  double value;
};

int compare_couple(const void *a, const void *b);

#endif
