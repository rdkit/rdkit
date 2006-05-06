
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

// $Id: misc.cc,v 1.1 2005/03/03 15:52:35 fleuret Exp $

#include "misc.h"

int compare_couple(const void *a, const void *b) {
  if(((Couple *) a)->value < ((Couple *) b)->value) return -1;
  else if(((Couple *) a)->value > ((Couple *) b)->value) return 1;
  else return 0;
}

