
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

// $Id: create_samples.cc,v 1.1 2005/03/03 15:52:35 fleuret Exp $

using namespace std;

#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>

// This defines the positive population

bool in_region(double x, double y) {
  return x*x + y*y <= 0.25;
}

int main(int argc, char **argv) {

  const int nb_samples = 1000;

  const int nb_features = 1000;
  double vxf[nb_features], vyf[nb_features], kf[nb_features];
  for(int j = 0; j < nb_features; j++) {
    double alpha = drand48() * 2 * M_PI;
    vxf[j] = sin(alpha);
    vyf[j] = cos(alpha);
    kf[j]  = - drand48() * vxf[j] - drand48() * vyf[j];
  }

  cout << "Saving the training set.\n";
  ofstream training("train.dat");
  training << nb_samples << " " << nb_features << "\n";
  for(int k = 0; k < nb_samples; k++) {
    double x = drand48(), y = drand48();
    for(int j = 0; j < nb_features; j++)
      training << ((x*vxf[j] + y*vyf[j] + kf[j] >= 0) ? 1 : 0) << ((j < nb_features-1) ? " " : "\n");
    training << (in_region(x, y) ? 1 : 0) << "\n";
  }

  int delta = 50;
  cout << "Saving the test set.\n";
  ofstream test("test.dat");
  test << delta*delta << " " << nb_features << "\n";
  for(int xx = 0; xx < delta; xx++) for(int yy = 0; yy < delta; yy++) {
    double x = double(xx)/double(delta-1), y = double(yy)/double(delta-1);
    for(int j = 0; j < nb_features; j++)
      test << ((x*vxf[j] + y*vyf[j] + kf[j] >= 0) ? 1 : 0) << ((j < nb_features-1) ? " " : "\n");
    test << (in_region(x, y) ? 1 : 0) << "\n";
  }
}
