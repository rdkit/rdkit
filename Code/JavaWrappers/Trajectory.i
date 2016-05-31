/*
*  Copyright (C) 2016 Sereina Riniker, Paolo Tosco
* 
*    @@ All Rights Reserved @@
*   This file is part of the RDKit.
*   The contents are covered by the terms of the BSD license
*   which is included in the file license.txt, found at the root
*   of the RDKit source tree.
*/

%{
#include <GraphMol/Trajectory/Snapshot.h>
#include <GraphMol/Trajectory/Trajectory.h>
#include <boost/shared_array.hpp>
%}
%include <GraphMol/Trajectory/Snapshot.h>
%include <GraphMol/Trajectory/Trajectory.h>
%include <boost/shared_array.hpp>

%extend RDKit::Snapshot {
  Snapshot(double *posArray, double energy = 0.0) {
    std::cerr << "ok1, posArray = " << posArray << std::endl;
    boost::shared_array<double> pos(posArray);
    std::cerr << "ok2, posArray = " << pos.get() << std::endl;
    return new RDKit::Snapshot(pos, energy);
  }
}
