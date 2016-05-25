/*
*  Copyright (C) 2016 Sereina Riniker, Paolo Tosco
* 
*    @@ All Rights Reserved @@
*   This file is part of the RDKit.
*   The contents are covered by the terms of the BSD license
*   which is included in the file license.txt, found at the root
*   of the RDKit source tree.
*/

%include <boost_shared_ptr.i>
%include <carrays.i>
%array_functions(double, SWIGArrayUtility)

%{
#include <GraphMol/Trajectory/Snapshot.h>
#include <GraphMol/Trajectory/Trajectory.h>
%}
%include <GraphMol/Trajectory/Snapshot.h>
%include <GraphMol/Trajectory/Trajectory.h>

%extend RDKit::Snapshot {
  Snapshot(double *pos, double energy = 0.0) :
    d_trajectory(NULL),
      d_energy(energy),
      d_pos(boost::shared_array<double>(pos)) {}
}
