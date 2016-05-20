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
%array_functions(double, double_array)
%shared_ptr(DoubleArray)

%inline %{
#include <boost/shared_array.hpp>

struct DoubleArray {
  double *value;
  DoubleArray(double *v) : value(v) {}
};
%}

%{
#include <GraphMol/Trajectory/Snapshot.h>
#include <GraphMol/Trajectory/Trajectory.h>
%}
%include <GraphMol/Trajectory/Snapshot.h>
%include <GraphMol/Trajectory/Trajectory.h>
/*
http://www.swig.org/Doc3.0/Library.html#Library_std_shared_ptr
http://www.swig.org/Doc1.3/Library.html#Library_carrays
*/
