// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <Shape/shape.h>
#include <Shape/solutionInfo.h>
#include <GraphMol/ROMol.h>
namespace python = boost::python;
  
namespace RDKit {
        float align(ROMol &refMol, ROMol &dbMol, 
                              bool scoreOnly=false, 
                              unsigned int nBestHits = 10,
                              unsigned int maxIterations=10,
                              std::string whichScore = "Tanimoto") { 
                              SolutionInfo si = shape(refMol, dbMol, nBestHits, scoreOnly, maxIterations, whichScore);
                              return si.score;
                              }
}

BOOST_PYTHON_MODULE(rdShape) {
  python::scope().attr("__doc__") =
    "Module containing functions to align molecules as groups of Gaussians. \n\n\
    The code was developed by Hans de Winter. \n\n\
    The original version of the code can be downloaded from : \n\n\
    http://www.silicos-it.com/software/shape-it/1.0.1/shape-it.html \n\n\
    "
    ;


  std::string docString = "Align dbMol to refMol by maximizing Gaussian volume overlap \n\n\
 \n\
 ARGUMENTS:\n\n\
    - refMol : reference molecule \n\
    - dbMol : database molecule \n\
    - nBestHits : store up to 2 best hits (default: 10)\n\
    - scoreOnly : only score, don't align (default: false)\n\
    - maxIteration : maximum number of iterations (default: 10) \n\
    - whichScore : possible options are Tanimoto, Tversky_Ref, Tversky_Db (default: Tanimoto) \n\
                  \n\
 RETURNS:\n\n\
    double Tanimoto volume overlap coefficient\n\
\n";
  python::def("Align", RDKit::align,
              (python::arg("refMol"), python::arg("dbMol"),
               python::arg("nBestHits")=0,
               python::arg("scoreOnly")=false, python::arg("maxIter")=0, python::arg("whichScore")="Tanimoto"),
              docString.c_str());
}
  