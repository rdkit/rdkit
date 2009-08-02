// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <boost/python.hpp>
#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>

namespace python = boost::python;
using namespace RDKit;


BOOST_PYTHON_MODULE(rdsubgraphs)
{

  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for subgraphs."
    ;
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  std::string docString;

  // ------------------------------------------------------------------------
  docString="Finds all subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - verbose: (optional, internal use) toggles verbosity in the search algorithm.\n\
      Defaults to 0.\n\
\n\
  RETURNS: a tuple of 2-tuples with bond IDs\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
  has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
  but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
  python::def("FindAllSubgraphsOfLengthN", &findAllSubgraphsOfLengthN,
              (python::arg("mol"),python::arg("length"),
               python::arg("useHs")=false),
              docString.c_str());
  // ------------------------------------------------------------------------
  docString="Finds unique subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from\n\
      another.\n\
      Defaults to 1.\n\
\n\
  RETURNS: a tuple of tuples with bond IDs\n\
\n\
\n";
  python::def("FindUniqueSubgraphsOfLengthN", &findUniqueSubgraphsOfLengthN, 
              (python::arg("mol"),python::arg("length"),
               python::arg("useHs")=false,python::arg("useBO")=true),
              docString.c_str());
                  
  // ------------------------------------------------------------------------
  docString="Finds all paths of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target length for the paths.\n\
\n\
    - useBonds: (optional) toggles the use of bond indices in the paths.\n\
      Otherwise atom indices are used.  *Note* this behavior is different\n\
      from that for subgraphs.\n\
      Defaults to 1.\n\
\n\
  RETURNS: a tuple of tuples with IDs for the bonds.\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
\n\
       has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
       but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
  python::def("FindAllPathsOfLengthN", &findAllPathsOfLengthN, 
              (python::arg("mol"),python::arg("length"),
               python::arg("useBonds")=true,python::arg("useHs")=false),
              docString.c_str());




}


