//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBGRAPHUTILS_H_
#define _RD_SUBGRAPHUTILS_H_

#include "Subgraphs.h"
#include <boost/tuple/tuple.hpp>

namespace RDKit{
  class ROMol;
  
  //typedef DiscrimTuple PathDiscrimTuple;
  typedef boost::tuples::tuple<double,double,double> PathDiscrimTuple;
  PathDiscrimTuple CalcPathDiscriminators(const ROMol &mol,const PATH_TYPE &path,
					  bool useBO=true);
  PATH_LIST uniquifyPaths (const ROMol &mol, const PATH_LIST &allPathsb,
			   bool useBO=true,double tol=1e-8);

  // Return the list of bond that connect a list of atoms
  // ASSUMPTION: the atoms specified in the list are connected
  PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE atomIds);

  // create a new molecule object from a part of molecule "mol". The part of
  // of the molecule is specified as a list of bonds in "path".
  // the optional argument "useQuery" will set all the bond and atoms in the 
  // the new molecule to "QueryAtoms" and "QueryBonds" instead of regular Atoms and Bonds
  //  atomIdxMap provides a mapping between the atomsIds in mol to the atomIds in
  // the newly created sub-molecule (the molecule that is returned)
  ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
		      bool useQuery,
		      std::map<int,int> &atomIdxMap);

  // same as PathToSubmol above except for a bunch of optional arguments
  ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
		      bool useQuery=false);

		      
}


#endif
