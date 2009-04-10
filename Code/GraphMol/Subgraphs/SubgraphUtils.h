//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBGRAPHUTILS_H_
#define _RD_SUBGRAPHUTILS_H_

#include "Subgraphs.h"
#include <boost/tuple/tuple.hpp>

namespace RDKit{
  class ROMol;
  namespace Subgraphs {
    //! used to return atomic discriminators (three doubles)
    typedef boost::tuples::tuple<double,double,double> DiscrimTuple;
    //! calculates a set of molecular discriminators from the distance matrix
    /*!
      Computes:
        -# BalabanJ 
        -# the first eigenvalue of the distance matrix 
        -# the last but one eigenvalue of the distance matrix 

      \param mol    the molecule of interest
      \param useBO  toggles inclusion of the bond order in the discriminators
                    (when false, the discriminators are purely topological)
      \param force  forces the calculation (instead of using cached results)
	
      \return a \c DiscrimTuple with the results
      
    */
    DiscrimTuple computeDiscriminators(const ROMol &mol, 
					      bool useBO=true, 
					      bool force=false);
    //! \brief Same as MolOps::computeDiscriminators(const ROMol &mol),
    //! except that this directly uses the user-supplied distance matrix
    DiscrimTuple computeDiscriminators(double *distMat, unsigned int nb, unsigned int na);

    DiscrimTuple CalcPathDiscriminators(const ROMol &mol,const PATH_TYPE &path,
                                           bool useBO=true);
    PATH_LIST uniquifyPaths (const ROMol &mol, const PATH_LIST &allPathsb,
                             bool useBO=true,double tol=1e-8);

    // Return the list of bond that connect a list of atoms
    // ASSUMPTION: the atoms specified in the list are connected
    PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE &atomIds);

    // create a new molecule object from a part of molecule "mol". The part of
    // of the molecule is specified as a list of bonds in "path".
    // the optional argument "useQuery" will set all the bond and atoms in the 
    // the new molecule to "QueryAtoms" and "QueryBonds" instead of regular Atoms and Bonds
    //  atomIdxMap provides a mapping between the atomsIds in mol to the atomIds in
    // the newly created sub-molecule (the molecule that is returned)
    ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
                        bool useQuery,
                        std::map<int,int> &atomIdxMap);
    ROMol *PathToSubmol(const ROMol &mol, const PATH_TYPE &path, 
                        bool useQuery=false);
  } // end of namespace Subgraphs
}


#endif
