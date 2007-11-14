//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

/*! \file Subgraphs.h

   \brief functionality for finding subgraphs and paths in molecules

   Difference between _subgraphs_ and _paths_ :
      Subgraphs are potentially branched, whereas paths (in our 
      terminology at least) cannot be.  So, the following graph:
 
            C--0--C--1--C--3--C
                  |
                  2
                  |
                  C
      has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
      but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
*/
#ifndef _RD_SUBGRAPHS_H_
#define _RD_SUBGRAPHS_H_

#include <vector>
#include <list>
#include <map>


namespace RDKit{
  class ROMol;
  // NOTE: before replacing the defn of PATH_TYPE: be aware that
  // we do occasionally use reverse iterators on these things, so
  // replacing with a slist would probably be a bad idea.
  typedef std::vector<int> PATH_TYPE;
  typedef std::list< PATH_TYPE > PATH_LIST;
  typedef PATH_LIST::const_iterator PATH_LIST_CI;

  typedef std::map<int, PATH_LIST> INT_PATH_LIST_MAP;
  typedef INT_PATH_LIST_MAP::const_iterator INT_PATH_LIST_MAP_CI;
  typedef INT_PATH_LIST_MAP::iterator INT_PATH_LIST_MAP_I;

  // --- --- --- --- --- --- --- --- --- --- --- --- --- 
  //
  //
  // --- --- --- --- --- --- --- --- --- --- --- --- --- 

  //! \brief find all subgraphs in a range of sizes
  /*!
   *   \param mol - the molecule to be considered
   *   \param lowerLen - the minimum subgraph size to find
   *   \param upperLen - the maximum subgraph size to find
   *   \param useHs     - if set, hydrogens in the graph will be considered
   *                      eligible to be in paths. NOTE: this will not add
   *                      Hs to the graph. 
  */
  INT_PATH_LIST_MAP findAllSubgraphsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
                                                  unsigned int upperLen, bool useHs=false);

  //! \brief find all subgraphs of a particular size
  /*!
   *   \param mol - the molecule to be considered
   *   \param targetLen - the length of the subgraphs to be returned
   *   \param useHs     - if set, hydrogens in the graph will be considered
   *                      eligible to be in paths. NOTE: this will not add
   *                      Hs to the graph. 
  */
  PATH_LIST findAllSubgraphsOfLengthN(const ROMol &mol,unsigned int targetLen,
                                      bool useHs=false);

  //! \brief find unique subgraphs of a particular size
  /*!
   *   \param mol - the molecule to be considered
   *   \param targetLen - the length of the subgraphs to be returned
   *   \param useHs     - if set, hydrogens in the graph will be considered
   *                      eligible to be in paths. NOTE: this will not add
   *                      Hs to the graph. 
   *   \param useBO     - if set, bond orders will be considered when uniquifying
   *                      the paths
  */
  PATH_LIST findUniqueSubgraphsOfLengthN(const ROMol &mol,unsigned int targetLen,
                                         bool useHs=false,bool useBO=true);
  //! \brief find all paths of a particular size
  /*!
   *   \param mol - the molecule to be considered
   *   \param targetLen - the length of the paths to be returned
   *   \param useBonds  - if set, the path indices will be bond indices,
   *                      not atom indices
   *   \param useHs     - if set, hydrogens in the graph will be considered
   *                      eligible to be in paths. NOTE: this will not add
   *                      Hs to the graph. 
  */
  PATH_LIST findAllPathsOfLengthN(const ROMol &mol,unsigned int targetLen,
                                  bool useBonds=true,bool useHs=false);
}

  
#endif
