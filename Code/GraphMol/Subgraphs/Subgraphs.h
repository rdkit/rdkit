//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
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

  INT_PATH_LIST_MAP findAllSubgraphsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
						 unsigned int upperLen, bool useHs=false,
						 bool verbose=false);

  PATH_LIST findAllSubgraphsOfLengthN(const ROMol &mol,unsigned int,
				      bool useHs=false,bool verbose=false);
  PATH_LIST findUniqueSubgraphsOfLengthN(const ROMol&,unsigned int,
					 bool useHs=false,bool useBO=true);
  PATH_LIST findAllPathsOfLengthN(const ROMol&,unsigned int,
				  bool useBonds=true,bool useHs=false);
}

  
#endif
