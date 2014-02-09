//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//! \file RankAtoms.h
/*!
    \brief Utility functionality used by atom rankers.
*/
#ifndef _RD_RANKATOMS_H_
#define _RD_RANKATOMS_H_

#include <queue>
#include <vector>
#include <list>
#include <algorithm>
#include <boost/foreach.hpp>

namespace RankAtoms {
  typedef std::vector<int> INT_VECT;
  typedef std::list<int> INT_LIST;

  //! utility function for ranking atoms
  void updateInPlayIndices(const INT_VECT &ranks,INT_LIST &indicesInPlay);

  //! returns the count of unique items in an std::vector
  template <typename T>
  unsigned int countClasses(const std::vector<T> &vect){
    std::vector<T> sortedVect = vect;
    std::sort(sortedVect.begin(),sortedVect.end(),std::less<T>());
    typename std::vector<T>::iterator newEnd=std::unique(sortedVect.begin(),sortedVect.end());
    return newEnd-sortedVect.begin();
  }

  //! functor for implementing > on two std::pairs.  The first entries are compared.
  template <typename T1, typename T2>
  struct pairGreater : public std::binary_function<std::pair<T1,T2>,std::pair<T1,T2>,bool> {
    bool operator() (const std::pair<T1,T2> &v1,const std::pair<T1,T2> &v2) const {
      return v1.first > v2.first;
    }
  };

  //! function for implementing < on two std::pairs.  The first entries are compared.
  template <typename T1, typename T2>
  struct pairLess : public std::binary_function<std::pair<T1,T2>,std::pair<T1,T2>,bool> {
    bool operator() (const std::pair<T1,T2> &v1,const std::pair<T1,T2> &v2) const {
      return v1.first < v2.first;
    }
  };

  template <typename T>
  class argless : public std::binary_function<T,T,bool> {
  public:
    argless(const T& c) : std::binary_function<T,T,bool>(), container(c) {};
    bool operator() (unsigned int v1,unsigned int v2) const {
      return container[v1]<container[v2];
    }
    const T &container;
  };

  
  //! ranks the entries in a vector
  /*!
    \param vect the vector to rank
    \param res  is used to return the ranks of each entry
  */
  template <typename T>  
  void rankVect(const std::vector<T> &vect,INT_VECT &res){
    PRECONDITION(res.size()>=vect.size(),"vector size mismatch");
    unsigned int nEntries = vect.size();

    std::vector< unsigned int > indices(nEntries);
    for(unsigned int i=0;i<nEntries;++i) indices[i]=i; 
    std::sort(indices.begin(),indices.end(),argless< std::vector<T> >(vect) );

    int currRank=0;
    T lastV = vect[indices[0]];
    BOOST_FOREACH(unsigned int idx,indices){
      T v = vect[idx];
      if(v==lastV){
        res[idx] = currRank;
      } else {
        res[idx] = ++currRank;
        lastV = v;
      }
    }
  }    

  //! finds the relative rankings of the entries in \c vals.
  /*!
    \param nAtoms         the number of points in play
    \param vals           the values to be ranked
    \param indicesInPlay  a list containing the indices that
           are being considered (only those entries in \c ranks
           that appear in \c indicesInPlay will be modified)
    \param ranks          the current ranks of entries, this is updated
           with new ranks
  */
  template <typename T>
  void sortAndRankVect(unsigned int nAtoms,
                       const std::vector<T> &vals,
                       const INT_LIST &indicesInPlay,
                       INT_VECT &ranks) {
    // --------------
    //
    // start by getting the internal ranking of the values passed in
    //
    // --------------
    INT_VECT newRanks(vals.size());
    rankVect(vals,newRanks);

    // --------------
    //  
    // At the end of this operation, ranks will contain the ranks
    // of atoms that are no longer active (-1 for active atoms).
    //
    // --------------
    BOOST_FOREACH(int idx,indicesInPlay){
      ranks[idx] = -1;
    }

    INT_VECT idxVect;
    idxVect.assign(indicesInPlay.begin(),indicesInPlay.end());

    // -------------
    //
    //  Loop over all the new ranks.  We'll know that we're done
    //  when currNewRank > maxNewRank
    //
    // -------------

    int currNewRank= *(std::min_element(newRanks.begin(),newRanks.end()));
    int maxNewRank = *(std::max_element(newRanks.begin(),newRanks.end()));
    while(currNewRank<=maxNewRank){
      //
      // If this rank is already present in ranks, increment
      //  this rank and all new ranks that are higher:
      //
      while(std::find(ranks.begin(),ranks.end(),currNewRank)!=ranks.end()){
        BOOST_FOREACH(int &rank,newRanks){
          if(rank>=currNewRank)
            ++rank;
        }
        // increment both the current rank *and* the maximum new rank
        ++currNewRank;
        ++maxNewRank;
      }

      //
      //  now grab all entries with this new rank and copy them into
      //  the ranks list
      //
      INT_VECT::iterator ivIt=std::find(newRanks.begin(),newRanks.end(),currNewRank);
      while(ivIt!=newRanks.end()){
        int offset=ivIt-newRanks.begin();
        int idx = idxVect[offset];
        ranks[idx] = currNewRank;
        ++ivIt;
        ivIt=std::find(ivIt,newRanks.end(),currNewRank);
      }
      ++currNewRank;
    }
  }
  template <typename T>
  void sortAndRankVect2(const std::vector<std::vector<T> > &vals,
                        const INT_LIST &indicesInPlay,
                        INT_VECT &ranks) {
    // --------------
    //
    // start by getting the internal ranking of the values passed in
    //
    // --------------
    INT_VECT newRanks(vals.size());
    rankVect(vals,newRanks);

    // --------------
    //  
    // At the end of this operation, ranks will contain the ranks
    // of atoms that are no longer active (-1 for active atoms).
    //
    // --------------
    BOOST_FOREACH(int idx,indicesInPlay){
      ranks[idx] = newRanks[idx];
    }
    return;

    INT_VECT idxVect;
    idxVect.assign(indicesInPlay.begin(),indicesInPlay.end());

    // -------------
    //
    //  Loop over all the new ranks.  We'll know that we're done
    //  when currNewRank > maxNewRank
    //
    // -------------

    int currNewRank= *(std::min_element(newRanks.begin(),newRanks.end()));
    int maxNewRank = *(std::max_element(newRanks.begin(),newRanks.end()));
    while(currNewRank<=maxNewRank){
      //
      // If this rank is already present in ranks, increment
      //  this rank and all new ranks that are higher:
      //
      while(std::find(ranks.begin(),ranks.end(),currNewRank)!=ranks.end()){
        BOOST_FOREACH(int &rank,newRanks){
          if(rank>=currNewRank)
            ++rank;
        }
        // increment both the current rank *and* the maximum new rank
        ++currNewRank;
        ++maxNewRank;
      }

      //
      //  now grab all entries with this new rank and copy them into
      //  the ranks list
      //
      INT_VECT::iterator ivIt=std::find(newRanks.begin(),newRanks.end(),currNewRank);
      while(ivIt!=newRanks.end()){
        int offset=ivIt-newRanks.begin();
        int idx = idxVect[offset];
        ranks[idx] = currNewRank;
        ++ivIt;
        ivIt=std::find(ivIt,newRanks.end(),currNewRank);
      }
      ++currNewRank;
    }
  }
}
#endif
