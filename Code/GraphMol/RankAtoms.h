//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
  template <typename T>
  struct pairGTFunctor {
    bool operator() (const std::pair<T,int> &v1,const std::pair<T,int> &v2){
      return v1.first > v2.first;
    }
  };

  //! function for implementing < on two std::pairs.  The first entries are compared.
  template <typename T>
  bool pairLess(const std::pair<T,int> &v1,const std::pair<T,int> &v2){
    return v1.first < v2.first;
  }

  //! ranks the entries in a vector
  /*!
    \param vect the vector to rank
    \param res  is used to return the ranks of each entry
  */
  template <typename T>  
  void rankVect(const std::vector<T> &vect,INT_VECT &res){
    PRECONDITION(res.size()>=vect.size(),"vector size mismatch");
    unsigned int nEntries = vect.size();

#if 0
    std::priority_queue< std::pair<T,int>,
                         std::vector< std::pair<T, int> >,
                         pairGTFunctor<T> > sortedVect;
    for(unsigned int i=0;i<nEntries;++i){
      sortedVect.push(std::make_pair(vect[i],i));
    }

    int currRank=0;
    T lastV = sortedVect.top().first;
    for(unsigned int i=0;i<nEntries;++i){
      const std::pair<T,int> &p = sortedVect.top();
      if(p.first==lastV){
        res[p.second] = currRank;
      } else {
        res[p.second] = ++currRank;
        lastV = p.first;
      }
      sortedVect.pop();
    }
#else
    std::vector< std::pair<T,int> > sortedVect;
    sortedVect.resize(nEntries);
    for(unsigned int i=0;i<nEntries;++i){
      sortedVect[i]=std::make_pair(vect[i],i);
    }
    std::sort(sortedVect.begin(),sortedVect.end(),pairLess<T>);
    int currRank=0;
    T lastV = sortedVect[0].first;
    for(unsigned int i=0;i<nEntries;++i){
      const std::pair<T,int> &p = sortedVect[i];
      if(p.first==lastV){
        res[p.second] = currRank;
      } else {
        res[p.second] = ++currRank;
        lastV = p.first;
      }
    }
#endif
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
    INT_VECT newRanks;
    newRanks.resize(vals.size());
    rankVect(vals,newRanks);

    // EFF: this double vector strategy is maybe not the speediest
    //   or most memory efficient.  The operation can probably be done
    //   in place (in ranks) pretty painlessly, but I want to get
    //   this damn thing WORKING.
    INT_VECT fixedRanks;
    fixedRanks.resize(nAtoms);

    // --------------
    //  
    // At the end of this operation, fixedRanks will contain the ranks
    // of atoms that are no longer active (-1 for active atoms).
    //
    // --------------
    fixedRanks = ranks;
    INT_LIST::const_iterator ilCIt;
    for(ilCIt=indicesInPlay.begin();ilCIt!=indicesInPlay.end();ilCIt++){
      fixedRanks[*ilCIt] = -1;
    }

#ifdef VERYVERBOSE_CANON
    std::cout << "new: ";
    debugVect(newRanks);
    std::cout << "fixed: ";
    debugVect(fixedRanks);
#endif

    INT_VECT idxVect;
    idxVect.reserve(indicesInPlay.size());
    std::copy(indicesInPlay.begin(),indicesInPlay.end(),idxVect.begin());

    // -------------
    //
    //  Loop over all the new ranks.  We'll know that we're done
    //  when currNewRank > maxNewRank
    //
    // -------------
    INT_VECT::iterator ivIt;
    int currNewRank= *(std::min_element(newRanks.begin(),newRanks.end()));
    int maxNewRank = *(std::max_element(newRanks.begin(),newRanks.end()));
    while(currNewRank<=maxNewRank){
      //
      // If this rank is already present in fixedRanks, increment
      //  this rank and all new ranks that are higher:
      //
      while(std::find(fixedRanks.begin(),fixedRanks.end(),currNewRank)!=fixedRanks.end()){
        for(ivIt=newRanks.begin();ivIt!=newRanks.end();ivIt++){
          if(*ivIt>=currNewRank)
            *ivIt += 1;
        }
        // increment both thie current rank *and* the maximum new rank
        currNewRank++;
        maxNewRank++;
      }

      //
      //  now grab all entries with this new rank and copy them into
      //  the fixedRanks list
      //
      ivIt=std::find(newRanks.begin(),newRanks.end(),currNewRank);
      while(ivIt!=newRanks.end()){
        int offset=ivIt-newRanks.begin();
        int idx = idxVect[offset];
        fixedRanks[idx] = currNewRank;
        ivIt++;
        ivIt=std::find(ivIt,newRanks.end(),currNewRank);
      }
      currNewRank++;
    }
    ranks = fixedRanks;
  }
}
#endif
