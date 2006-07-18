//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SUBSTRUCTMATCH_H__
#define _RD_SUBSTRUCTMATCH_H__

// std bits
#include <vector>

// vflib bits
#include <argraph.h>
#include <vf_mono_state.h>
#include <vf2_mono_state.h>

namespace RDKit{
  class ROMol;
  class Atom;
  class Bond;

  //! \brief used to return matches from substructure searching,
  //!   The format is (queryAtomIdx, molAtomIdx)
  typedef std::vector< std::pair<int,int> > MatchVectType; 
  //typedef VFMonoState MatcherState;
  typedef VF2MonoState MatcherState;

  typedef ARGraph<const Atom, const Bond> AR_MOLGRAPH;
  //! Internal Use Only
  AR_MOLGRAPH *getMolGraph(const ROMol &mol);

  //! Find a substructure match for a query in a molecule
  /*!
      \param mol       The ROMol to be searched
      \param query     The query ROMol
      \param matchVect Used to return the match
                       (pre-existing contents will be deleted)
      \param recursionPossible  flags whether or not recursive matches are allowed

      \return whether or not a match was found
    
  */
  bool SubstructMatch(const ROMol &mol,const ROMol &query,
		      MatchVectType &matchVect,
		      bool recursionPossible=true);
  //! \overload
  bool SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,
		      MatchVectType &matchVect,
		      bool recursionPossible=true);

  //! Find all substructure matches for a query in a molecule
  /*!
      \param mol       The ROMol to be searched
      \param query     The query ROMol
      \param matchVect Used to return the matches
                       (pre-existing contents will be deleted)
      \param uniquify  Toggles uniquification (by atom index) of the results
      \param recursionPossible  flags whether or not recursive matches are allowed

      \return the number of matches found
    
  */
  unsigned int SubstructMatch(const ROMol &mol,const ROMol &query,
		     std::vector< MatchVectType > &matchVect,
		     bool uniquify=true,bool recursionPossible=true);
  //! \overload
  unsigned int SubstructMatch(AR_MOLGRAPH *molG,const ROMol &query,
		     std::vector< MatchVectType > &matchVect,
		     bool uniquify=true,bool recursionPossible=true);
}

#endif




