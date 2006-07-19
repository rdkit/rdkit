//
//  Copyright (C) 2006 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_CHEMTRANSFORMS_H__
#define _RD_CHEMTRANSFORMS_H__

#include <boost/smart_ptr.hpp>
#include <vector>

namespace RDKit{
  class ROMol;
  typedef boost::shared_ptr<ROMol>    ROMOL_SPTR;

  //! \brief Returns a copy of an ROMol with the atoms and bonds that 
  //!      match a pattern removed.
  /*!
      \param mol       the ROMol of interest
      \param query     the query ROMol
      \param onlyFrags if this is set, only matches that correspond to an
                       entire fragment will be removed.  This is useful for
		       salt removal

      \return a copy of \c mol with the matching atoms and bonds (if any)
              removed.		       
  */
  ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query,
			  bool replaceAll=false);

  //! \brief Returns a list of copies of an ROMol with the atoms and bonds that 
  //!      match a pattern replaced with the atoms contained in another molecule.
  /*!
     Bonds are created between the joining atom in the existing molecule
     and the atoms in the new molecule. So, using SMILES instead of molecules:
            replaceSubstructs('OC(=O)NCCNC(=O)O','C(=O)O','[X]') ->
	          ['[X]NCCNC(=O)O','OC(=O)NCCN[X]']
            replaceSubstructs('OC(=O)NCCNC(=O)O','C(=O)O','[X]',true) ->
	          ['[X]NCCN[X]']
     Chains should be handled "correctly":
            replaceSubstructs('CC(=O)C','C(=O)','[X]') ->
	          ['C[X]C']
     As should rings:
            replaceSubstructs('C1C(=O)C1','C(=O)','[X]') ->
	          ['C1[X]C1']
     And higher order branches:
            replaceSubstructs('CC(=O)(C)C','C(=O)','[X]') ->
	          ['C[X](C)C']
     Note that the client is responsible for making sure that the
       resulting molecule actually makes sense - this function does not
       perform sanitization.

      \param mol         the ROMol of interest
      \param query       the query ROMol
      \param replacement the ROMol to be inserted
      \param replaceAll  if this is true, only a single result, with all occurances
                         of the substructure replaced, will be returned.

      \return a vector of pointers to copies of \c mol with the matching atoms
          and bonds (if any) replaced

  */
  std::vector<ROMOL_SPTR> replaceSubstructs(const ROMol &mol, const ROMol &query,
					    const ROMol &replacement,
					    bool replaceAll=false);
}

#endif




