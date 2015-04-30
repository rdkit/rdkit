//
//  Copyright (C) 2006-2012 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_CHEMTRANSFORMS_H__
#define _RD_CHEMTRANSFORMS_H__

#include <boost/smart_ptr.hpp>
#include <vector>
#include <iostream>

#include "MolFragmenter.h"

namespace RDKit{
  class ROMol;
  typedef boost::shared_ptr<ROMol>    ROMOL_SPTR;

  //! \brief Returns a copy of an ROMol with the atoms and bonds that 
  //!      match a pattern removed.
  /*!
      \param mol       the ROMol of interest
      \param query     the query ROMol
      \param onlyFrags  if this is set, atoms will only be removed if
                        the entire fragment in which they are found is
                        matched by the query.

      \return a copy of \c mol with the matching atoms and bonds (if any)
              removed.		       
  */
  ROMol *deleteSubstructs(const ROMol &mol, const ROMol &query,
			  bool onlyFrags=false);

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
      \param replacementConnectionPoint   index of the atom in the replacement that
                                          the bond should made to

      \return a vector of pointers to copies of \c mol with the matching atoms
          and bonds (if any) replaced

  */
  std::vector<ROMOL_SPTR> replaceSubstructs(const ROMol &mol, const ROMol &query,
					    const ROMol &replacement,
					    bool replaceAll=false,
                                            unsigned int replacementConnectionPoint=0);

  //! \brief Returns a copy of an ROMol with the atoms and bonds that 
  //!      don't fall within a substructure match removed.
  //!
  //!   dummy atoms are left to indicate attachment points.
  //!
  /*!
      \param mol       the ROMol of interest
      \param coreQuery a query ROMol to be used to match the core

      \return a copy of \c mol with the non-matching atoms and bonds (if any)
              removed and dummies at the connection points.		       
  */
  ROMol *replaceSidechains(const ROMol &mol, const ROMol &coreQuery);
			  
  //! \brief Returns a copy of an ROMol with the atoms and bonds that 
  //!      do fall within a substructure match removed.
  //!
  //!   dummy atoms are left to indicate attachment points.
  //!
  /*!
      Note that this is essentially identical to the replaceSidechains function, except we
      invert the query and replace the atoms that *do* match the query.

      \param mol            - the ROMol of interest
      \param coreQuery      - a query ROMol to be used to match the core
      \param replaceDummies - if set, atoms matching dummies in the core will also be replaced
      \param labelByIndex  - if set, the dummy atoms at attachment points are labelled with the
                             index+1 of the corresponding atom in the core
      \param requireDummyMatch - if set, only side chains that are connected to atoms in
                                 the core that have attached dummies will be considered.
                                 Molecules that have sidechains that are attached
                                 at other points will be rejected (NULL returned).

      \return a copy of \c mol with the non-matching atoms and bonds (if any)
              removed and dummies at the connection points. The client is responsible
              for deleting this molecule. If the core query is not matched, NULL is returned.
  */
  ROMol *replaceCore(const ROMol &mol, const ROMol &coreQuery,
                     bool replaceDummies=true,bool labelByIndex=false,
                     bool requireDummyMatch=false);

  //! \brief Carries out a Murcko decomposition on the molecule provided
  //!
  /*!

      \param mol    - the ROMol of interest

      \return a new ROMol with the Murcko scaffold
              The client is responsible for deleting this molecule.
  */
  ROMol *MurckoDecompose(const ROMol &mol);

  //! \brief Combined two molecules to create a new one
  //!
  /*!

      \param mol1           - the first ROMol to be combined
      \param mol2           - the second ROMol to be combined
      \param offset         - a constant offset to be added to every
                              atom position in mol2 

      \return a new ROMol with the two molecules combined.
              The new molecule has not been sanitized.
              The client is responsible for deleting this molecule.
  */
  ROMol *combineMols(const ROMol &mol1, const ROMol &mol2,
                     RDGeom::Point3D offset=RDGeom::Point3D(0,0,0));
  
  //! \brief Adds named recursive queries to a molecule's atoms based on atom labels
  //!
  /*!

      \param mol            - the molecule to be modified
      \param queries        - the dictionary of named queries to add
      \param propName       - the atom property to use to get query names
      \param reactantLabels - to store pairs of (atom index, query string)


      NOTES:
        - existing query information, if present, will be supplemented (AND logic)
        - non-query atoms will be replaced with query atoms using only the query
          logic
        - query names can be present as comma separated lists, they will then
          be combined using OR logic.
        - throws a KeyErrorException if a particular query name is not present
          in \c queries
      
  */
  void addRecursiveQueries(ROMol &mol,const std::map<std::string,ROMOL_SPTR> &queries,std::string propName,
                           std::vector<std::pair<unsigned int, std::string> > *reactantLabels=NULL);



  //! \brief parses a query definition file and sets up a set of definitions
  //!  suitable for use by addRecursiveQueries()
  /*!

      \param filename         - the name of the file to be read
      \param queryDefs        - the dictionary of named queries (return value)
      \param standardize      - if true, query names will be converted to lower case
      \param delimiter        - the line delimiter in the file
      \param comment          - text used to recognize comment lines
      \param nameColumn       - column with the names of queries
      \param smartsColumn     - column with the SMARTS definitions of the queries

  */
  void parseQueryDefFile(std::string filename,std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize=true,std::string delimiter="\t",std::string comment="//",
                         unsigned int nameColumn=0,unsigned int smartsColumn=1);
  //! \overload
  void parseQueryDefFile(std::istream *inStream,std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize=true,std::string delimiter="\t",std::string comment="//",
                         unsigned int nameColumn=0,unsigned int smartsColumn=1);
  //! \brief equivalent to parseQueryDefFile() but the query definitions are explicitly passed in
  void parseQueryDefText(const std::string &queryDefText,std::map<std::string,ROMOL_SPTR> &queryDefs,
                         bool standardize=true,std::string delimiter="\t",std::string comment="//",
                         unsigned int nameColumn=0,unsigned int smartsColumn=1);
  
}

#endif




