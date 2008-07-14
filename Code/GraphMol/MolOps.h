//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_MOL_OPS_H_
#define _RD_MOL_OPS_H_

#include <RDGeneral/types.h>
#include <boost/tuple/tuple.hpp>
#include <boost/smart_ptr.hpp>

extern const int ci_LOCAL_INF;
namespace RDKit{
  class ROMol;
  class RWMol;
  class Atom;
  class Bond;
  typedef std::vector<double> INVAR_VECT;
  typedef INVAR_VECT::iterator INVAR_VECT_I;
  typedef INVAR_VECT::const_iterator INVAR_VECT_CI;

  //! used to request particular types of atom invariants
  typedef enum {
    ComprehensiveInvariants, //!< for general atomic ranking
    ChiralSearchInvariants,  //!< for ranking atoms for calculating CIP codes
  } AtomInvariantType;

  //! used to return atomic discriminators (three doubles)
  typedef boost::tuples::tuple<double,double,double> DiscrimTuple;


  //! \brief Groups a variety of molecular query and transformation operations.
  namespace MolOps {

    //! return the number of electrons available on an atom to donate for aromaticity
    /*!
       The result is determined using the default valency, number of lone pairs,
       number of bonds and the formal charge. Note that the atom may not donate
       all of these electrons to a ring for aromaticity (also used in Conjugation
       and hybridization code).

       \param at the atom of interest

       \return the number of electrons
    */
    int countAtomElec(const Atom *at);

    //! sums up all atomic formal charges and returns the result
    int getFormalCharge(const ROMol &mol);

    //! returns whether or not the given Atom is involved in a conjugated bond
    bool atomHasConjugatedBond(const Atom *at);

    //! find fragments (disconnected components of the molecular graph)
    /*!

      \param mol     the molecule of interest
      \param mapping used to return the mapping of Atoms->fragments.
         On return \c mapping will be <tt>mol->getNumAtoms()</tt> long
	 and will contain the fragment assignment for each Atom

      \return the number of fragments found.	 
      
    */
    unsigned int getMolFrags(const ROMol &mol,INT_VECT &mapping);
    //! find fragments (disconnected components of the molecular graph)
    /*!

      \param mol    the molecule of interest
      \param frags  used to return the Atoms in each fragment
         On return \c mapping will be \c numFrags long, and each entry
	 will contain the indices of the Atoms in that fragment.

      \return the number of fragments found.	 
      
    */
    unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags);

    //! splits a molecule into its component fragments
    //  (disconnected components of the molecular graph)
    /*!

      \param mol     the molecule of interest
      \param sanitizeFrags  toggles sanitization of the fragments after
                            they are built
      \param frags used to return the mapping of Atoms->fragments.
         if provided, \c frags will be <tt>mol->getNumAtoms()</tt> long 
	 on return and will contain the fragment assignment for each Atom

      \return a vector of the fragments as smart pointers to ROMols
      
    */
    std::vector<boost::shared_ptr<ROMol> > getMolFrags(const ROMol &mol,
                                                       bool sanitizeFrags=true,
                                                       INT_VECT *frags=0);

    //! finds a molecule's minimium spanning tree (MST)
    /*!
      \param mol  the molecule of interest
      \param mst  used to return the MST as a vector of bond indices
    */
    void findSpanningTree(const ROMol &mol,INT_VECT &mst);

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
    DiscrimTuple computeDiscriminators(double *distMat, int nb, int na);
    

    //! calculates Balaban's J index for the molecule
    /*!
      \param mol      the molecule of interest
      \param useBO    toggles inclusion of the bond order in the calculation
                      (when false, we're not really calculating the J value)
      \param force    forces the calculation (instead of using cached results)
      \param bondPath when included, only paths using bonds whose indices occur
                      in this vector will be included in the calculation
      \param cacheIt  If this is true, the calculated value will be cached
                      as a property on the molecule
      \return the J index
      
    */
    double computeBalabanJ(const ROMol &mol, 
				  bool useBO=true,
				  bool force=false,
				  const std::vector<int> *bondPath=0,
				  bool cacheIt=true);
				


    //! \name Dealing with hydrogens
    //{@

    //! returns a copy of a molecule with hydrogens added in as explicit Atoms
    /*!
        \param mol          the molecule to add Hs to
        \param explicitOnly (optional) if this \c true, only explicit Hs will be added
        \param addCoords    (optional) If this is true, estimates for the atomic coordinates
                    of the added Hs will be used.
     
        \return the new molecule 

        <b>Notes:</b>
	   - it makes no sense to use the \c addCoords option if the molecule's heavy
	     atoms don't already have coordinates.
	   - the caller is responsible for <tt>delete</tt>ing the pointer this returns.
     */
    ROMol *addHs(const ROMol &mol,bool explicitOnly=false,bool addCoords=false);


    //! returns a copy of a molecule with hydrogens removed
    /*!
        \param mol          the molecule to remove Hs from
        \param implicitOnly (optional) if this \c true, only implicit Hs will be removed
        \param updateExplicitCount  (optional) If this is \c true, when explicit Hs are removed
	     from the graph, the heavy atom to which they are bound will have its counter of
	     explicit Hs increased.
        \param sanitize:  (optional) If this is \c true, the final molecule will be 
       sanitized
     
        \return the new molecule 

        <b>Notes:</b>
           - Hydrogens which aren't connected to a heavy atom will not be
             removed.  This prevents molecules like <tt>"[H][H]"</tt> from having
             all atoms removed.
           - Labelled hydrogen (e.g. atoms with atomic number=1, but mass > 1),
             will not be removed.

	   - the caller is responsible for <tt>delete</tt>ing the pointer this returns.
    */
    ROMol *removeHs(const ROMol &mol,bool implicitOnly=false,
			   bool updateExplicitCount=false,bool sanitize=true);

    //! returns a copy of a molecule with hydrogens removed and added as queries
    //!  to the heavy atoms to which they are bound.
    /*!
      This is really intended to be used with molecules that contain QueryAtoms
      
      \param mol the molecule to remove Hs from
     
      \return the new molecule 

      <b>Notes:</b>
        - Atoms that do not already have hydrogen count queries will have one
	        added, other H-related queries will not be touched. Examples:
	      - C[H] -> [C;!H0]
	      - [C;H1][H] -> [C;H1]
	      - [C;H2][H] -> [C;H2]
        - Hydrogens which aren't connected to a heavy atom will not be
          removed.  This prevents molecules like <tt>"[H][H]"</tt> from having
          all atoms removed.
	      - the caller is responsible for <tt>delete</tt>ing the pointer this returns.
	
    */
    ROMol *mergeQueryHs(const ROMol &mol);
    //@}

    //! \name Sanitization
    //@{

    //! \brief carries out a collection of tasks for cleaning up a molecule and ensuring
    //! that it makes "chemical sense"
    /*!
       This functions calls the following in sequence
         -# MolOps::cleanUp()
         -# MolOps::Kekulize()
         -# MolOps::setAromaticity()
         -# MolOps::setConjugation()
         -# MolOps::setHybridization()
         -# MolOps::cleanupChirality()
         -# MolOps::adjustHs()
	 
       \param mol the RWMol to be cleaned

       <b>Notes:</b>
        - If there is a failure in the sanitization, a \c SanitException
	  will be thrown.
        - in general the user of this function should cast the molecule following this 
          function to a ROMol, so that new atoms and bonds cannot be added to the 
          molecule and screw up the sanitizing that has been done here
    */
    void sanitizeMol(RWMol &mol);

    //! Sets up the aromaticity for a molecule
    /*!

      This is what happens here:
         -# find all the simple rings by calling the findSSSR function
         -# loop over all the Atoms in each ring and mark them if they are candidates
            for aromaticity. A ring atom is a candidate if it can spare electrons
            to the ring and if it's from the first two rows of the periodic table.
         -# ased on the candidate atoms, mark the rings to be either candidates 
            or non-candidates. A ring is a candidate only if all its atoms are candidates
         -# apply Hueckel rule to each of the candidate rings to check if the ring can be
	    aromatic

      \param mol the RWMol of interest
      
      \return 1 on succes, 0 otherwise
      
      <b>Assumptions:</b>
        - Kekulization has been done (i.e. \c MolOps::Kekulize() has already
	  been called)
       
    */
    int setAromaticity(RWMol &mol);


    //! Designed to be called by the sanitizer to handle special cases before anything is done.
    /*!

        Currently this:
	 - modifies nitro groups, so that the nitrogen does not have a unreasonable
  	   valence of 5, as follows:
             - the nitrogen gets a positve charge
             - one of the oxygens gets a negative chage and the double bond to this 
               oxygen is changed to a single bond
	   The net result is that nitro groups can be counted on to be:
	     \c "[N+](=O)[O-]"

       \param mol    the molecule of interest
     
    */
    void cleanUp(RWMol &mol);


    //! adjust the number of implicit and explicit Hs for special cases
    /*!

        Currently this:
         - modifies aromatic nitrogens so that, when appropriate, they have an
  	   explicit H marked (e.g. so that we get things like \c "c1cc[nH]cc1"

        \param mol    the molecule of interest
	
        <b>Assumptions</b>
           - this is called after the molecule has been sanitized,
             aromaticity has been perceived, and the implicit valence of
             everything has been calculated.

    */
    void adjustHs(RWMol &mol);
    
    //! Kekulizes the molecule
    /*!

       \param mol             the molecule of interest
       \param markAtomsBonds  if this is set to true, \c isAromatic boolean settings
                              on both the Bonds and Atoms are turned to false following
                              the Kekulization, otherwise they are left alone in their 
                              original state.
       \param maxBackTracks   the maximum number of attempts at back-tracking. The algorithm 
                              uses a back-tracking procedure to revist a previous setting of 
                              double bond if we hit a wall in the kekulization process
                              
       <b>Notes:</b>
         - even if \c markAtomsBonds is \c false the \c BondType for all aromatic
	   bonds will be changed from \c RDKit::Bond::AROMATIC to \c RDKit::Bond::SINGLE
	   or RDKit::Bond::DOUBLE during Kekulization.

    */
    void Kekulize(RWMol &mol, bool markAtomsBonds=true, unsigned int maxBackTracks=100);

    //! flags the molecule's conjugated bonds
    void setConjugation(ROMol &mol);

    //! calculates and sets the hybridization of all a molecule's Stoms
    void setHybridization(ROMol &mol);


    // @}

    //! \name Ring finding and SSSR
    //@{

    //! finds a molecule's Smallest Set of Smallest Rings
    /*!
      Currently this implements a modified form of Figueras algorithm
        (JCICS - Vol. 36, No. 5, 1996, 986-991)

      \param mol the molecule of interest
      \param res used to return the vector of rings. Each entry is a vector with
          atom indices.  This information is also stored in the molecule's
	  RingInfo structure, so this argument is optional (see overload)
	  
      \return number of smallest rings found

      Base algorithm:
        - The original algorithm starts by finding representative degree 2
	  nodes. 
	- Representative because if a series of deg 2 nodes are found only
          one of them is picked.
        - The smallest ring around each of them is found. 
        - The bonds that connect to this degree 2 node are them chopped off, yielding
          new deg two nodes 
        - The process is repeated on the new deg 2 nodes.
        - If no deg 2 nodes are found, a deg 3 node is picked. The smallest ring
          with it is found. A bond from this is "carefully" (look in the paper)
          selected and chopped, yielding deg 2 nodes. The process is same as 
          above once this is done.
      
      Our Modifications:
        - If available, more than one smallest ring around a representative deg 2
	  node will be computed and stored
        - Typically 3 rings are found around a degree 3 node (when no deg 2s are available)
	  and all the bond to that node are chopped.
        - The extra rings that were found in this process are removed after all the nodes
	  have been covered.

      These changes were motivated by several factors:
        - We believe the original algorithm fails to find the correct SSSR
          (finds the correct number of them but the wrong ones) on some sample mols
        - Since SSSR may not be unique, a post-SSSR step to symmetrize may be done.
          The extra rings this process adds can be quite useful.
    */
    int findSSSR(const ROMol &mol, VECT_INT_VECT &res);
    //! \overload
    int findSSSR(const ROMol &mol, VECT_INT_VECT *res=0);

    //! symmetrize the molecule's Smallest Set of Smallest Rings
    /*!
       SSSR rings obatined from "findSSSR" can be non-unique in some case.
       For example, cubane has five SSSR rings, not six as one would hope.
       
       This function adds additional rings to the SSSR list if necessary 
       to make the list symmetric, e.g. all atoms in cubane will be part of the same number 
       of SSSRs. This function choses these extra rings from the extra rings computed
       and discarded during findSSSR. The new ring are chosen such that:
        - replacing a same sized ring in the SSSR list with an extra ring yields
          the same union of bond IDs as the orignal SSSR list
     
      \param mol - the molecule of interest
      \param res used to return the vector of rings. Each entry is a vector with
          atom indices.  This information is also stored in the molecule's
	  RingInfo structure, so this argument is optional (see overload)
      
      \return the total number of rings = (new rings + old SSSRs)
     
      <b>Notes:</b>
       - if no SSSR rings are found on the molecule - MolOps::findSSSR() is called first
    */
    int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res);

    //@}

    //! \name Shortest paths and other matrices
    //@{

    //! returns a molecule's adjacency matrix
    /*!
      \param mol             the molecule of interest
      \param useBO           toggles use of bond orders in the matrix
      \param emptyVal        sets the empty value (for non-adjacent atoms)
      \param force           forces calculation of the matrix, even if already computed
      \param propNamePrefix  used to set the cached property name

      \return the adjacency matrix.

      <b>Notes</b>
        - The result of this is cached in the molecule's local property dictionary,
	  which will handle deallocation. Do the caller should <b>not</b> \c delete
	  this pointer.
	  
    */
    double * getAdjacencyMatrix(const ROMol &mol,
				       bool useBO=false,
				       int emptyVal=0,
				       bool force=false,
				       const char *propNamePrefix=0);

    //! Computes the molecule's topological distance matrix
    /*!
       Uses the Floyd-Warshall all-pairs-shortest-paths algorithm.
      
      \param mol             the molecule of interest
      \param useBO           toggles use of bond orders in the matrix
      \param useAtomWts      sets the diagonal elements of the result to
               6.0/(atomic number) so that the matrix can be used to calculate
	       Balaban J values.  This does not affect the bond weights. 
      \param force           forces calculation of the matrix, even if already computed
      \param propNamePrefix  used to set the cached property name

      \return the distance matrix.

      <b>Notes</b>
        - The result of this is cached in the molecule's local property dictionary,
	  which will handle deallocation. Do the caller should <b>not</b> \c delete
	  this pointer.
	  
     
    */
    double *getDistanceMat(const ROMol &mol,
				  bool useBO=false,
				  bool useAtomWts=false,
				  bool force=false,
				  const char *propNamePrefix=0);


    //! Computes the molecule's topological distance matrix
    /*!
       Uses the Floyd-Warshall all-pairs-shortest-paths algorithm.
      
      \param mol             the molecule of interest
      \param activeAtoms     only elements corresponding to these atom indices
                             will be included in the calculation
      \param bonds           only bonds found in this list will be included in the
                             calculation
      \param useBO           toggles use of bond orders in the matrix
      \param useAtomWts      sets the diagonal elements of the result to
               6.0/(atomic number) so that the matrix can be used to calculate
	       Balaban J values.  This does not affect the bond weights. 

      \return the distance matrix.

      <b>Notes</b>
        - The results of this call are not cached, the caller <b>should</b> \c delete
	  this pointer.
	  
     
    */
    double *getDistanceMat(const ROMol &mol,
				  const std::vector<int> &activeAtoms,
				  const std::vector<const Bond *> &bonds,
				  bool useBO=false,
				  bool useAtomWts=false);


    //! Find the shortest path between two atoms
    /*!
      Uses the Bellman-Ford algorithm

     \param mol  molecule of interest
     \param aid1 index of the first atom
     \param aid2 index of the second atom

     \return an std::list with the indices of the atoms along the shortest
        path

     <b>Notes:</b>
       - the starting and end atoms are not include in the path

    */
    INT_LIST getShortestPath(const ROMol &mol, int aid1, int aid2);

    //@}
    
    //! \name Canonicalization
    //@{

    //! assign a canonical ordering to a molecule's atoms
    /*!
      The algorithm used here is a modification of the published Daylight canonical
      smiles algorithm (i.e. it uses atom invariants and products of primes).
      
      \param mol               the molecule of interest
      \param ranks             used to return the ranks
      \param invariantType     the type of invariant to be computed 
      \param breakTies         toggles breaking of ties (see below)
      \param rankHistory       used to return the rank history (see below)
      \param includeChirality  toggles inclusion of chirality in the atom ranking
                               (see below)

      <b>Notes:</b>
        - Tie breaking should be done when it's important to have a full ordering
          of the atoms (e.g. when generating canonical traversal trees). If it's
	        acceptable to have ties between symmetry-equivalent atoms (e.g. when
	        generating CIP codes), tie breaking can/should be skipped.
	      - if the \c rankHistory argument is provided, the evolution of the ranks of
	        individual atoms will be tracked.  The \c rankHistory pointer should be
	        to a VECT_INT_VECT that has at least \c mol.getNumAtoms() elements.
        - the \c includeChirality argument only makes sense when generating full
          invariants (i.e. ComprehensiveInvariant)
    */
    void rankAtoms(const ROMol &mol,INT_VECT &ranks,
			  AtomInvariantType invariantType=ComprehensiveInvariants,
			  bool breakTies=true,VECT_INT_VECT *rankHistory=0,
			  bool includeChirality=false);

    //! calculates a set of atom invariants
    /*!
      These invariants include terms for:
        -# degree
        -# explicit valence
        -# atomic number
        -# isotope
        -# number of hs
        -# size of rings an atom is involved in
        -# formal charge
        -# chiral code (optional)
        -# whether of not the atom determines the stereochemistry of a
        	   double bond. (optional)
      
      \param mol               the molecule of interest
      \param res               used to return the results
      \param includeChirality  toggles inclusion of chirality, see MolOps::rankAtoms()
                                for an explanation.
    */
    void buildAtomInvariants(const ROMol &mol,INVAR_VECT &res,
				    bool includeChirality=false);
    @}

    //! \name Stereochemistry
    //@{

    //! removes bogus chirality markers (those on non-sp3 centers):
    void cleanupChirality(RWMol &mol);

    //! \brief Uses a conformer to assign ChiralType to a molecule's atoms
    /*!
      \param mol                  the molecule of interest
      \param confId               the conformer to use
      \param replaceExistingTags  if this flag is true, any existing atomic chiral
                                  tags will be replaced

      If the conformer provided is not a 3D conformer, nothing will be done.
    */
    void assignChiralTypesFrom3D(ROMol &mol,int confId=-1,bool replaceExistingTags=true);

    //! calculates a set of chiral atom invariants
    /*!
      These are based upon Labute's proposal in the article:
      "An Efficient Algorithm for the Determination of Topological
      RS Chirality" <i>Journal of the CCG</i> (1996)
      and include terms for:
        -# atomic number
        -# isotope
        -# number of multiple bonds from/to this atom
        -# degree
        -# number of hydrogens

      \param mol               the molecule of interest
      \param res               used to return the results
    */
    void buildChiralAtomInvariants(const ROMol &mol,
					  INVAR_VECT &res);


    //! Figure out the CIP ranks for the atoms of a molecule
    /*!
      \param mol the molecule to be altered
      \param ranks  used to return the set of ranks.  
                    Should be at least mol.getNumAtoms() long.
    
      <b>Notes:</b>
         - All atoms gain a property "_CIPRank" with their overall
           CIP ranking.
    
    */
    void assignAtomCIPRanks(const ROMol &mol,INT_VECT &ranks);

    //! calculate CIP (R/S) codes for a molecule's chiral centers
    /*!
      \param mol     the molecule of interest
      \param cleanIt toggles removal of chiral flags from atoms that do
                     not have four unique neighbors
      \param force   forces the calculation to be repeated even if it has 
                     already been done 

      <b>Notes:</b>
        - The codes are stored on atoms as a string property named
          "_CIPCode".
        - Throughout we assume that we're working with a hydrogen-suppressed
          graph.
    */
    void assignAtomChiralCodes(ROMol &mol,bool cleanIt=false,bool force=false);
    //! calculate CIP (Z/E) codes for a molecule's double bonds
    /*!
      \param mol     the molecule of interest
      \param cleanIt toggles removal of stereo flags from double bonds that can
                     not have stereochemistry
      \param force   forces the calculation to be repeated even if it has 
                     already been done 

      <b>Notes:M</b>
        - Throughout we assume that we're working with a hydrogen-suppressed
          graph.

    */
    void assignBondStereoCodes(ROMol &mol,bool cleanIt=false,bool force=false);
    //! \brief finds bonds that could be cis/trans in a molecule and mark them as
    //!  Bond::STEREONONE
    /*!
      \param mol     the molecule of interest
      \param cleanIt toggles removal of stereo flags from double bonds that can
                     not have stereochemistry

      This function is usefuly in two situations
        - when parsing a mol file; for the bonds marked here, coordinate informations 
	  on the neighbors can be used to indentify cis or trans states
        - when writing a mol file; bonds that can be cis/trans but not marked as either 
	  need to be specially marked in the mol file
    */
    void findPotentialStereoBonds(ROMol &mol,bool cleanIt=false);
    //@}

  }; // end of namespace MolOps
}; // end of namespace RDKit

#endif
