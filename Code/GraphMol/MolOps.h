//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOL_OPS_H
#define RD_MOL_OPS_H

#include <vector>
#include <map>
#include <list>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/smart_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/types.h>
#include "SanitException.h"
#include <RDGeneral/FileParseException.h>

RDKIT_GRAPHMOL_EXPORT extern const int ci_LOCAL_INF;
namespace RDKit {
class ROMol;
class RWMol;
class Atom;
class Bond;
class Conformer;
typedef std::vector<double> INVAR_VECT;
typedef INVAR_VECT::iterator INVAR_VECT_I;
typedef INVAR_VECT::const_iterator INVAR_VECT_CI;

//! \brief Groups a variety of molecular query and transformation operations.
namespace MolOps {

//! return the number of electrons available on an atom to donate for
/// aromaticity
/*!
   The result is determined using the default valency, number of lone pairs,
   number of bonds and the formal charge. Note that the atom may not donate
   all of these electrons to a ring for aromaticity (also used in Conjugation
   and hybridization code).

   \param at the atom of interest

   \return the number of electrons
*/
RDKIT_GRAPHMOL_EXPORT int countAtomElec(const Atom *at);

//! sums up all atomic formal charges and returns the result
RDKIT_GRAPHMOL_EXPORT int getFormalCharge(const ROMol &mol);

//! returns whether or not the given Atom is involved in a conjugated bond
RDKIT_GRAPHMOL_EXPORT bool atomHasConjugatedBond(const Atom *at);

//! find fragments (disconnected components of the molecular graph)
/*!

  \param mol     the molecule of interest
  \param mapping used to return the mapping of Atoms->fragments.
     On return \c mapping will be <tt>mol->getNumAtoms()</tt> long
     and will contain the fragment assignment for each Atom

  \return the number of fragments found.

*/
RDKIT_GRAPHMOL_EXPORT unsigned int getMolFrags(const ROMol &mol,
                                               std::vector<int> &mapping);
//! find fragments (disconnected components of the molecular graph)
/*!

  \param mol    the molecule of interest
  \param frags  used to return the Atoms in each fragment
     On return \c mapping will be \c numFrags long, and each entry
     will contain the indices of the Atoms in that fragment.

  \return the number of fragments found.

*/
RDKIT_GRAPHMOL_EXPORT unsigned int getMolFrags(
    const ROMol &mol, std::vector<std::vector<int>> &frags);

//! splits a molecule into its component fragments
/// (disconnected components of the molecular graph)
/*!

  \param mol     the molecule of interest
  \param sanitizeFrags  toggles sanitization of the fragments after
                        they are built
  \param frags used to return the mapping of Atoms->fragments.
     if provided, \c frags will be <tt>mol->getNumAtoms()</tt> long
         on return and will contain the fragment assignment for each Atom
  \param fragsMolAtomMapping  used to return the Atoms in each fragment
     On return \c mapping will be \c numFrags long, and each entry
     will contain the indices of the Atoms in that fragment.
   \param copyConformers  toggles copying conformers of the fragments after
                        they are built
  \return a vector of the fragments as smart pointers to ROMols

*/
RDKIT_GRAPHMOL_EXPORT std::vector<boost::shared_ptr<ROMol>> getMolFrags(
    const ROMol &mol, bool sanitizeFrags = true,
    std::vector<int> *frags = nullptr,
    std::vector<std::vector<int>> *fragsMolAtomMapping = nullptr,
    bool copyConformers = true);

//! splits a molecule into pieces based on labels assigned using a query
/*!

  \param mol     the molecule of interest
  \param query   the query used to "label" the molecule for fragmentation
  \param sanitizeFrags  toggles sanitization of the fragments after
                        they are built
  \param whiteList  if provided, only labels in the list will be kept
  \param negateList if true, the white list logic will be inverted: only labels
                    not in the list will be kept

  \return a map of the fragments and their labels

*/
template <typename T>
RDKIT_GRAPHMOL_EXPORT std::map<T, boost::shared_ptr<ROMol>>
getMolFragsWithQuery(const ROMol &mol, T (*query)(const ROMol &, const Atom *),
                     bool sanitizeFrags = true,
                     const std::vector<T> *whiteList = nullptr,
                     bool negateList = false);

#if 0
    //! finds a molecule's minimum spanning tree (MST)
    /*!
      \param mol  the molecule of interest
      \param mst  used to return the MST as a vector of bond indices
    */
    RDKIT_GRAPHMOL_EXPORT void findSpanningTree(const ROMol &mol,std::vector<int> &mst);
#endif

//! \name Dealing with hydrogens
//{@

//! returns a copy of a molecule with hydrogens added in as explicit Atoms
/*!
    \param mol          the molecule to add Hs to
    \param explicitOnly (optional) if this \c true, only explicit Hs will be
   added
    \param addCoords    (optional) If this is true, estimates for the atomic
   coordinates
                of the added Hs will be used.
    \param onlyOnAtoms   (optional) if provided, this should be a vector of
                IDs of the atoms that will be considered for H addition.
    \param addResidueInfo   (optional) if this is true, add residue info to
                hydrogen atoms (useful for PDB files).

    \return the new molecule

    <b>Notes:</b>
       - it makes no sense to use the \c addCoords option if the molecule's
   heavy
         atoms don't already have coordinates.
       - the caller is responsible for <tt>delete</tt>ing the pointer this
   returns.
 */
RDKIT_GRAPHMOL_EXPORT ROMol *addHs(const ROMol &mol, bool explicitOnly = false,
                                   bool addCoords = false,
                                   const UINT_VECT *onlyOnAtoms = nullptr,
                                   bool addResidueInfo = false);
//! \overload
/// modifies the molecule in place
RDKIT_GRAPHMOL_EXPORT void addHs(RWMol &mol, bool explicitOnly = false,
                                 bool addCoords = false,
                                 const UINT_VECT *onlyOnAtoms = nullptr,
                                 bool addResidueInfo = false);

//! Sets Cartesian coordinates for a terminal atom.
//! Useful for growing an atom off a molecule with sensible
//! coordinates based on the geometry of the neighbor.
/*!
    NOTE: this sets appropriate coordinates in all of the molecule's conformers.
    \param mol       the molecule the atoms belong to
    \param idx       index of the terminal atom whose coordinates are set
    \param otherIdx  index of the bonded neighbor atom
*/

RDKIT_GRAPHMOL_EXPORT void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
                                                 unsigned int otherIdx);

//! returns a copy of a molecule with hydrogens removed
/*!
    \param mol          the molecule to remove Hs from
    \param implicitOnly (optional) if this \c true, only implicit Hs will be
   removed
    \param updateExplicitCount  (optional) If this is \c true, when explicit Hs
   are removed
         from the graph, the heavy atom to which they are bound will have its
   counter of
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
       - two coordinate Hs, like the central H in C[H-]C, will not be removed
       - Hs connected to dummy atoms will not be removed
       - Hs that are part of the definition of double bond Stereochemistry
         will not be removed
       - Hs that are not connected to anything else will not be removed
       - Hs that have a query defined (i.e. hasQuery() returns true) will not
         be removed

       - the caller is responsible for <tt>delete</tt>ing the pointer this
   returns.
*/

RDKIT_GRAPHMOL_EXPORT ROMol *removeHs(const ROMol &mol,
                                      bool implicitOnly = false,
                                      bool updateExplicitCount = false,
                                      bool sanitize = true);
//! \overload
/// modifies the molecule in place
RDKIT_GRAPHMOL_EXPORT void removeHs(RWMol &mol, bool implicitOnly = false,
                                    bool updateExplicitCount = false,
                                    bool sanitize = true);
struct RDKIT_GRAPHMOL_EXPORT RemoveHsParameters {
  bool removeDegreeZero = false;    /**< hydrogens that have no bonds */
  bool removeHigherDegrees = false; /**< hydrogens with two (or more) bonds */
  bool removeOnlyHNeighbors =
      false; /**< hydrogens with bonds only to other hydrogens */
  bool removeIsotopes = false; /**< hydrogens with non-default isotopes */
  bool removeAndTrackIsotopes = false; /**< removes hydrogens with non-default
   isotopes and keeps track of the heavy atom the isotopes were attached to in
   the private _isotopicHs atom property, so they are re-added by AddHs() as the
   original isotopes if possible*/
  bool removeDummyNeighbors =
      false; /**< hydrogens with at least one dummy-atom neighbor */
  bool removeDefiningBondStereo =
      false; /**< hydrogens defining bond stereochemistry */
  bool removeWithWedgedBond = true; /**< hydrogens with wedged bonds to them */
  bool removeWithQuery = false;     /**< hydrogens with queries defined */
  bool removeMapped = true;         /**< mapped hydrogens */
  bool removeInSGroups = true;      /**< part of a SubstanceGroup.
    An H atom will only be removed if it doesn't cause any SGroup to become empty,
    and if it doesn't play a special role in the SGroup (XBOND, attach point
    or a CState) */
  bool showWarnings = true; /**< display warnings for Hs that are not removed */
  bool removeNonimplicit = true; /**< DEPRECATED equivalent of !implicitOnly */
  bool updateExplicitCount =
      false; /**< DEPRECATED equivalent of updateExplicitCount */
  bool removeHydrides = true; /**< Removing Hydrides */
  bool removeNontetrahedralNeighbors =
      false; /**<  remove Hs which are bonded to atoms with specified
                non-tetrahedral stereochemistry */
};
//! \overload
/// modifies the molecule in place
RDKIT_GRAPHMOL_EXPORT void removeHs(RWMol &mol, const RemoveHsParameters &ps,
                                    bool sanitize = true);
//! \overload
/// The caller owns the pointer this returns
RDKIT_GRAPHMOL_EXPORT ROMol *removeHs(const ROMol &mol,
                                      const RemoveHsParameters &ps,
                                      bool sanitize = true);

//! removes all Hs from a molecule
RDKIT_GRAPHMOL_EXPORT void removeAllHs(RWMol &mol, bool sanitize = true);
//! \overload
/// The caller owns the pointer this returns
RDKIT_GRAPHMOL_EXPORT ROMol *removeAllHs(const ROMol &mol,
                                         bool sanitize = true);

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
    - the caller is responsible for <tt>delete</tt>ing the pointer this
  returns.
    - By default all hydrogens are removed, however if
      mergeUnmappedOnly is true, any hydrogen participating
      in an atom map will be retained

*/
RDKIT_GRAPHMOL_EXPORT ROMol *mergeQueryHs(const ROMol &mol,
                                          bool mergeUnmappedOnly = false,
                                          bool mergeIsotopes = false);
//! \overload
/// modifies the molecule in place
RDKIT_GRAPHMOL_EXPORT void mergeQueryHs(RWMol &mol,
                                        bool mergeUnmappedOnly = false,
                                        bool mergeIsotopes = false);

//! returns a pair of booleans (hasQueryHs, hasUnmergaebleQueryHs)
/*!
  This is really intended to be used with molecules that contain QueryAtoms
  such as when checking smarts patterns for explicit hydrogens


  \param mol the molecule to check for query Hs from
  \return std::pair  if pair.first is true if the molecule has query hydrogens,
  if pair.second is true, the queryHs cannot be removed my mergeQueryHs
*/
RDKIT_GRAPHMOL_EXPORT std::pair<bool, bool> hasQueryHs(const ROMol &mol);

typedef enum {
  ADJUST_IGNORENONE = 0x0,
  ADJUST_IGNORECHAINS = 0x1,
  ADJUST_IGNORERINGS = 0x4,
  ADJUST_IGNOREDUMMIES = 0x2,
  ADJUST_IGNORENONDUMMIES = 0x8,
  ADJUST_IGNOREMAPPED = 0x10,
  ADJUST_IGNOREALL = 0xFFFFFFF
} AdjustQueryWhichFlags;

//! Parameters controlling the behavior of MolOps::adjustQueryProperties
/*!

  Note that some of the options here are either directly contradictory or make
  no sense when combined with each other. We generally assume that client code
  is doing something sensible and don't attempt to detect possible conflicts or
  problems.

*/
struct RDKIT_GRAPHMOL_EXPORT AdjustQueryParameters {
  bool adjustDegree = true; /**< add degree queries */
  std::uint32_t adjustDegreeFlags = ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS;

  bool adjustRingCount = false; /**< add ring-count queries */
  std::uint32_t adjustRingCountFlags =
      ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS;

  bool makeDummiesQueries = true; /**< convert dummy atoms without isotope
                                labels to any-atom queries */

  bool aromatizeIfPossible = true; /**< perceive and set aromaticity */

  bool makeBondsGeneric =
      false; /**< convert bonds to generic queries (any bonds) */
  std::uint32_t makeBondsGenericFlags = ADJUST_IGNORENONE;

  bool makeAtomsGeneric =
      false; /**< convert atoms to generic queries (any atoms) */
  std::uint32_t makeAtomsGenericFlags = ADJUST_IGNORENONE;

  bool adjustHeavyDegree = false; /**< adjust the heavy-atom degree instead of
                               overall degree */
  std::uint32_t adjustHeavyDegreeFlags =
      ADJUST_IGNOREDUMMIES | ADJUST_IGNORECHAINS;

  bool adjustRingChain = false; /**< add ring-chain queries */
  std::uint32_t adjustRingChainFlags = ADJUST_IGNORENONE;

  bool useStereoCareForBonds =
      false; /**< remove stereochemistry info from double bonds that do not have
                the stereoCare property set */

  bool adjustConjugatedFiveRings =
      false; /**< sets bond queries in conjugated five-rings to
                SINGLE|DOUBLE|AROMATIC */

  bool setMDLFiveRingAromaticity =
      false; /**< uses the 5-ring aromaticity behavior of the (former) MDL
                software as documented in the Chemical Representation Guide */

  bool adjustSingleBondsToDegreeOneNeighbors =
      false; /**<  sets single bonds between aromatic or conjugated atoms and
                degree one neighbors to SINGLE|AROMATIC */

  bool adjustSingleBondsBetweenAromaticAtoms =
      false; /**<  sets non-ring single bonds between two aromatic or conjugated
                atoms to SINGLE|AROMATIC */

  //! \brief returns an AdjustQueryParameters object with all adjustments
  //! disabled
  static AdjustQueryParameters noAdjustments() {
    AdjustQueryParameters res;
    res.adjustDegree = false;
    res.makeDummiesQueries = false;
    res.aromatizeIfPossible = false;
    return res;
  }
  AdjustQueryParameters() {}
};

//! updates an AdjustQueryParameters object from a JSON string
RDKIT_GRAPHMOL_EXPORT void parseAdjustQueryParametersFromJSON(
    MolOps::AdjustQueryParameters &p, const std::string &json);

//! returns a copy of a molecule with query properties adjusted
/*!
  \param mol the molecule to adjust
  \param params controls the adjustments made

  \return the new molecule, the caller owns the memory
*/
RDKIT_GRAPHMOL_EXPORT ROMol *adjustQueryProperties(
    const ROMol &mol, const AdjustQueryParameters *params = nullptr);
//! \overload
/// modifies the molecule in place
RDKIT_GRAPHMOL_EXPORT void adjustQueryProperties(
    RWMol &mol, const AdjustQueryParameters *params = nullptr);

//! returns a copy of a molecule with the atoms renumbered
/*!

  \param mol the molecule to work with
  \param newOrder the new ordering of the atoms (should be numAtoms long)
     for example: if newOrder is [3,2,0,1], then atom 3 in the original
     molecule will be atom 0 in the new one

  \return the new molecule

  <b>Notes:</b>
    - the caller is responsible for <tt>delete</tt>ing the pointer this
  returns.

*/
RDKIT_GRAPHMOL_EXPORT ROMol *renumberAtoms(
    const ROMol &mol, const std::vector<unsigned int> &newOrder);

//! @}

//! \name Sanitization
/// {

typedef enum {
  SANITIZE_NONE = 0x0,
  SANITIZE_CLEANUP = 0x1,
  SANITIZE_PROPERTIES = 0x2,
  SANITIZE_SYMMRINGS = 0x4,
  SANITIZE_KEKULIZE = 0x8,
  SANITIZE_FINDRADICALS = 0x10,
  SANITIZE_SETAROMATICITY = 0x20,
  SANITIZE_SETCONJUGATION = 0x40,
  SANITIZE_SETHYBRIDIZATION = 0x80,
  SANITIZE_CLEANUPCHIRALITY = 0x100,
  SANITIZE_ADJUSTHS = 0x200,
  SANITIZE_CLEANUP_ORGANOMETALLICS = 0x400,
  SANITIZE_CLEANUPATROPISOMERS = 0x800,
  SANITIZE_ALL = 0xFFFFFFF
} SanitizeFlags;

//! \brief carries out a collection of tasks for cleaning up a molecule and
//! ensuring that it makes "chemical sense"
/*!
   This functions calls the following in sequence
     -# MolOps::cleanUp()
     -# mol.updatePropertyCache()
     -# MolOps::symmetrizeSSSR()
     -# MolOps::Kekulize()
     -# MolOps::assignRadicals()
     -# MolOps::setAromaticity()
     -# MolOps::setConjugation()
     -# MolOps::setHybridization()
     -# MolOps::cleanupChirality()
     -# MolOps::adjustHs()
     -# mol.updatePropertyCache()

   \param mol : the RWMol to be cleaned

   \param operationThatFailed : the first (if any) sanitization operation that
                                fails is set here.
                                The values are taken from the \c SanitizeFlags
                                enum. On success, the value is \c
                                SanitizeFlags::SANITIZE_NONE

   \param sanitizeOps : the bits here are used to set which sanitization
                        operations are carried out. The elements of the \c
                        SanitizeFlags enum define the operations.

   <b>Notes:</b>
    - If there is a failure in the sanitization, a \c MolSanitizeException
      will be thrown.
    - in general the user of this function should cast the molecule following
      this function to a ROMol, so that new atoms and bonds cannot be added to
      the molecule and screw up the sanitizing that has been done here
*/
RDKIT_GRAPHMOL_EXPORT void sanitizeMol(RWMol &mol,
                                       unsigned int &operationThatFailed,
                                       unsigned int sanitizeOps = SANITIZE_ALL);
//! \overload
RDKIT_GRAPHMOL_EXPORT void sanitizeMol(RWMol &mol);

//! \brief Identifies chemistry problems (things that don't make chemical
//! sense) in a molecule
/*!
   This functions uses the operations in sanitizeMol but does not change
   the input structure and returns a list of the problems encountered instead
   of stopping at the first failure,

   The problems this looks for come from the sanitization operations:
     -# mol.updatePropertyCache()  : Unreasonable valences
     -# MolOps::Kekulize()  : Unkekulizable ring systems, aromatic atoms not
   in rings, aromatic bonds to non-aromatic atoms.

   \param mol : the ROMol to be cleaned

   \param sanitizeOps : the bits here are used to set which sanitization
                        operations are carried out. The elements of the \c
                        SanitizeFlags enum define the operations.

   \return a vector of \c MolSanitizeException values that indicate what
           problems were encountered

*/
RDKIT_GRAPHMOL_EXPORT
std::vector<std::unique_ptr<MolSanitizeException>> detectChemistryProblems(
    const ROMol &mol, unsigned int sanitizeOps = SANITIZE_ALL);

//! Possible aromaticity models
/*!
- \c AROMATICITY_DEFAULT at the moment always uses \c AROMATICITY_RDKIT
- \c AROMATICITY_RDKIT is the standard RDKit model (as documented in the RDKit
Book)
- \c AROMATICITY_SIMPLE only considers 5- and 6-membered simple rings (it
does not consider the outer envelope of fused rings)
- \c AROMATICITY_MDL
- \c AROMATICITY_CUSTOM uses a caller-provided function
*/
typedef enum {
  AROMATICITY_DEFAULT = 0x0,  ///< future proofing
  AROMATICITY_RDKIT = 0x1,
  AROMATICITY_SIMPLE = 0x2,
  AROMATICITY_MDL = 0x4,
  AROMATICITY_CUSTOM = 0xFFFFFFF  ///< use a function
} AromaticityModel;

//! Sets up the aromaticity for a molecule
/*!

  This is what happens here:
     -# find all the simple rings by calling the findSSSR function
     -# loop over all the Atoms in each ring and mark them if they are
  candidates
        for aromaticity. A ring atom is a candidate if it can spare electrons
        to the ring and if it's from the first two rows of the periodic table.
     -# based on the candidate atoms, mark the rings to be either candidates
        or non-candidates. A ring is a candidate only if all its atoms are
  candidates
     -# apply Hueckel rule to each of the candidate rings to check if the ring
  can be
        aromatic

  \param mol the RWMol of interest
  \param model the aromaticity model to use
  \param func a custom function for assigning aromaticity (only used when
  model=\c AROMATICITY_CUSTOM)

  \return >0 on success, <= 0 otherwise

  <b>Assumptions:</b>
    - Kekulization has been done (i.e. \c MolOps::Kekulize() has already
      been called)

*/
RDKIT_GRAPHMOL_EXPORT int setAromaticity(
    RWMol &mol, AromaticityModel model = AROMATICITY_DEFAULT,
    int (*func)(RWMol &) = nullptr);

//! Designed to be called by the sanitizer to handle special cases before
/// anything is done.
/*!

    Currently this:
     - modifies nitro groups, so that the nitrogen does not have an
   unreasonable valence of 5, as follows:
         - the nitrogen gets a positive charge
         - one of the oxygens gets a negative chage and the double bond to
   this oxygen is changed to a single bond The net result is that nitro groups
   can be counted on to be: \c "[N+](=O)[O-]"
     - modifies halogen-oxygen containing species as follows:
        \c [Cl,Br,I](=O)(=O)(=O)O -> [X+3]([O-])([O-])([O-])O
        \c [Cl,Br,I](=O)(=O)O -> [X+3]([O-])([O-])O
        \c [Cl,Br,I](=O)O -> [X+]([O-])O
     - converts the substructure [N,C]=P(=O)-* to [N,C]=[P+](-[O-])-*

   \param mol    the molecule of interest

*/
RDKIT_GRAPHMOL_EXPORT void cleanUp(RWMol &mol);

//! Designed to be called by the sanitizer to handle special cases for
//! organometallic species before valence is perceived
/*!

    \b Note that this function is experimental and may either change in behavior
   or be replaced with something else in future releases.

    Currently this:
     - replaces single bonds between "hypervalent" organic atoms and metals with
       dative bonds (this is following an IUPAC recommendation:
       https://iupac.qmul.ac.uk/tetrapyrrole/TP8.html)

   \param mol    the molecule of interest

*/
RDKIT_GRAPHMOL_EXPORT void cleanUpOrganometallics(RWMol &mol);

//! Called by the sanitizer to assign radical counts to atoms
RDKIT_GRAPHMOL_EXPORT void assignRadicals(RWMol &mol);

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
RDKIT_GRAPHMOL_EXPORT void adjustHs(RWMol &mol);

//! Kekulizes the molecule
/*!

   \param mol             the molecule of interest

   \param markAtomsBonds  if this is set to true, \c isAromatic boolean
   settings on both the Bonds and Atoms are turned to false following the
   Kekulization, otherwise they are left alone in their original state.

   \param maxBackTracks   the maximum number of attempts at back-tracking. The
   algorithm uses a back-tracking procedure to revisit a previous setting of
   double bond if we hit a wall in the kekulization process

   <b>Notes:</b>
     - this does not modify query bonds which have bond type queries (like
   those which come from SMARTS) or rings containing them.
     - even if \c markAtomsBonds is \c false the \c BondType for all modified
       aromatic bonds will be changed from \c RDKit::Bond::AROMATIC to \c
       RDKit::Bond::SINGLE or RDKit::Bond::DOUBLE during Kekulization.

*/
RDKIT_GRAPHMOL_EXPORT void Kekulize(RWMol &mol, bool markAtomsBonds = true,
                                    unsigned int maxBackTracks = 100);
//! Kekulizes the molecule if possible. If the kekulization fails the molecule
//! will not be modified
/*!

   \param mol             the molecule of interest

   \param markAtomsBonds  if this is set to true, \c isAromatic boolean
   settings on both the Bonds and Atoms are turned to false following the
   Kekulization, otherwise they are left alone in their original state.

   \param maxBackTracks   the maximum number of attempts at back-tracking. The
   algorithm uses a back-tracking procedure to revisit a previous setting of
   double bond if we hit a wall in the kekulization process

   \returns whether or not the kekulization succeeded

   <b>Notes:</b>
     - even if \c markAtomsBonds is \c false the \c BondType for all aromatic
       bonds will be changed from \c RDKit::Bond::AROMATIC to \c
       RDKit::Bond::SINGLE or RDKit::Bond::DOUBLE during Kekulization.

*/
RDKIT_GRAPHMOL_EXPORT bool KekulizeIfPossible(RWMol &mol,
                                              bool markAtomsBonds = true,
                                              unsigned int maxBackTracks = 100);

//! flags the molecule's conjugated bonds
RDKIT_GRAPHMOL_EXPORT void setConjugation(ROMol &mol);

//! calculates and sets the hybridization of all a molecule's Stoms
RDKIT_GRAPHMOL_EXPORT void setHybridization(ROMol &mol);

//!  @}

//! \name Ring finding and SSSR
//! @{

//! finds a molecule's Smallest Set of Smallest Rings
/*!
  Currently this implements a modified form of Figueras algorithm
    (JCICS - Vol. 36, No. 5, 1996, 986-991)

  \param mol the molecule of interest
  \param res used to return the vector of rings. Each entry is a vector with
      atom indices.  This information is also stored in the molecule's
      RingInfo structure, so this argument is optional (see overload)
  \param includeDativeBonds - determines whether or not dative bonds are used in
      the ring finding.

  \return number of smallest rings found

  Base algorithm:
    - The original algorithm starts by finding representative degree 2
      nodes.
    - Representative because if a series of deg 2 nodes are found only
      one of them is picked.
    - The smallest ring around each of them is found.
    - The bonds that connect to this degree 2 node are them chopped off,
  yielding
      new deg two nodes
    - The process is repeated on the new deg 2 nodes.
    - If no deg 2 nodes are found, a deg 3 node is picked. The smallest ring
      with it is found. A bond from this is "carefully" (look in the paper)
      selected and chopped, yielding deg 2 nodes. The process is same as
      above once this is done.

  Our Modifications:
    - If available, more than one smallest ring around a representative deg 2
      node will be computed and stored
    - Typically 3 rings are found around a degree 3 node (when no deg 2s are
  available)
      and all the bond to that node are chopped.
    - The extra rings that were found in this process are removed after all
  the nodes have been covered.

  These changes were motivated by several factors:
    - We believe the original algorithm fails to find the correct SSSR
      (finds the correct number of them but the wrong ones) on some sample
  mols
    - Since SSSR may not be unique, a post-SSSR step to symmetrize may be
  done. The extra rings this process adds can be quite useful.
*/
RDKIT_GRAPHMOL_EXPORT int findSSSR(const ROMol &mol,
                                   std::vector<std::vector<int>> &res,
                                   bool includeDativeBonds = false);
//! \overload
RDKIT_GRAPHMOL_EXPORT int findSSSR(const ROMol &mol,
                                   std::vector<std::vector<int>> *res = nullptr,
                                   bool includeDativeBonds = false);

//! use a DFS algorithm to identify ring bonds and atoms in a molecule
/*!
  \b NOTE: though the RingInfo structure is populated by this function,
  the only really reliable calls that can be made are to check if
  mol.getRingInfo().numAtomRings(idx) or mol.getRingInfo().numBondRings(idx)
  return values >0
*/
RDKIT_GRAPHMOL_EXPORT void fastFindRings(const ROMol &mol);

RDKIT_GRAPHMOL_EXPORT void findRingFamilies(const ROMol &mol);

//! symmetrize the molecule's Smallest Set of Smallest Rings
/*!
   SSSR rings obatined from "findSSSR" can be non-unique in some case.
   For example, cubane has five SSSR rings, not six as one would hope.

   This function adds additional rings to the SSSR list if necessary
   to make the list symmetric, e.g. all atoms in cubane will be part of the
  same number of SSSRs. This function choses these extra rings from the extra
  rings computed and discarded during findSSSR. The new ring are chosen such
  that:
    - replacing a same sized ring in the SSSR list with an extra ring yields
      the same union of bond IDs as the original SSSR list

  \param mol - the molecule of interest
  \param res used to return the vector of rings. Each entry is a vector with
      atom indices.  This information is also stored in the molecule's
      RingInfo structure, so this argument is optional (see overload)
  \param includeDativeBonds - determines whether or not dative bonds are used in
      the ring finding.

  \return the total number of rings = (new rings + old SSSRs)

  <b>Notes:</b>
   - if no SSSR rings are found on the molecule - MolOps::findSSSR() is called
  first
*/
RDKIT_GRAPHMOL_EXPORT int symmetrizeSSSR(ROMol &mol,
                                         std::vector<std::vector<int>> &res,
                                         bool includeDativeBonds = false);
//! \overload
RDKIT_GRAPHMOL_EXPORT int symmetrizeSSSR(ROMol &mol,
                                         bool includeDativeBonds = false);

//! @}

//! \name Shortest paths and other matrices
//! @{

//! returns a molecule's adjacency matrix
/*!
  \param mol             the molecule of interest
  \param useBO           toggles use of bond orders in the matrix
  \param emptyVal        sets the empty value (for non-adjacent atoms)
  \param force           forces calculation of the matrix, even if already
  computed
  \param propNamePrefix  used to set the cached property name

  \return the adjacency matrix.

  <b>Notes</b>
    - The result of this is cached in the molecule's local property
  dictionary, which will handle deallocation. The caller should <b>not</b> \c
  delete this pointer.

*/
RDKIT_GRAPHMOL_EXPORT double *getAdjacencyMatrix(
    const ROMol &mol, bool useBO = false, int emptyVal = 0, bool force = false,
    const char *propNamePrefix = nullptr,
    const boost::dynamic_bitset<> *bondsToUse = nullptr);

//! Computes the molecule's topological distance matrix
/*!
   Uses the Floyd-Warshall all-pairs-shortest-paths algorithm.

  \param mol             the molecule of interest
  \param useBO           toggles use of bond orders in the matrix
  \param useAtomWts      sets the diagonal elements of the result to
           6.0/(atomic number) so that the matrix can be used to calculate
           Balaban J values.  This does not affect the bond weights.
  \param force           forces calculation of the matrix, even if already
  computed
  \param propNamePrefix  used to set the cached property name

  \return the distance matrix.

  <b>Notes</b>
    - The result of this is cached in the molecule's local property
  dictionary, which will handle deallocation. The caller should <b>not</b> \c
  delete this pointer.


*/
RDKIT_GRAPHMOL_EXPORT double *getDistanceMat(
    const ROMol &mol, bool useBO = false, bool useAtomWts = false,
    bool force = false, const char *propNamePrefix = nullptr);

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
    - The results of this call are not cached, the caller <b>should</b> \c
  delete
      this pointer.


*/
RDKIT_GRAPHMOL_EXPORT double *getDistanceMat(
    const ROMol &mol, const std::vector<int> &activeAtoms,
    const std::vector<const Bond *> &bonds, bool useBO = false,
    bool useAtomWts = false);

//! Computes the molecule's 3D distance matrix
/*!

  \param mol             the molecule of interest
  \param confId          the conformer to use
  \param useAtomWts      sets the diagonal elements of the result to
           6.0/(atomic number)
  \param force           forces calculation of the matrix, even if already
  computed
  \param propNamePrefix  used to set the cached property name
                         (if set to an empty string, the matrix will not be
  cached)

  \return the distance matrix.

  <b>Notes</b>
    - If propNamePrefix is not empty the result of this is cached in the
      molecule's local property dictionary, which will handle deallocation.
      In other cases the caller is responsible for freeing the memory.

*/
RDKIT_GRAPHMOL_EXPORT double *get3DDistanceMat(
    const ROMol &mol, int confId = -1, bool useAtomWts = false,
    bool force = false, const char *propNamePrefix = nullptr);
//! Find the shortest path between two atoms
/*!
  Uses the Bellman-Ford algorithm

 \param mol  molecule of interest
 \param aid1 index of the first atom
 \param aid2 index of the second atom

 \return an std::list with the indices of the atoms along the shortest
    path

 <b>Notes:</b>
   - the starting and end atoms are included in the path
   - if no path is found, an empty path is returned

*/
RDKIT_GRAPHMOL_EXPORT std::list<int> getShortestPath(const ROMol &mol, int aid1,
                                                     int aid2);

//! @}

//! \name Stereochemistry
//! @{

// class to hold hybridizations

class Hybridizations {
 public:
  Hybridizations() {
    throw FileParseException("not to be called without a mol parameter");
  };
  Hybridizations(const ROMol &mol);
  Hybridizations(const Hybridizations &) {
    throw FileParseException("not to be called without a mol parameter");
  };

  ~Hybridizations() = default;

  Atom::HybridizationType operator[](int idx) {
    return static_cast<Atom::HybridizationType>(d_hybridizations[idx]);
  }
  // Atom::HybridizationType &operator[](unsigned int idx) {
  //      return static_cast<Atom::HybridizationType>(d_hybridizations[idx]);
  //   d_hybridizations[d_hybridizations[idx]];
  // }

  // // void clear() { d_hybridizations.clear(); }
  // // void resize(unsigned int sz) { d_hybridizations.resize(sz); }
  unsigned int size() const { return d_hybridizations.size(); }

 private:
  std::vector<int> d_hybridizations;
};

//! removes bogus chirality markers (e.g. tetrahedral flags on non-sp3 centers):
RDKIT_GRAPHMOL_EXPORT void cleanupChirality(RWMol &mol);

//! \overload
RDKIT_GRAPHMOL_EXPORT void cleanupAtropisomers(RWMol &);
//! removes bogus atropisomeric markers (e.g. those without sp2 begin and end
//! atoms):
RDKIT_GRAPHMOL_EXPORT void cleanupAtropisomers(RWMol &mol,
                                               Hybridizations &hybridizations);

//! \brief Uses a conformer to assign ChiralTypes to a molecule's atoms
/*!
  \param mol                  the molecule of interest
  \param confId               the conformer to use
  \param replaceExistingTags  if this flag is true, any existing atomic chiral
                              tags will be replaced

  If the conformer provided is not a 3D conformer, nothing will be done.


  NOTE that this does not check to see if atoms are chiral centers (i.e. all
  substituents are different), it merely sets the chiral type flags based on
  the coordinates and atom ordering. Use \c assignStereochemistryFrom3D() if
  you want chiral flags only on actual stereocenters.
*/
RDKIT_GRAPHMOL_EXPORT void assignChiralTypesFrom3D(
    ROMol &mol, int confId = -1, bool replaceExistingTags = true);

//! \brief Uses a conformer to assign ChiralTypes to a molecule's atoms and
//! stereo flags to its bonds
/*!

  \param mol                  the molecule of interest
  \param confId               the conformer to use
  \param replaceExistingTags  if this flag is true, any existing info about
                              stereochemistry will be replaced

  If the conformer provided is not a 3D conformer, nothing will be done.
*/
RDKIT_GRAPHMOL_EXPORT void assignStereochemistryFrom3D(
    ROMol &mol, int confId = -1, bool replaceExistingTags = true);

//! \brief Use bond directions to assign ChiralTypes to a molecule's atoms and
//! stereo flags to its bonds
/*!

  \param mol                  the molecule of interest
  \param confId               the conformer to use
  \param replaceExistingTags  if this flag is true, any existing info about
                              stereochemistry will be replaced
*/
RDKIT_GRAPHMOL_EXPORT void assignChiralTypesFromBondDirs(
    ROMol &mol, int confId = -1, bool replaceExistingTags = true);

//! \deprecated: this function will be removed in a future release. Use
//! setDoubleBondNeighborDirections() instead
RDKIT_GRAPHMOL_EXPORT void detectBondStereochemistry(ROMol &mol,
                                                     int confId = -1);
//! Sets bond directions based on double bond stereochemistry
RDKIT_GRAPHMOL_EXPORT void setDoubleBondNeighborDirections(
    ROMol &mol, const Conformer *conf = nullptr);
//! removes directions from single bonds. Wiggly bonds will have the property
//! _UnknownStereo set on them
RDKIT_GRAPHMOL_EXPORT void clearSingleBondDirFlags(ROMol &mol,
                                                   bool onlyWedgeFlags = false);

//! removes directions from all bonds. Wiggly bonds and cross bonds will have
//! the property _UnknownStereo set on them
RDKIT_GRAPHMOL_EXPORT void clearAllBondDirFlags(ROMol &mol);
RDKIT_GRAPHMOL_EXPORT void clearDirFlags(ROMol &mol,
                                         bool onlyWedgeFlags = false);

//! Assign CIS/TRANS bond stereochemistry tags based on neighboring
//! directions
RDKIT_GRAPHMOL_EXPORT void setBondStereoFromDirections(ROMol &mol);

//! Assign stereochemistry tags to atoms and bonds.
/*!
  If useLegacyStereoPerception is true, it also does the CIP stereochemistry
  assignment for the molecule's atoms (R/S) and double bonds (Z/E).
  This assignment is based on legacy code which is fast, but is
  known to incorrectly assign CIP labels in some cases.
  instead, to assign CIP labels based on an accurate, though slower,
  implementation of the CIP rules, call CIPLabeler::assignCIPLabels().
  Chiral atoms will have a property '_CIPCode' indicating their chiral code.

  \param mol     the molecule to use
  \param cleanIt if true, any existing values of the property `_CIPCode`
                 will be cleared, atoms with a chiral specifier that aren't
                 actually chiral (e.g. atoms with duplicate
                 substituents or only 2 substituents, etc.) will have
                 their chiral code set to CHI_UNSPECIFIED. Bonds with
                 STEREOCIS/STEREOTRANS specified that have duplicate
                 substituents based upon the CIP atom ranks will be
                 marked STEREONONE.
  \param force   causes the calculation to be repeated even if it has
                 already been done
  \param flagPossibleStereoCenters   set the _ChiralityPossible property on
                                     atoms that are possible stereocenters

  <b>Notes:M</b>
    - Throughout we assume that we're working with a hydrogen-suppressed
      graph.

*/
RDKIT_GRAPHMOL_EXPORT void assignStereochemistry(
    ROMol &mol, bool cleanIt = false, bool force = false,
    bool flagPossibleStereoCenters = false);
//! Removes all stereochemistry information from atoms (i.e. R/S) and bonds
/// i.e. Z/E)
/*!

  \param mol     the molecule of interest
*/
RDKIT_GRAPHMOL_EXPORT void removeStereochemistry(ROMol &mol);

//! \brief finds bonds that could be cis/trans in a molecule and mark them as
//!  Bond::STEREOANY.
/*!
  \param mol     the molecule of interest
  \param cleanIt toggles removal of stereo flags from double bonds that can
                 not have stereochemistry

  This function finds any double bonds that can potentially be part of
  a cis/trans system. No attempt is made here to mark them cis or
  trans. No attempt is made to detect double bond stereo in ring systems.

  This function is useful in the following situations:
    - when parsing a mol file; for the bonds marked here, coordinate
      information on the neighbors can be used to indentify cis or trans
  states
    - when writing a mol file; bonds that can be cis/trans but not marked as
      either need to be specially marked in the mol file
    - finding double bonds with unspecified stereochemistry so they
      can be enumerated for downstream 3D tools

  The CIPranks on the neighboring atoms are checked in this function. The
  _CIPCode property if set to any on the double bond.
*/
RDKIT_GRAPHMOL_EXPORT void findPotentialStereoBonds(ROMol &mol,
                                                    bool cleanIt = false);
//! \brief Uses the molParity atom property to assign ChiralType to a
//! molecule's atoms
/*!
  \param mol                  the molecule of interest
  \param replaceExistingTags  if this flag is true, any existing atomic chiral
                              tags will be replaced
*/
RDKIT_GRAPHMOL_EXPORT void assignChiralTypesFromMolParity(
    ROMol &mol, bool replaceExistingTags = true);

//! @}

//! returns the number of atoms which have a particular property set
RDKIT_GRAPHMOL_EXPORT unsigned getNumAtomsWithDistinctProperty(
    const ROMol &mol, std::string prop);

//! returns whether or not a molecule needs to have Hs added to it.
RDKIT_GRAPHMOL_EXPORT bool needsHs(const ROMol &mol);

//! \brief Replaces haptic bond with explicit dative bonds.
/*!
 *
 * @param mol the molecule of interest
 *
 * One way of showing haptic bonds (such as cyclopentadiene to iron in
 * ferrocene) is to use a dummy atom with a dative bond to the iron atom with
 * the bond labelled with the atoms involved in the organic end of the bond.
 * Another way is to have explicit dative bonds from the atoms of the haptic
 * group to the metal atom.  This function converts the former representation to
 * the latter.
 */
RDKIT_GRAPHMOL_EXPORT ROMol *hapticBondsToDative(const ROMol &mol);

//! \overload modifies molecule in place.
RDKIT_GRAPHMOL_EXPORT void hapticBondsToDative(RWMol &mol);

//! \brief Replaces explicit dative bonds with haptic.
/*!
 *
 * @param mol the molecule of interest
 *
 * Does the reverse of hapticBondsToDative.  If there are multiple contiguous
 * atoms attached by dative bonds to an atom (probably a metal atom), the dative
 * bonds will be replaced by a dummy atom in their centre attached to the
 * (metal) atom by a dative bond, which is labelled with ENDPTS of the atoms
 * that had the original dative bonds.
 */
RDKIT_GRAPHMOL_EXPORT ROMol *dativeBondsToHaptic(const ROMol &mol);

//! \overload modifies molecule in place.
RDKIT_GRAPHMOL_EXPORT void dativeBondsToHaptic(RWMol &mol);

namespace details {
//! not recommended for use in other code
RDKIT_GRAPHMOL_EXPORT void KekulizeFragment(
    RWMol &mol, const boost::dynamic_bitset<> &atomsToUse,
    boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds = true,
    unsigned int maxBackTracks = 100);

// If the bond is dative, and it has a common_properties::MolFileBondEndPts
// prop, returns a vector of the indices of the atoms mentioned in the prop.
RDKIT_GRAPHMOL_EXPORT std::vector<int> hapticBondEndpoints(const Bond *bond);

}  // namespace details

//! attachment points encoded as attachPt properties are added to the graph as
/// dummy atoms
/*!
 *
 * @param mol the molecule of interest
 * @param addAsQueries if true, the dummy atoms will be added as null queries
 *       (i.e. they will match any atom in a substructure search)
 * @param addCoords if true and the molecule has one or more conformers,
 *    positions for the attachment points will be added to the conformer(s).
 *
 */
RDKIT_GRAPHMOL_EXPORT void expandAttachmentPoints(RWMol &mol,
                                                  bool addAsQueries = true,
                                                  bool addCoords = true);
//! dummy atoms in the graph are removed and replaced with attachment point
//! annotations on the attached atoms
/*!
 *
 * @param mol the molecule of interest
 * @param markedOnly if true, only dummy atoms with the _fromAttachPoint
 *    property will be collapsed
 *
 * In order for a dummy atom to be considered for collapsing it must have:
 * - degree 1 with a single or unspecified bond
 * - the bond to it can not be wedged
 * - either no query or be an AtomNullQuery
 *
 */
RDKIT_GRAPHMOL_EXPORT void collapseAttachmentPoints(RWMol &mol,
                                                    bool markedOnly = true);

namespace details {
//! attachment points encoded as attachPt properties are added to the graph as
/// dummy atoms
/*!
 *
 * @param mol the molecule of interest
 * @param atomIdx the index of the atom to which the attachment point should be
 *       added
 * @param val the attachment point value. Should be 1 or 2
 * @param addAsQueries if true, the dummy atoms will be added as null queries
 *       (i.e. they will match any atom in a substructure search)
 * @param addCoords if true and the molecule has one or more conformers,
 *    positions for the attachment points will be added to the conformer(s).
 *
 */
RDKIT_GRAPHMOL_EXPORT unsigned int addExplicitAttachmentPoint(
    RWMol &mol, unsigned int atomIdx, unsigned int val, bool addAsQuery = true,
    bool addCoords = true);

//! returns whether or not an atom is an attachment point
/*!
 *
 * @param mol the molecule of interest
 * @param markedOnly if true, only dummy atoms with the _fromAttachPoint
 *    property will be collapsed
 *
 * In order for a dummy atom to be considered for collapsing it must have:
 * - degree 1 with a single or unspecified bond
 * - the bond to it can not be wedged
 * - either no query or be an AtomNullQuery
 *
 */
RDKIT_GRAPHMOL_EXPORT bool isAttachmentPoint(const Atom *atom,
                                             bool markedOnly = true);

}  // namespace details

}  // namespace MolOps
}  // namespace RDKit

#endif
