/* 
* $Id$
*
*  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
*  All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met: 
*
*     * Redistributions of source code must retain the above copyright 
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above
*       copyright notice, this list of conditions and the following 
*       disclaimer in the documentation and/or other materials provided 
*       with the distribution.
*     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
*       nor the names of its contributors may be used to endorse or promote 
*       products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

%typemap(javaimports) RDKit::MolOps "
/** 
Groups a variety of molecular query and transformation operations. */"

%javamethodmodifiers RDKit::MolOps::addHs 	( 	const ROMol &  	mol, 		bool  	explicitOnly = false, 		bool  	addCoords = false	  	) 			"
/**
<p>
returns a copy of a molecule with hydrogens added in as explicit Atoms
<p>
<p>
@param
mol 	the molecule to add Hs to
explicitOnly 	(optional) if this true, only explicit Hs will be added
addCoords 	(optional) If this is true, estimates for the atomic coordinates of the added Hs will be used.
<p>
@return
the new molecule
<p>
@notes
<li>it makes no sense to use the addCoords option if the molecule's heavy atoms don't already have coordinates.
<li>the caller is responsible for deleteing the pointer this returns.

*/
public";

%javamethodmodifiers RDKit::MolOps::adjustHs 	( 	RWMol &  	mol 	 )  	"
/**
<p>
adjust the number of implicit and explicit Hs for special cases
<p>
Currently this:
<p>
    * modifies aromatic nitrogens so that, when appropriate, they have an explicit H marked (e.g. so that we get things like 'c1cc[nH]cc1'
<p>
<p>
@param
mol 	the molecule of interest
Assumptions
<p>
    * this is called after the molecule has been sanitized, aromaticity has been perceived, and the implicit valence of everything has been calculated.

*/
public";

%javamethodmodifiers RDKit::MolOps::assignChiralTypesFrom3D 	( 	ROMol &  	mol, 		int  	confId = -1, 		bool  	replaceExistingTags = true	  	) 			"
/**
<p>
Uses a conformer to assign ChiralType to a molecule's atoms.
<p>
<p>
@param
mol 	the molecule of interest
confId 	the conformer to use
replaceExistingTags 	if this flag is true, any existing atomic chiral tags will be replaced
If the conformer provided is not a 3D conformer, nothing will be done.
*/
public";

%javamethodmodifiers RDKit::MolOps::assignStereochemistry 	( 	ROMol &  	mol, 		bool  	cleanIt = false, 		bool  	force = false	  	) 			"
/**
<p>
Assign stereochemistry tags to atoms (i.e. R/S) and bonds (i.e. Z/E).
<p>
<p>
@param
mol 	the molecule of interest
cleanIt 	toggles removal of stereo flags from double bonds that can not have stereochemistry
force 	forces the calculation to be repeated even if it has already been done
<p>
@notes
<li>Throughout we assume that we're working with a hydrogen-suppressed graph.

*/
public";

%javamethodmodifiers RDKit::MolOps::cleanUp 	( 	RWMol &  	mol 	 )  	"
/**
<p>
Designed to be called by the sanitizer to handle special cases before anything is done.
<p>
Currently this:
<p>
    * modifies nitro groups, so that the nitrogen does not have a unreasonable valence of 5, as follows:
          o the nitrogen gets a positve charge
          o one of the oxygens gets a negative chage and the double bond to this oxygen is changed to a single bond The net result is that nitro groups can be counted on to be: '[N+](=O)[O-]'
<p>
<p>
@param
mol 	the molecule of interest

*/
public";

%javamethodmodifiers RDKit::MolOps::computeBalabanJ 	( 	const ROMol &  	mol, 		bool  	useBO = true, 		bool  	force = false, 		const std::vector< int > *  	bondPath = 0, 		bool  	cacheIt = true	  	) 			"
/**
<p>
calculates Balaban's J index for the molecule
<p>
<p>
@param
mol 	the molecule of interest
useBO 	toggles inclusion of the bond order in the calculation (when false, we're not really calculating the J value)
force 	forces the calculation (instead of using cached results)
bondPath 	when included, only paths using bonds whose indices occur in this vector will be included in the calculation
cacheIt 	If this is true, the calculated value will be cached as a property on the molecule
<p>
@return
the J index

*/
public";

%javamethodmodifiers RDKit::MolOps::computeBalabanJ 	( 	double *  	distMat, 		int  	nb, 		int  	nAts	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::MolOps::countAtomElec 	( 	const Atom *  	at 	 )  	"
/**
<p>
return the number of electrons available on an atom to donate for aromaticity
<p>
The result is determined using the default valency, number of lone pairs, number of bonds and the formal charge. Note that the atom may not donate all of these electrons to a ring for aromaticity (also used in Conjugation and hybridization code).
<p>
<p>
@param
at 	the atom of interest
<p>
@return
the number of electrons

*/
public";

%javamethodmodifiers RDKit::MolOps::findPotentialStereoBonds 	( 	ROMol &  	mol, 		bool  	cleanIt = false	  	) 			"
/**
<p>
finds bonds that could be cis/trans in a molecule and mark them as Bond::STEREONONE
<p>
<p>
@param
mol 	the molecule of interest
cleanIt 	toggles removal of stereo flags from double bonds that can not have stereochemistry
This function is usefuly in two situations
<p>
    * when parsing a mol file; for the bonds marked here, coordinate informations on the neighbors can be used to indentify cis or trans states
    * when writing a mol file; bonds that can be cis/trans but not marked as either need to be specially marked in the mol file

*/
public";

%javamethodmodifiers RDKit::MolOps::findSSSR 	( 	const ROMol &  	mol, 		VECT_INT_VECT &  	res	  	) 			"
/**
<p>
finds a molecule's Smallest Set of Smallest Rings
<p>
Currently this implements a modified form of Figueras algorithm (JCICS - Vol. 36, No. 5, 1996, 986-991)
<p>
<p>
@param
mol 	the molecule of interest
res 	used to return the vector of rings. Each entry is a vector with atom indices. This information is also stored in the molecule's RingInfo structure, so this argument is optional (see overload)
<p>
@return
number of smallest rings found
Base algorithm:
<p>
    * The original algorithm starts by finding representative degree 2 nodes.
    * Representative because if a series of deg 2 nodes are found only one of them is picked.
    * The smallest ring around each of them is found.
    * The bonds that connect to this degree 2 node are them chopped off, yielding new deg two nodes
    * The process is repeated on the new deg 2 nodes.
    * If no deg 2 nodes are found, a deg 3 node is picked. The smallest ring with it is found. A bond from this is 'carefully' (look in the paper) selected and chopped, yielding deg 2 nodes. The process is same as above once this is done.
<p>
Our Modifications:
<p>
    * If available, more than one smallest ring around a representative deg 2 node will be computed and stored
    * Typically 3 rings are found around a degree 3 node (when no deg 2s are available) and all the bond to that node are chopped.
    * The extra rings that were found in this process are removed after all the nodes have been covered.
<p>
These changes were motivated by several factors:
<p>
    * We believe the original algorithm fails to find the correct SSSR (finds the correct number of them but the wrong ones) on some sample mols
    * Since SSSR may not be unique, a post-SSSR step to symmetrize may be done. The extra rings this process adds can be quite useful.

*/
public";

%javamethodmodifiers RDKit::MolOps::findSSSR 	( 	const ROMol &  	mol, 		VECT_INT_VECT *  	res = 0	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::MolOps::getAdjacencyMatrix 	( 	const ROMol &  	mol, 		bool  	useBO = false, 		int  	emptyVal = 0, 		bool  	force = false, 		const char *  	propNamePrefix = 0	  	) 			"
/**
<p>
returns a molecule's adjacency matrix
<p>
<p>
@param
mol 	the molecule of interest
useBO 	toggles use of bond orders in the matrix
emptyVal 	sets the empty value (for non-adjacent atoms)
force 	forces calculation of the matrix, even if already computed
propNamePrefix 	used to set the cached property name
<p>
@return
the adjacency matrix.
<p>
@notes
<li>The result of this is cached in the molecule's local property dictionary, which will handle deallocation. Do the caller should not delete this pointer.

*/
public";

%javamethodmodifiers RDKit::MolOps::getDistanceMat 	( 	const ROMol &  	mol, 		const std::vector< int > &  	activeAtoms, 		const std::vector< const Bond * > &  	bonds, 		bool  	useBO = false, 		bool  	useAtomWts = false	  	) 			"
/**
<p>
Computes the molecule's topological distance matrix.
<p>
Uses the Floyd-Warshall all-pairs-shortest-paths algorithm.
<p>
<p>
@param
mol 	the molecule of interest
activeAtoms 	only elements corresponding to these atom indices will be included in the calculation
bonds 	only bonds found in this list will be included in the calculation
useBO 	toggles use of bond orders in the matrix
useAtomWts 	sets the diagonal elements of the result to 6.0/(atomic number) so that the matrix can be used to calculate Balaban J values. This does not affect the bond weights.
<p>
@return
the distance matrix.
<p>
@notes
<li>The results of this call are not cached, the caller should delete this pointer.

*/
public";

%javamethodmodifiers RDKit::MolOps::getDistanceMat 	( 	const ROMol &  	mol, 		bool  	useBO = false, 		bool  	useAtomWts = false, 		bool  	force = false, 		const char *  	propNamePrefix = 0	  	) 			"
/**
<p>
Computes the molecule's topological distance matrix.
<p>
Uses the Floyd-Warshall all-pairs-shortest-paths algorithm.
<p>
<p>
@param
mol 	the molecule of interest
useBO 	toggles use of bond orders in the matrix
useAtomWts 	sets the diagonal elements of the result to 6.0/(atomic number) so that the matrix can be used to calculate Balaban J values. This does not affect the bond weights.
force 	forces calculation of the matrix, even if already computed
propNamePrefix 	used to set the cached property name
<p>
@return
the distance matrix.
<p>
@notes
<li>The result of this is cached in the molecule's local property dictionary, which will handle deallocation. Do the caller should not delete this pointer.

*/
public";

%javamethodmodifiers RDKit::MolOps::getMolFrags 	( 	const ROMol &  	mol, 		VECT_INT_VECT &  	frags	  	) 			"
/**
<p>
find fragments (disconnected components of the molecular graph)
<p>
<p>
@param
mol 	the molecule of interest
frags 	used to return the Atoms in each fragment On return mapping will be numFrags long, and each entry will contain the indices of the Atoms in that fragment.
<p>
@return
the number of fragments found.

*/
public";

%javamethodmodifiers RDKit::MolOps::getMolFrags 	( 	const ROMol &  	mol, 		INT_VECT &  	mapping	  	) 			"
/**
<p>
find fragments (disconnected components of the molecular graph)
<p>
<p>
@param
mol 	the molecule of interest
mapping 	used to return the mapping of Atoms->fragments. On return mapping will be mol->getNumAtoms() long and will contain the fragment assignment for each Atom
<p>
@return
the number of fragments found.

*/
public";

%javamethodmodifiers RDKit::MolOps::getMolFrags 	( 	const ROMol &  	mol, 		bool  	sanitizeFrags = true, 		INT_VECT *  	frags = 0	  	) 			"
/**
<p>
splits a molecule into its component fragments
<p>
<p>
@param
mol 	the molecule of interest
sanitizeFrags 	toggles sanitization of the fragments after they are built
frags 	used to return the mapping of Atoms->fragments. if provided, frags will be mol->getNumAtoms() long on return and will contain the fragment assignment for each Atom
<p>
@return
a vector of the fragments as smart pointers to ROMols

*/
public";

%javamethodmodifiers RDKit::MolOps::getShortestPath 	( 	const ROMol &  	mol, 		int  	aid1, 		int  	aid2	  	) 			"
/**
<p>
Find the shortest path between two atoms.
<p>
Uses the Bellman-Ford algorithm
<p>
<p>
@param
mol 	molecule of interest
aid1 	index of the first atom
aid2 	index of the second atom
<p>
@return
an std::list with the indices of the atoms along the shortest path
<p>
@notes
<li>the starting and end atoms are included in the path
<li>if no path is found, an empty path is returned

*/
public";

%javamethodmodifiers RDKit::MolOps::Kekulize 	( 	RWMol &  	mol, 		bool  	markAtomsBonds = true, 		unsigned int  	maxBackTracks = 100	  	) 			"
/**
<p>
Kekulizes the molecule.
<p>
<p>
@param
mol 	the molecule of interest
markAtomsBonds 	if this is set to true, isAromatic boolean settings on both the Bonds and Atoms are turned to false following the Kekulization, otherwise they are left alone in their original state.
maxBackTracks 	the maximum number of attempts at back-tracking. The algorithm uses a back-tracking procedure to revist a previous setting of double bond if we hit a wall in the kekulization process
<p>
@notes
<li>even if markAtomsBonds is false the BondType for all aromatic bonds will be changed from RDKit::Bond::AROMATIC to RDKit::Bond::SINGLE or RDKit::Bond::DOUBLE during Kekulization.

*/
public";

%javamethodmodifiers RDKit::MolOps::mergeQueryHs 	( 	const ROMol &  	mol 	 )  	"
/**
<p>
returns a copy of a molecule with hydrogens removed and added as queries to the heavy atoms to which they are bound.
<p>
This is really intended to be used with molecules that contain QueryAtoms
<p>
<p>
@param
mol 	the molecule to remove Hs from
<p>
@return
the new molecule
<p>
@notes
<li>Atoms that do not already have hydrogen count queries will have one added, other H-related queries will not be touched. Examples:
<li>o C[H] -> [C;!H0]
<li>o [C;H1][H] -> [C;H1]
<li>o [C;H2][H] -> [C;H2]
<li>Hydrogens which aren't connected to a heavy atom will not be removed. This prevents molecules like '[H][H]' from having all atoms removed.
<li>o the caller is responsible for deleteing the pointer this returns.

*/
public";

%javamethodmodifiers RDKit::MolOps::rankAtoms 	( 	const ROMol &  	mol, 		INT_VECT &  	ranks, 		bool  	breakTies = true, 		VECT_INT_VECT *  	rankHistory = 0	  	) 			"
/**
<p>
assign a canonical ordering to a molecule's atoms
<p>
The algorithm used here is a modification of the published Daylight canonical smiles algorithm (i.e. it uses atom invariants and products of primes).
<p>
<p>
@param
mol 	the molecule of interest
ranks 	used to return the ranks
breakTies 	toggles breaking of ties (see below)
rankHistory 	used to return the rank history (see below)
<p>
@notes
<li>Tie breaking should be done when it's important to have a full ordering of the atoms (e.g. when generating canonical traversal trees). If it's acceptable to have ties between symmetry-equivalent atoms (e.g. when generating CIP codes), tie breaking can/should be skipped.
<li>o if the rankHistory argument is provided, the evolution of the ranks of individual atoms will be tracked. The rankHistory pointer should be to a VECT_INT_VECT that has at least mol.getNumAtoms() elements.

*/
public";

%javamethodmodifiers RDKit::MolOps::removeHs 	( 	const ROMol &  	mol, 		bool  	implicitOnly = false, 		bool  	updateExplicitCount = false, 		bool  	sanitize = true	  	) 			"
/**
<p>
returns a copy of a molecule with hydrogens removed
<p>
<p>
@param
mol 	the molecule to remove Hs from
implicitOnly 	(optional) if this true, only implicit Hs will be removed
updateExplicitCount 	(optional) If this is true, when explicit Hs are removed from the graph, the heavy atom to which they are bound will have its counter of explicit Hs increased.
sanitize,: 	(optional) If this is true, the final molecule will be sanitized
<p>
@return
the new molecule
<p>
@notes
<li>Hydrogens which aren't connected to a heavy atom will not be removed. This prevents molecules like '[H][H]' from having all atoms removed.
<li>Labelled hydrogen (e.g. atoms with atomic number=1, but mass > 1), will not be removed.
<li>the caller is responsible for deleteing the pointer this returns.

*/
public";

%javamethodmodifiers RDKit::MolOps::removeStereochemistry 	( 	ROMol &  	mol 	 )  	"
/**
<p>
Removes all stereochemistry information from atoms (i.e. R/S) and bonds (i.e. Z/E).
<p>
<p>
@param
mol 	the molecule of interest

*/
public";

%javamethodmodifiers RDKit::MolOps::sanitizeMol 	( 	RWMol &  	mol 	 )  	"
/**
<p>
carries out a collection of tasks for cleaning up a molecule and ensuring that it makes 'chemical sense'
<p>
This functions calls the following in sequence
<p>
   1. MolOps::cleanUp()
   2. MolOps::Kekulize()
   3. MolOps::setAromaticity()
   4. MolOps::setConjugation()
   5. MolOps::setHybridization()
   6. MolOps::cleanupChirality()
   7. MolOps::adjustHs()
<p>
<p>
@param
mol 	the RWMol to be cleaned
<p>
@notes
<li>If there is a failure in the sanitization, a SanitException will be thrown.
<li>in general the user of this function should cast the molecule following this function to a ROMol, so that new atoms and bonds cannot be added to the molecule and screw up the sanitizing that has been done here

*/
public";

%javamethodmodifiers RDKit::MolOps::setAromaticity 	( 	RWMol &  	mol 	 )  	"
/**
<p>
Sets up the aromaticity for a molecule.
<p>
This is what happens here:
<p>
   1. find all the simple rings by calling the findSSSR function
   2. loop over all the Atoms in each ring and mark them if they are candidates for aromaticity. A ring atom is a candidate if it can spare electrons to the ring and if it's from the first two rows of the periodic table.
   3. ased on the candidate atoms, mark the rings to be either candidates or non-candidates. A ring is a candidate only if all its atoms are candidates
   4. apply Hueckel rule to each of the candidate rings to check if the ring can be aromatic
<p>
<p>
@param
mol 	the RWMol of interest
<p>
@return
1 on succes, 0 otherwise
Assumptions:
<p>
    * Kekulization has been done (i.e. MolOps::Kekulize() has already been called)

*/
public";

%javamethodmodifiers RDKit::MolOps::symmetrizeSSSR 	( 	ROMol &  	mol, 		VECT_INT_VECT &  	res	  	) 			"
/**
<p>
symmetrize the molecule's Smallest Set of Smallest Rings
<p>
SSSR rings obatined from 'findSSSR' can be non-unique in some case. For example, cubane has five SSSR rings, not six as one would hope.
<p>
This function adds additional rings to the SSSR list if necessary to make the list symmetric, e.g. all atoms in cubane will be part of the same number of SSSRs. This function choses these extra rings from the extra rings computed and discarded during findSSSR. The new ring are chosen such that:
<p>
    * replacing a same sized ring in the SSSR list with an extra ring yields the same union of bond IDs as the orignal SSSR list
<p>
<p>
@param
mol 	- the molecule of interest
res 	used to return the vector of rings. Each entry is a vector with atom indices. This information is also stored in the molecule's RingInfo structure, so this argument is optional (see overload)
<p>
@return
the total number of rings = (new rings + old SSSRs)
<p>
@notes
<li>if no SSSR rings are found on the molecule - MolOps::findSSSR() is called first

*/
public";

%javamethodmodifiers RDKit::MolOps::symmetrizeSSSR 	( 	ROMol &  	mol 	 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>

*/
public";

