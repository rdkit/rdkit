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

%typemap(javaimports) RDKit::ROMol "
/** 
ROMol is a molecule class that is intended to have a fixed topology.
<p>
This is the primary class for most molecule operations.
<p>
If you need to be manipulating the molecule (e.g. adding or deleting atoms or bonds, use an RWMol instead.
<p>
<p>
@notes
<li>each ROMol maintains a Dict of properties:
<li>o Each property is keyed by name and can store an arbitrary type.
<li>o Properties can be marked as calculated, in which case they will be cleared when the clearComputedProps() method is called.
<li>o Because they have no impact upon chemistry, all property operations are const, this allows extra flexibility for clients who need to store extra data on ROMol objects.
<li>each ROMol has collections of bookmarks for Atoms and Bonds:
<li>o the Atom bookmarks and Bond bookmarks are stored separately from each other
<li>o each bookmark, an integer, can map to more than one Atom or Bond
<li>o these are currently used in molecule construction, but could also be useful for reaction mapping and the like
<li>information about rings (SSSR and the like) is stored in the molecule's RingInfo pointer.
 */"

%javamethodmodifiers RDKit::ROMol::addConformer 	( 	Conformer *  	conf, 		bool  	assignId = false	  	) 			"
/**
<p>
Add a new conformation to the molecule.
<p>
<p>
@param
conf 	- conformation to be added to the molecule, this molecule takes ownership of the conformer
assignId 	- a unique ID will be assigned to the the conformation if true otherwise it is assumed that the conformation already has an (unique) ID set

*/
public";

%javamethodmodifiers RDKit::ROMol::beginAromaticAtoms 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::beginAtoms 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::beginBonds 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::beginHeteros 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::beginQueryAtoms 	( 	QueryAtom const *  	what 	 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::clearAtomBookmark 	( 	const int  	mark, 		ATOM_SPTR  	atom	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::clearBondBookmark 	( 	int  	mark, 		BOND_SPTR  	bond	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::clearProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::clearProp 	( 	const char *  	key 	 )  	const "
/**
<p>
clears the value of a property
<p>
<p>
@notes
<li>if no property with name key exists, a KeyErrorException will be thrown.
<li>if the property is marked as computed, it will also be removed from our list of computedProperties

*/
public";

%javamethodmodifiers RDKit::ROMol::endAromaticAtoms 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::endAtoms 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::endBonds 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::endHeteros 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::endQueryAtoms 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getAtomBonds 	( 	Atom const *  	at 	 )  	const"
/**
<p>
provides access to all Bond objects connected to an Atom
<p>
<p>
@param
at 	the atom whose neighbors we are looking for
<p>
@example
<pre><code>
... molPtr is a const ROMol * ...
... atomPtr is a const Atom * ...
ROMol::OEDGE_ITER beg,end;
boost::tie(beg,end) = molPtr->getAtomBonds(atomPtr);
while(beg!=end){
const BOND_SPTR bond=(*molPtr)[*beg];
... do something with the Bond ...
++beg;
}
</code></pre>
or, if you need a non-const Bond *:
<p>
        ... molPtr is a ROMol * ...
        ... atomPtr is a const Atom * ...
        ROMol::OEDGE_ITER beg,end;
        boost::tie(beg,end) = molPtr->getAtomBonds(atomPtr);
        while(beg!=end){
          BOND_SPTR bond=(*molPtr)[*beg];
          ... do something with the Bond ...
          ++beg;
        }

*/
public";

%javamethodmodifiers RDKit::ROMol::getAtomDegree 	( 	Atom::ATOM_SPTR  	at 	 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getAtomNeighbors 	( 	Atom::ATOM_SPTR  	at 	 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getAtomNeighbors 	( 	Atom const *  	at 	 )  	const"
/**
<p>
provides access to all neighbors around an Atom
<p>
<p>
@param
at 	the atom whose neighbors we are looking for
<p>
@example
<pre><code>
... molPtr is a const ROMol & ...
... atomPtr is a const Atom * ...
ROMol::ADJ_ITER nbrIdx,endNbrs;
boost::tie(nbrIdx,endNbrs) = molPtr.getAtomNeighbors(atomPtr);
while(nbrIdx!=endNbrs){
const ATOM_SPTR at=molPtr[*nbrIdx];
... do something with the Atom ...
++nbrIdx;
}
</code></pre>

*/
public";

%javamethodmodifiers RDKit::ROMol::getAtomWithIdx 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getBondBetweenAtoms 	( 	unsigned int  	idx1, 		unsigned int  	idx2	  	) 			const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getBondWithIdx 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getConformer 	( 	int  	id = -1 	 )  	"
/**
<p>
return the conformer with a specified ID if the ID is negative the first conformation will be returned
*/
public";

%javamethodmodifiers RDKit::ROMol::getConformer 	( 	int  	id = -1 	 )  	const"
/**
<p>
return the conformer with a specified ID if the ID is negative the first conformation will be returned
*/
public";

%javamethodmodifiers RDKit::ROMol::getEdges 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getEdges 	( 		 )  	"
/**
<p>
<p>
@return
an iterator pair for looping over all Bonds
<p>
@example
<pre><code>
ROMol::EDGE_ITER firstB,lastB;
boost::tie(firstB,lastB) = mol.getEdges();
while(firstB!=lastB){
BOND_SPTR bond = mol[*firstB];
... do something with the Bond ...
++firstB;
}
</code></pre>
template<typename T >
*/
public";

%javamethodmodifiers RDKit::ROMol::getProp 	( 	const std::string  	key, 		T &  	res	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::ROMol::getProp 	( 	const char *  	key, 		T &  	res	  	) 			const "
/**
<p>
allows retrieval of a particular property value
<p>
<p>
@param
key 	the name under which the property should be stored. If a property is already stored under this name, it will be replaced.
res 	a reference to the storage location for the value.
<p>
@notes
<li>if no property with name key exists, a KeyErrorException will be thrown.
<li>the boost::lexical_cast machinery is used to attempt type conversions. If this fails, a boost::bad_lexical_cast exception will be thrown.

*/
public";

%javamethodmodifiers RDKit::ROMol::getRingInfo 	( 		 )  	const "
/**
<p>
<p>
@return
a pointer to our RingInfo structure Note: the client should not delete this.
.
*/
public";

%javamethodmodifiers RDKit::ROMol::getTopology 	( 		 )  	const "
/**
<p>
brief returns a pointer to our underlying BGL object
<p>
This can be useful if you need to call other BGL algorithms:
<p>
Here's an example:
<p>
           ... mol is a const ROMol ...
           ... mapping is an INT_VECT ...
           mapping.resize(mol.getNumAtoms());
           const MolGraph &G_p = mol.getTopology();
           int res = boost::connected_components(G_p,&mapping[0]);

*/
public";

%javamethodmodifiers RDKit::ROMol::getVertices 	( 		 )  	const"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::getVertices 	( 		 )  	"
/**
<p>
<p>
@return
an iterator pair for looping over all Atoms
<p>
@example
<pre><code>
ROMol::VERTEX_ITER atBegin,atEnd;
boost::tie(atBegin,atEnd) = mol.getVertices();
while(atBegin!=atEnd){
ATOM_SPTR at2=mol[*atBegin];
... do something with the Atom ...
++atBegin;
}
</code></pre>

*/
public";

%javamethodmodifiers RDKit::ROMol::hasProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::setAtomBookmark 	( 	Atom *  	at, 		int  	mark	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::ROMol::setBondBookmark 	( 	Bond *  	bond, 		int  	mark	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::ROMol::setProp 	( 	const std::string  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::ROMol::setProp 	( 	const char *  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
/**
<p>
sets a property value
<p>
<p>
@param
key 	the name under which the property should be stored. If a property is already stored under this name, it will be replaced.
val 	the value to be stored
computed 	(optional) allows the property to be flagged computed.

*/
public";

%javamethodmodifiers RDKit::ROMol::updatePropertyCache 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates any of our lazy properties
<p>
<p>
@notes
<li>this calls updatePropertyCache() on each of our Atoms and Bonds

*/
public";

