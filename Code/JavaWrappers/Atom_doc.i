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

%typemap(javaimports) RDKit::Atom "
/** 
The class for representing atoms.
<p>
<p>
@notes
<li>many of the methods of Atom require that the Atom be associated with a molecule (an ROMol).
<li>each Atom maintains a Dict of properties:
<li>o Each property is keyed by name and can store an arbitrary type.
<li>o Properties can be marked as calculated, in which case they will be cleared when the clearComputedProps() method is called.
<li>o Because they have no impact upon chemistry, all property operations are const, this allows extra flexibility for clients who need to store extra data on Atom objects.
<li>Atom objects are lazy about computing their explicit and implicit valence values. These will not be computed until their values are requested.
Chirality:
<p>
The chirality of an Atom is determined by two things:
<p>
    * its chiralTag
    * the input order of its bonds (see note below for handling of implicit Hs)
<p>
For tetrahedral coordination, the chiralTag tells you what direction you have to rotate to get from bond 2 to bond 3 while looking down bond 1. This is pretty much identical to the SMILES representation of chirality.
<p>
NOTE: if an atom has an implicit H, the bond to that H is considered to be at the *end* of the list of other bonds. */"

%javamethodmodifiers RDKit::Atom::calcExplicitValence 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates and returns our explicit valence
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::calcImplicitValence 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates and returns our implicit valence
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::clearProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Atom::clearProp 	( 	const char *  	key 	 )  	const "
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

%javamethodmodifiers RDKit::Atom::copy 	( 		 )  	const "
/**
<p>
makes a copy of this Atom and returns a pointer to it.
<p>
<p>
@notes
Reimplemented in RDKit::QueryAtom.
*/
public";

%javamethodmodifiers RDKit::Atom::getDegree 	( 		 )  	const"
/**
<p>
<p>
@return
the explicit degree of the Atom (number of bonded neighbors in the graph)
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::getImplicitValence 	( 		 )  	const"
/**
<p>
<p>
@return
the implicit valence for this Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::getNumImplicitHs 	( 		 )  	const"
/**
<p>
<p>
@return
the number of implicit Hs this Atom is bound to
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::getNumRadicalElectrons 	( 		 )  	const "
/**
<p>
<p>
@return
the number of radical electrons for this Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::getPerturbationOrder 	( 	INT_LIST  	probe 	 )  	const"
/**
<p>
returns the perturbation order for a list of integers
<p>
This value is associated with chirality.
<p>
<p>
@param
probe 	a list of bond indices. This must be the same length as our number of incoming bonds (our degree).
<p>
@return
the number of swaps required to convert the ordering of the probe list to match the order of our incoming bonds: e.g. if our incoming bond order is: [0,1,2,3]
getPerturbationOrder([1,0,2,3]) = 1
getPerturbationOrder([1,2,3,0]) = 3
getPerturbationOrder([1,2,0,3]) = 2
See the class documentation for a more detailed description of our representation of chirality.
<p>
<p>
@notes
<li>requires an owning molecule
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Atom::getProp 	( 	const std::string  	key, 		T &  	res	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Atom::getProp 	( 	const char *  	key, 		T &  	res	  	) 			const "
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

%javamethodmodifiers RDKit::Atom::getTotalDegree 	( 		 )  	const"
/**
<p>
<p>
@return
the total degree of the Atom (number of bonded neighbors + number of Hs)
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::getTotalNumHs 	( 	bool  	includeNeighbors = false 	 )  	const"
/**
<p>
<p>
@return
the total number of Hs (implicit and explicit) that this Atom is bound to
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Atom::hasProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Atom::Match 	( 	const ATOM_SPTR  	what 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>
Reimplemented in RDKit::QueryAtom.
*/
public";

%javamethodmodifiers RDKit::Atom::Match 	( 	Atom const *  	what 	 )  	const "
/**
<p>
<p>
@return
whether or not we match the argument
<p>
@notes
<li>for Atom objects, 'match' means that atomic numbers are the same.
Reimplemented in RDKit::QueryAtom.
*/
public";

%javamethodmodifiers RDKit::Atom::setIdx 	( 	unsigned int  	index 	 )  	"
/**
<p>
sets our index within the ROMol
<p>
<p>
@notes
<li>this makes no sense if we do not have an owning molecule
<li>the index should be < this->getOwningMol()->getNumAtoms()
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Atom::setProp 	( 	const std::string  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Atom::setProp 	( 	const char *  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
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

%javamethodmodifiers RDKit::Atom::updatePropertyCache 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates any of our lazy properties
<p>
<p>
@notes
<li>requires an owning molecule
<li>the current lazy properties are implicit and explicit valence

*/
public";

