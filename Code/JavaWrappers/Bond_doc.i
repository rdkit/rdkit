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

%typemap(javaimports) RDKit::Bond "
/** 
class for representing a bond
<p>
<p>
@notes
<li>many of the methods of Atom require that the Atom be associated with a molecule (an ROMol).
<li>each Bond maintains a Dict of properties:
<li>o Each property is keyed by name and can store an arbitrary type.
<li>o Properties can be marked as calculated, in which case they will be cleared when the clearComputedProps() method is called.
<li>o Because they have no impact upon chemistry, all property operations are const, this allows extra flexibility for clients who need to store extra data on Bond objects.
 */"

%javamethodmodifiers RDKit::Bond::clearProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Bond::clearProp 	( 	const char *  	key 	 )  	const "
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

%javamethodmodifiers RDKit::Bond::copy 	( 		 )  	const "
/**
<p>
<p>
@return
a copy
<p>
@notes
Reimplemented in RDKit::QueryBond.
*/
public";

%javamethodmodifiers RDKit::Bond::getBeginAtom 	( 		 )  	const"
/**
<p>
<p>
@return
a pointer to our begin Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getBeginAtomIdx 	( 		 )  	const "
/**
<p>
<p>
@return
the index of our begin Atom
<p>
@notes
<li>this makes no sense if we do not have an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getEndAtom 	( 		 )  	const"
/**
<p>
<p>
@return
a pointer to our end Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getEndAtomIdx 	( 		 )  	const "
/**
<p>
<p>
@return
the index of our end Atom
<p>
@notes
<li>this makes no sense if we do not have an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getIdx 	( 		 )  	const "
/**
<p>
<p>
@return
our index within the ROMol
<p>
@notes
<li>this makes no sense if we do not have an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getOtherAtom 	( 	Atom const *  	what 	 )  	const"
/**
<p>
<p>
@return
a pointer to the other Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::getOtherAtomIdx 	( 	unsigned int  	thisIdx 	 )  	const"
/**
<p>
given the index of one Atom, returns the index of the other
<p>
<p>
@notes
<li>this makes no sense if we do not have an owning molecule
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Bond::getProp 	( 	const std::string  	key, 		T &  	res	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Bond::getProp 	( 	const char *  	key, 		T &  	res	  	) 			const "
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

%javamethodmodifiers RDKit::Bond::getStereoAtoms 	( 		 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Bond::getValenceContrib 	( 	const Atom *  	at 	 )  	const"
/**
<p>
<p>
@return
our contribution to the explicit valence of an Atom
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::hasProp 	( 	const std::string  	key 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Bond::Match 	( 	const Bond::BOND_SPTR  	what 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>
Reimplemented in RDKit::QueryBond.
*/
public";

%javamethodmodifiers RDKit::Bond::Match 	( 	Bond const *  	what 	 )  	const "
/**
<p>
<p>
@return
whether or not we match the argument
<p>
@notes
<li>for Bond objects, 'match' means that either one of the Bonds has bondType Bond::UNSPECIFIED or both Bonds have the same bondType.
Reimplemented in RDKit::QueryBond.
*/
public";

%javamethodmodifiers RDKit::Bond::setBeginAtom 	( 	Atom::ATOM_SPTR  	at 	 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Bond::setBeginAtom 	( 	Atom *  	at 	 )  	"
/**
<p>
sets our begin Atom
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::setBeginAtomIdx 	( 	unsigned int  	what 	 )  	"
/**
<p>
sets the index of our begin Atom
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::setEndAtom 	( 	Atom::ATOM_SPTR  	at 	 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::Bond::setEndAtom 	( 	Atom *  	at 	 )  	"
/**
<p>
sets our end Atom
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::setEndAtomIdx 	( 	unsigned int  	what 	 )  	"
/**
<p>
sets the index of our end Atom
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::setIdx 	( 	unsigned int  	index 	 )  	"
/**
<p>
sets our index within the ROMol
<p>
<p>
@notes
<li>this makes no sense if we do not have an owning molecule
<li>the index should be < this->getOwningMol()->getNumBonds()
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Bond::setProp 	( 	const std::string  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
template<typename T >
*/
public";

%javamethodmodifiers RDKit::Bond::setProp 	( 	const char *  	key, 		T  	val, 		bool  	computed = false	  	) 			const "
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

%javamethodmodifiers RDKit::Bond::updatePropertyCache 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates any of our lazy properties
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

%javamethodmodifiers RDKit::Bond::updatePropertyCache 	( 	bool  	strict = true 	 )  	"
/**
<p>
calculates any of our lazy properties
<p>
<p>
@notes
<li>requires an owning molecule

*/
public";

