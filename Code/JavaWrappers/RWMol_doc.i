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

%typemap(javaimports) RDKit::RWMol "
/** 
RWMol is a molecule class that is intended to be edited.
<p>
See documentation for ROMol for general remarks */"

%javamethodmodifiers RDKit::RWMol::addAtom 	( 	ATOM_SPTR  	atom, 		bool  	updateLabel = true	  	) 			"
/**
<p>
adds an Atom to our collection
<p>
<p>
@param
atom 	pointer to the Atom to add
updateLabel 	(optional) if this is true, the new Atom will be our activeAtom
<p>
@return
the new number of atoms
<p>
@notes
Reimplemented from RDKit::ROMol.
*/
public";

%javamethodmodifiers RDKit::RWMol::addAtom 	( 	Atom *  	atom, 		bool  	updateLabel = true, 		bool  	takeOwnership = false	  	) 			"
/**
<p>
adds an Atom to our collection
<p>
<p>
@param
atom 	pointer to the Atom to add
updateLabel 	(optional) if this is true, the new Atom will be our activeAtom
takeOwnership 	(optional) if this is true, we take ownership of atom instead of copying it.
<p>
@return
the new number of atoms
Reimplemented from RDKit::ROMol.
*/
public";

%javamethodmodifiers RDKit::RWMol::addAtom 	( 	bool  	updateLabel = true 	 )  	"
/**
<p>
adds an empty Atom to our collection
<p>
<p>
@param
updateLabel 	(optional) if this is true, the new Atom will be our activeAtom
<p>
@return
the new number of atoms

*/
public";

%javamethodmodifiers RDKit::RWMol::addBond 	( 	BOND_SPTR  	bsp 	 )  	"
/**
<p>
adds a Bond to our collection
<p>
<p>
@param
bond 	pointer to the Bond to add
<p>
@return
the new number of bonds
<p>
@notes
Reimplemented from RDKit::ROMol.
*/
public";

%javamethodmodifiers RDKit::RWMol::addBond 	( 	Bond *  	bond, 		bool  	takeOwnership = false	  	) 			"
/**
<p>
adds a Bond to our collection
<p>
<p>
@param
bond 	pointer to the Bond to add
takeOwnership 	(optional) if this is true, we take ownership of bond instead of copying it.
<p>
@return
the new number of bonds
Reimplemented from RDKit::ROMol.
*/
public";

%javamethodmodifiers RDKit::RWMol::addBond 	( 	Atom *  	beginAtom, 		Atom *  	endAtom, 		Bond::BondType  	order = Bond::UNSPECIFIED	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::RWMol::addBond 	( 	Atom::ATOM_SPTR  	atom1, 		Atom::ATOM_SPTR  	atom2, 		Bond::BondType  	order = Bond::UNSPECIFIED	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::RWMol::addBond 	( 	unsigned int  	beginAtomIdx, 		unsigned int  	endAtomIdx, 		Bond::BondType  	order = Bond::UNSPECIFIED	  	) 			"
/**
<p>
adds a Bond between the indicated Atoms
<p>
<p>
@return
the number of Bonds

*/
public";

%javamethodmodifiers RDKit::RWMol::createPartialBond 	( 	unsigned int  	beginAtomIdx, 		Bond::BondType  	order = Bond::UNSPECIFIED	  	) 			"
/**
<p>
starts a Bond and sets its beginAtomIdx
<p>
<p>
@return
a pointer to the new bond
The caller should set a bookmark to the returned Bond in order to be able to later complete it:
<p>
        Bond *pBond = mol->createPartialBond(1);
	mol->setBondBookmark(pBond,666);
	... do some other stuff ...
	mol->finishPartialBond(2,666,Bond::SINGLE);
	mol->clearBondBookmark(666,pBond);
      
<p>
or, if we want to set the BondType initially:
<p>
        Bond *pBond = mol->createPartialBond(1,Bond::DOUBLE);
	mol->setBondBookmark(pBond,666);
	... do some other stuff ...
	mol->finishPartialBond(2,666);
	mol->clearBondBookmark(666,pBond);
      
<p>
the call to finishPartialBond() will take priority if you set the BondType in both calls.
*/
public";

%javamethodmodifiers RDKit::RWMol::finishPartialBond 	( 	unsigned int  	endAtomIdx, 		int  	bondBookmark, 		Bond::BondType  	order = Bond::UNSPECIFIED	  	) 			"
/**
<p>
finishes a partially constructed bond
<p>
<p>
@return
the final number of Bonds
See the documentation for createPartialBond() for more details
*/
public";

%javamethodmodifiers RDKit::RWMol::getActiveAtom 	( 		 )  	"
/**
<p>
<p>
@return
a pointer to the 'active' Atom
If we have an activeAtom, it will be returned, otherwise the results of getLastAtom() will be returned.
*/
public";

%javamethodmodifiers RDKit::RWMol::removeAtom 	( 	Atom *  	atom 	 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::RWMol::replaceAtom 	( 	unsigned int  	idx, 		Atom *  	atom, 		bool  	updateLabel = false	  	) 			"
/**
<p>
replaces a particular Atom
<p>
<p>
@param
idx 	the index of the Atom to replace
atom 	the new atom, which will be copied.
updateLabel 	(optional) if this is true, the new Atom will be our activeAtom

*/
public";

%javamethodmodifiers RDKit::RWMol::setActiveAtom 	( 	unsigned int  	idx 	 )  	"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>

*/
public";

