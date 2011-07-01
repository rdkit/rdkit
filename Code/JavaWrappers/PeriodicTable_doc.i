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

%typemap(javaimports) RDKit::PeriodicTable "
/** 
singleton class for retrieving information about atoms
<p>
Use the singleton like this:
<p>
    const PeriodicTable *tbl = PeriodicTable::getTable();
    tbl->getAtomicWeight(6); // get atomic weight for Carbon
    tbl->getAtomicWeight('C'); // get atomic weight for Carbon
    
 */"

%javamethodmodifiers RDKit::PeriodicTable::getAtomicWeight 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getAtomicWeight 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getDefaultValence 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getDefaultValence 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getNouterElecs 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getNouterElecs 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRb0 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRb0 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRcovalent 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRcovalent 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRvdw 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getRvdw 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getValenceList 	( 	char *  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getValenceList 	( 	const std::string &  	elementSymbol 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::getValenceList 	( 	UINT  	atomicNumber 	 )  	const "
/**
<p>
<p>
@return
a vector of all stable valences. For atoms where we really don't have any idea what a reasonable maximum valence is (like transition metals), the vector ends with -1
1
*/
public";

%javamethodmodifiers RDKit::PeriodicTable::moreElectroNegative 	( 	UINT  	anum1, 		UINT  	anum2	  	) 			const "
/**
<p>
convenience function to determine which atom is more electronegative
<p>
check if atom with atomic number anum1 is more electronegative than the one with anum2 this is rather lame but here is how we do it
<p>
    * the atom with the higher number of outer shell electrons is considered more electronegative
    * if the # of outer shell elecs are the same the atom with the lower atomic weight is more electronegative
<p>


*/
public";

