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

%typemap(javaimports) RDKit::RingInfo "
/** 
A class to store information about a molecule's rings. */"

%javamethodmodifiers RDKit::RingInfo::addRing 	( 	const INT_VECT &  	atomIndices, 		const INT_VECT &  	bondIndices	  	) 			"
/**
<p>
adds a ring to our data
<p>
<p>
@param
atomIndices 	the integer indices of the atoms involved in the ring
bondIndices 	the integer indices of the bonds involved in the ring, this must be the same size as atomIndices.
<p>
@return
the number of rings
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::atomRings 	( 		 )  	const "
/**
<p>
<p>
@return
our atom-rings vectors
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::bondRings 	( 		 )  	const "
/**
<p>
<p>
@return
our bond-rings vectors
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::isAtomInRingOfSize 	( 	unsigned int  	idx, 		unsigned int  	size	  	) 			const"
/**
<p>
<p>
@return
whether or not the atom with index idx is in a size - ring.
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::isBondInRingOfSize 	( 	unsigned int  	idx, 		unsigned int  	size	  	) 			const"
/**
<p>
<p>
@return
whether or not the bond with index idx is in a size - ring.
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::minAtomRingSize 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
<p>
@return
the size of the smallest ring atom idx is involved in
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::minBondRingSize 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
<p>
@return
the size of the smallest ring bond idx is involved in
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::numAtomRings 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
<p>
@return
the number of rings atom idx is involved in
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::numBondRings 	( 	unsigned int  	idx 	 )  	const"
/**
<p>
<p>
@return
the number of rings bond idx is involved in
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

%javamethodmodifiers RDKit::RingInfo::numRings 	( 		 )  	const"
/**
<p>
<p>
@return
the total number of rings
<p>
@notes
<li>the object must be initialized before calling this

*/
public";

