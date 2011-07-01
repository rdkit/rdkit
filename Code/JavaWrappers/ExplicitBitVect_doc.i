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

%typemap(javaimports) RDKit::ExplicitBitVect "
/** 
a class for bit vectors that are densely occupied
<p>
RDKit::ExplicitBitVect objects store all of their bits using a boost::dynamic_bitset
<p>
These are very fast, but can require large amounts of memory for large, sparsely occupied vectors.  */"

%javamethodmodifiers RDKit::ExplicitBitVect (unsigned int size)"
/**
 	initialize with a particular size;
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect (const RDKit::ExplicitBitVect &other)"
/**

*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect (const std::string &)"
/**
 	construct from a string pickle
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect (const char *, const unsigned int)"
/**
 	construct from a text pickle
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::setBit (const unsigned int which)"
/**
 	sets a particular bit and returns its original value
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::unsetBit (const unsigned int which)"
/**
 	unsets a particular bit and returns its original value
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::getBit (const unsigned int which) const"
/**
 	returns the value of a particular bit
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::getNumBits () const"
/**
 	returns the number of bits (the length of the BitVect)
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::getNumOnBits () const"
/**
 	returns the number of on bits
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::getNumOffBits () const"
/**
 	returns the number of off bits
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::getOnBits (IntVect &v) const"
/**
 	replaces the contents of v with indices of our on bits
*/
public";

%javamethodmodifiers RDKit::ExplicitBitVect::clearBits ()"
/**
 	clears (sets to off) all of our bits
<p>

*/
public";

