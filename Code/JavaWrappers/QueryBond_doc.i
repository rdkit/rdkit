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

%typemap(javaimports) RDKit::QueryBond "
/** 
Class for storing Bond queries.
<p>
QueryBond objects are derived from Bond objects, so they can be added to molecules and the like, but they have much fancier querying capabilities. */"

%javamethodmodifiers RDKit::QueryBond::expandQuery 	( 	QUERYBOND_QUERY *  	what, 		Queries::CompositeQueryType  	how = Queries::COMPOSITE_AND, 		bool  	maintainOrder = true	  	) 			"
/**
<p>
expands our current query
<p>
<p>
@param
what 	the Queries::Query to be added
how 	the operator to be used in the expansion
maintainOrder 	(optional) flags whether the relative order of the queries needs to be maintained, if this is false, the order is reversed
<p>
@notes
<li>what should probably be constructed using one of the functions defined in QueryOps.h
<li>the maintainOrder option can be useful because the combination operators short circuit when possible.
Reimplemented from RDKit::Bond.
*/
public";

%javamethodmodifiers RDKit::QueryBond::Match 	( 	Bond const *  	what 	 )  	const "
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>
Reimplemented from RDKit::Bond.
<p>

*/
public";

