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

%typemap(javaimports) RDKit::TDTMolSupplier "
/** 
lazy file parser for TDT files */"

%javamethodmodifiers RDKit::TDTMolSupplier::TDTMolSupplier 	( 	const std::string &  	fileName, 		const std::string &  	nameRecord = "", 		int  	confId2D = -1, 		int  	confId3D = 0, 		bool  	sanitize = true	  	) 			"
/**
<p>
<p>
@param
fileName 	- the name of the TDT file
nameRecord 	- property name for the molecule name. If empty (the default), the name defaults to be empty
confId2D 	- if >=0 and 2D coordinates are provided, the 2D structure (depiction) in the input will be read into the corresponding conformer id.
confId3D 	- if >=0 and 3D coordinates are provided, the 3D structure (depiction) in the input will be read into the corresponding conformer id.
sanitize 	- if true sanitize the molecule before returning it

*/
public";

%javamethodmodifiers RDKit::TDTMolSupplier::getItemText 	( 	unsigned int  	idx 	 )  	"
/**
<p>
<p>
@return
the text block for a particular item
<p>
@param
idx 	- which item to return

*/
public";

