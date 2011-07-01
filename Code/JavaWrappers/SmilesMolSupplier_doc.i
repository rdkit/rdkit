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

%typemap(javaimports) RDKit::SmilesMolSupplier "
/** 
lazy file parser for Smiles tables */"

%javamethodmodifiers RDKit::SmilesMolSupplier::SmilesMolSupplier 	( 	const std::string &  	fileName, 		const std::string &  	delimiter = ' \t', 		int  	smilesColumn = 0, 		int  	nameColumn = 1, 		bool  	titleLine = true, 		bool  	sanitize = true	  	) 			"
/**
<p>
<p>
@param
fileName 	- the name of smiles table file
delimiter 	- delimiting characters between records on a each line NOTE that this is not a string, the tokenizer looks for the individual characters in delimiter, not the full string itself. So the default delimiter: ' \t', means ' ' or '\t'.
smilesColumn 	- column number for the SMILES string (defaults to the first column)
nameColumn 	- column number for the molecule name (defaults to the second column) If set to -1 we assume that no name is available for the molecule and the name is defaulted to the smiles string
titleLine 	- if true, the first line is assumed to list the names of properties in order seperated by 'delimiter'. It is also assume that the 'SMILES' column and the 'name' column are not specified here if false - no title line is assumed and the properties are recorded as the 'columnX' where 'X' is the cloumn number
sanitize 	- if true sanitize the molecule before returning it

*/
public";

%javamethodmodifiers RDKit::SmilesMolSupplier::getItemText 	( 	unsigned int  	idx 	 )  	"
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

