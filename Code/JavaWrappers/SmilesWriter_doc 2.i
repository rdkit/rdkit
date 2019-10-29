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

%typemap(javaimports) RDKit::SmilesWriter "
/** 
The SmilesWriter is for writing molecules and properties to delimited text files. */"

%javamethodmodifiers RDKit::SmilesWriter::SmilesWriter 	( 	std::string  	fileName, 		std::string  	delimiter = ' ', 		std::string  	nameHeader = 'Name', 		bool  	includeHeader = true, 		bool  	isomericSmiles = false, 		bool  	kekuleSmiles = false	  	) 			"
/**
<p>
<p>
@param
fileName 	: filename to write to ('-' to write to stdout)
delimiter 	: delimiter to use in the text file
nameHeader 	: used to label the name column in the output. If this is provided as the empty string, no names will be written.
includeHeader 	: toggles inclusion of a header line in the output
isomericSmiles 	: toggles generation of isomeric SMILES
kekuleSmiles 	: toggles the generation of kekule SMILES

*/
public";

%javamethodmodifiers RDKit::SmilesWriter::SmilesWriter 	( 	std::ostream *  	outStream, 		std::string  	delimiter = ' ', 		std::string  	nameHeader = 'Name', 		bool  	includeHeader = true, 		bool  	takeOwnership = false, 		bool  	isomericSmiles = false, 		bool  	kekuleSmiles = false	  	) 			"
/**
<p>
This is an overloaded member function, provided for convenience. It differs from the above function only in what argument(s) it accepts.
<p>

*/
public";

