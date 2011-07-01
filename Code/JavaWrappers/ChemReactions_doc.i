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

%typemap(javaimports) RDKit::ChemicalReaction "
/** 
This is a class for storing and applying general chemical reactions.
<p>
<p>
@example
<pre><code>
<p>

</code></pre>
basic usage will be something like:
<p>
     ChemicalReaction rxn;
     rxn.addReactantTemplate(r1);
     rxn.addReactantTemplate(r2);
     rxn.addProductTemplate(p1);
     rxn.initReactantMatchers();
     
     MOL_SPTR_VECT prods;
     for(MOL_SPTR_VECT::const_iterator r1It=reactantSet1.begin();
         r1It!=reactantSet1.end();++r1It;){
       for(MOL_SPTR_VECT::const_iterator r2It=reactantSet2.begin();
	   r2It!=reactantSet2.end();++r2It;){
	 MOL_SPTR_VECT rVect(2);
         rVect[0] = *r1It;
         rVect[1] = *r2It;
             
         std::vector<MOL_SPTR_VECT> lprods;
         lprods = rxn.runReactants(rVect);
         for(std::vector<MOL_SPTR_VECT>::const_iterator lpIt=lprods.begin();
            lpIt!=lprods.end();++lpIt){
            // we know this is a single-product reaction:
            prods.push_back((*lpIt)[0]);
         }
       }     
     }
     
 */"

%javamethodmodifiers RDKit::ChemicalReaction::addProductTemplate 	( 	ROMOL_SPTR  	mol 	 )  	"
/**
<p>
Adds a new product template.
<p>
<p>
@return
the number of products

*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::addReactantTemplate 	( 	ROMOL_SPTR  	mol 	 )  	"
/**
<p>
Adds a new reactant template.
<p>
<p>
@return
the number of reactants

*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::getImplicitPropertiesFlag 	( 		 )  	const "
/**
<p>
<p>
@return
whether or not the reaction uses implicit properties on the product atoms
This toggles whether or not unspecified atomic properties in the products are considered to be implicit and should be copied from the actual reactants. This is necessary due to a semantic difference between the 'reaction SMARTS' approach and the MDL RXN approach: In 'reaction SMARTS', this reaction: [C:1]-[Br:2].[O-:3]>>[C:1]-[O:3].[Br-:2] applied to [CH4+]Br should yield [CH4+]O Something similar drawn in an rxn file, and applied to [CH4+]Br should yield [CH3]O. In rxn there is no charge on the product C because nothing is specified in the rxn file; in 'SMARTS' the charge from the actual reactants is not *removed* because no charge is specified in the reaction.
*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::initReactantMatchers 	( 		 )  	"
/**
<p>
initializes our internal reactant-matching datastructures.
<p>
This must be called after adding reactants and before calling runReactants.
*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::runReactants 	( 	const MOL_SPTR_VECT  	reactants 	 )  	const"
/**
<p>
Runs the reaction on a set of reactants.
<p>
<p>
@param
reactants,: 	the reactants to be used. The length of this must be equal to this->getNumReactantTemplates()
<p>
@return
a vector of vectors of products. Each subvector will be this->getNumProductTemplates() long.
We return a vector of vectors of products because each individual template may map multiple times onto its reactant. This leads to multiple possible result sets.
*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::setImplicitPropertiesFlag 	( 	bool  	val 	 )  	"
/**
<p>
sets the implicit properties flag. See the documentation for getImplicitProertiesFlag() for a discussion of what this means.
*/
public";

%javamethodmodifiers RDKit::ChemicalReaction::validate 	( 	unsigned int &  	numWarnings, 		unsigned int &  	numErrors, 		bool  	silent = false	  	) 			const"
/**
<p>
validates the reactants and products to make sure the reaction seems 'reasonable'
<p>
<p>
@return
true if the reaction validates without errors (warnings do not stop validation)
<p>
@param
numWarnings,: 	used to return the number of validation warnings
numErrors,: 	used to return the number of validation errors
silent,: 	If this bool is true, no messages will be logged during the validation. By default, validation problems are reported to the warning and error logs depending on their severity.

*/
public";

