/* 
*
*  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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

%typemap(javaimports) RDKit::FilterCatalogParams "
/**
<p>
Parameter Class for initializing a FilterCatalog with the specified set of FilterCatalogs. Current
catalogs include:

<pre><code>
  FilterCatalogParams.FilterCatalogs.ALL
  FilterCatalogParams.FilterCatalogs.BRENK
  FilterCatalogParams.FilterCatalogs.PAINS_A
  FilterCatalogParams.FilterCatalogs.PAINS_B
  FilterCatalogParams.FilterCatalogs.PAINS_C
  FilterCatalogParams.FilterCatalogs.PAINS
  FilterCatalogParams.FilterCatalogs.ZINC
</pre></code>

Details:

<h4>Brenk</h4>
<ul>
<li>Reference: <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/18064617\">Brenk</a> Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases</li>
<li>Scope: unwanted functionality due to potential tox reasons or unfavourable pharmacokinetic properties</li>
</ul>

<h4>NIH</h4>
<ul>
<li>Reference: <a href=\"http://pubs.rsc.org/en/Content/ArticleLanding/2015/OB/c4ob02287d#!divAbstract\">NIH</a> A Unified Lead-oriented Synthesis of over Fifty Molecular Scaffolds. </li>
<li>Reference: <a href=\"http://pubs.acs.org/doi/abs/10.1021/jm901070c\">NIH</a> Quantitative Analyses of Aggregation, Autofluorescence, and Reactivity Artifacts in a Screen for Inhibitors of a Thiol Protease</li>
<li>Scope: uannotate compounds with problematic functional groups </li>
</ul>

<h4>PAINS</h4>
<ul>
<li>Reference: <a href=\"http://pubs.rsc.org/en/Content/ArticleLanding/2015/OB/c4ob02287d#!divAbstract\">NIH</a> New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays </li>
<li>Scope: PAINS filters </li>
</ul>

<h4>ZINC</h4>
<ul>
<li>Reference: <a href=\"http://blaster.docking.org/filtering/\">ZINC</a></li>
<li>Scope: drug-likeness and unwanted functional group filters </li>
</ul>

**/"
