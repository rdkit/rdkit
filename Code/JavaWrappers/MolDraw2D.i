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
%include "boost_tuple.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_pair.i"
%include "boost_tuple.i"
%{
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <map>
#include <vector>
#include <utility>
%}

%template(Int_String_Map) std::map< int, std::string >;
%template(DrawColour) boost::tuple<float,float,float>;

#ifdef SWIGJAVA
%template(ColourPalette) std::map< int, RDKit::DrawColour >;
#endif

%template(Int_Double_Map) std::map< int, double >;
%template(Float_Pair) std::pair<float,float>;
%template(Float_Pair_Vect) std::vector< std::pair<float,float> >;

%ignore RDKit::MolDraw2DSVG::MolDraw2DSVG(int,int,std::ostream &);

%include <GraphMol/MolDraw2D/MolDraw2D.h>
%include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
%include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
