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


%{
#include <boost/cstdint.hpp>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
%}

#if SWIGCSHARP
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(UInt_Pair, std::pair<unsigned int,int> );
SWIG_STD_VECTOR_SPECIALIZE_MINIMUM(UInt_Pair_Vect, std::vector< std::pair<unsigned int,int> >);
#endif

//%include <DataStructs/SparseIntVect.h>


%newobject RDKit::MorganFingerprints::getFingerprint;
%rename(MorganFingerprintMol) RDKit::MorganFingerprints::getFingerprint;
%newobject RDKit::MorganFingerprints::getFingerprintAsBitVect;
%rename(getMorganFingerprintAsBitVect) RDKit::MorganFingerprints::getFingerprintAsBitVect;
%include <GraphMol/Fingerprints/MorganFingerprints.h>

%include <DataStructs/BitOps.h>
%template(TanimotoSimilarityEBV) TanimotoSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilarity) DiceSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilarity) RDKit::DiceSimilarity<boost::uint32_t>;
%template(DiceSimilarity) RDKit::DiceSimilarity<boost::int32_t>;
%template(DiceSimilarity) RDKit::DiceSimilarity<boost::int64_t>;
%template(TanimotoSimilaritySIVu32) RDKit::TanimotoSimilarity<boost::uint32_t>;
%template(TanimotoSimilaritySIVi32) RDKit::TanimotoSimilarity<boost::int32_t>;
%template(TanimotoSimilaritySIVi64) RDKit::TanimotoSimilarity<boost::int64_t>;
