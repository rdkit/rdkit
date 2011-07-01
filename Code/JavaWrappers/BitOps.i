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
#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>
%}

%include <DataStructs/BitVect.h>
%include <DataStructs/BitOps.h>

/* Only the following are tested */
/* We expect that SWIG can figure out the overloading so we don't need to add a suffix like EBV to the name */
%template(OnBitProjSimilarity) OnBitProjSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(FoldFingerprint) FoldFingerprint<ExplicitBitVect>;

%template(AllProbeBitsMatch) AllProbeBitsMatch<ExplicitBitVect>;
%template(AllProbeBitsMatch) AllProbeBitsMatch<ExplicitBitVect>;
%template(TanimotoSimilarity) TanimotoSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(CosineSimilarity) CosineSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(KulczynskiSimilarity) KulczynskiSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(DiceSimilarity) DiceSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(TverskySimilarity) TverskySimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(SokalSimilarity) SokalSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(McConnaugheySimilarity) McConnaugheySimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(AsymmetricSimilarity) AsymmetricSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(BraunBlanquetSimilarity) BraunBlanquetSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(RusselSimilarity) RusselSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(OnBitSimilarity) OnBitSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(NumBitsInCommon) NumBitsInCommon<ExplicitBitVect,ExplicitBitVect>;
%template(AllBitSimilarity) AllBitSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(OnBitsInCommon) OnBitsInCommon<ExplicitBitVect,ExplicitBitVect>;
%template(OffBitsInCommon) OffBitsInCommon<ExplicitBitVect,ExplicitBitVect>;
%template(OnBitProjSimilarity) OnBitProjSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(OffBitProjSimilarity) OffBitProjSimilarity<ExplicitBitVect,ExplicitBitVect>;
%template(BitVectToText) BitVectToText<ExplicitBitVect>;

