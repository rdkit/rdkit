//
//  Copyright (c) 2010-2021 Novartis Institutes for BioMedical Research Inc.
//    and other RDKit contributors
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <RDGeneral/export.h>
#ifndef RDKIT_GUC_H_
#define RDKIT_GUC_H_

#ifdef __cplusplus
extern "C" {
#endif

double getTanimotoLimit(void);
double getDiceLimit(void);
bool getDoChiralSSS(void);
bool getDoEnhancedStereoSSS(void);
int getSubstructFpSize(void);
int getMorganFpSize(void);
int getFeatMorganFpSize(void);
int getLayeredFpSize(void);
int getRDKitFpSize(void);
int getHashedTorsionFpSize(void);
int getHashedAtomPairFpSize(void);
int getAvalonFpSize(void);
int getReactionSubstructFpSize(void);
int getReactionSubstructFpType(void);
int getReactionDifferenceFpSize(void);
int getReactionDifferenceFpType(void);
bool getIgnoreReactionAgents(void);
double getReactionStructuralFPAgentBitRatio(void);
bool getMoveUnmappedReactantsToAgents(void);
double getThresholdUnmappedReactantAtoms(void);
bool getInitReaction(void);
int getReactionDifferenceFPWeightNonagents(void);
int getReactionDifferenceFPWeightAgents(void);

#ifdef __cplusplus
}
#endif
#endif
