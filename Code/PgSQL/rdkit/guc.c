// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
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
//       products derived from this software without specific prior written permission.
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
#include "postgres.h"
#include "fmgr.h"
#include "utils/guc.h"

#include "rdkit.h"

static double rdkit_tanimoto_smlar_limit = 0.5;
static double rdkit_dice_smlar_limit = 0.5;
static double rdkit_agent_FP_bit_ratio = 0.2;
static double rdkit_threshold_unmapped_reactant_atoms = 0.2;
static bool rdkit_do_chiral_sss = false;
static bool rdkit_ignore_reaction_agents = false;
static bool rdkit_move_unmmapped_reactants_to_agents = true;
static bool rdkit_init_reaction = true;
static bool rdkit_guc_inited = false;

#define SSS_FP_SIZE 2048
#define LAYERED_FP_SIZE 1024
#define RDKIT_FP_SIZE 1024
#define MORGAN_FP_SIZE 512
#define FEATMORGAN_FP_SIZE 512
#define HASHED_TORSION_FP_SIZE 1024
#define HASHED_PAIR_FP_SIZE 2048
#define AVALON_FP_SIZE 512
#define REACTION_SSS_FP_SIZE 4096
#define REACTION_DIFFERENCE_FP_SIZE 2048
#define REACTION_SSS_FP_TYPE 5
#define REACTION_DIFFERENCE_FP_TYPE 1
#define REACTION_DFP_WEIGHT_AGENTS 1
#define REACTION_DFP_WEIGHT_NONAGENTS 10

static int rdkit_sss_fp_size = SSS_FP_SIZE;
static int rdkit_morgan_fp_size = MORGAN_FP_SIZE;
static int rdkit_featmorgan_fp_size = FEATMORGAN_FP_SIZE;
static int rdkit_layered_fp_size = LAYERED_FP_SIZE;
static int rdkit_rdkit_fp_size = RDKIT_FP_SIZE;
static int rdkit_hashed_torsion_fp_size = HASHED_TORSION_FP_SIZE;
static int rdkit_hashed_atompair_fp_size = HASHED_PAIR_FP_SIZE;
static int rdkit_avalon_fp_size = AVALON_FP_SIZE;
static int rdkit_reaction_sss_fp_size = REACTION_SSS_FP_SIZE;
static int rdkit_reaction_difference_fp_size = REACTION_DIFFERENCE_FP_SIZE;
static int rdkit_reaction_sss_fp_type = REACTION_SSS_FP_TYPE;
static int rdkit_reaction_difference_fp_type = REACTION_DIFFERENCE_FP_TYPE;
static int rdkit_difference_FP_weight_agents = REACTION_DFP_WEIGHT_AGENTS;
static int rdkit_difference_FP_weight_nonagents = REACTION_DFP_WEIGHT_NONAGENTS;

#if PG_VERSION_NUM < 80400
#error The earliest supported postgresql version is 8.4
#endif

static void
initRDKitGUC()
{
  if (rdkit_guc_inited)
    return;

  DefineCustomRealVariable(
                           "rdkit.tanimoto_threshold",
                           "Lower threshold of Tanimoto similarity",
                           "Molecules with similarity lower than threshold are not similar by % operation",
                           &rdkit_tanimoto_smlar_limit,
                           0.5,
                           0.0,
                           1.0,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomRealVariable(
                           "rdkit.dice_threshold",
                           "Lower threshold of Dice similarity",
                           "Molecules with similarity lower than threshold are not similar by # operation",
                           &rdkit_dice_smlar_limit,
                           0.5,
                           0.0,
                           1.0,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomBoolVariable(
                           "rdkit.do_chiral_sss",
                           "Should stereochemistry be taken into account in substructure matching",
                           "If false (the default), no stereochemistry information is used in substructure matches.",
                           &rdkit_do_chiral_sss,
                           false,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );

  DefineCustomIntVariable(
                           "rdkit.sss_fp_size",
                           "Size (in bits) of the fingerprint used for substructure screening",
                           "Size (in bits) of the fingerprint used for substructure screening",
                           &rdkit_sss_fp_size,
                           SSS_FP_SIZE,
                           64,
                           4096,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.morgan_fp_size",
                           "Size (in bits) of morgan fingerprints",
                           "Size (in bits) of morgan fingerprints",
                           &rdkit_morgan_fp_size,
                           MORGAN_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.featmorgan_fp_size",
                           "Size (in bits) of featmorgan fingerprints",
                           "Size (in bits) of featmorgan fingerprints",
                           &rdkit_featmorgan_fp_size,
                           FEATMORGAN_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.layered_fp_size",
                           "Size (in bits) of layered fingerprints",
                           "Size (in bits) of layered fingerprints",
                           &rdkit_layered_fp_size,
                           LAYERED_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.rdkit_fp_size",
                           "Size (in bits) of RDKit fingerprints",
                           "Size (in bits) of RDKit fingerprints",
                           &rdkit_rdkit_fp_size,
                           RDKIT_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.hashed_torsion_fp_size",
                           "Size (in bits) of topological torsion bit vector fingerprints",
                           "Size (in bits) of topological torsion bit vector fingerprints",
                           &rdkit_hashed_torsion_fp_size,
                           HASHED_TORSION_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.hashed_atompair_fp_size",
                           "Size (in bits) of atom pair bit vector fingerprints",
                           "Size (in bits) of atom pair torsion bit vector fingerprints",
                           &rdkit_hashed_atompair_fp_size,
                           HASHED_PAIR_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.reaction_sss_fp_size",
                           "Size (in bits) of structural chemical reaction fingerprint",
                           "Size (in bits) of structural chemical reaction fingerprint",
                           &rdkit_reaction_sss_fp_size,
                           REACTION_SSS_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.reaction_difference_fp_size",
                           "Size (in bits) of difference chemical reaction fingerprint",
                           "Size (in bits) of difference chemical reaction fingerprint",
                           &rdkit_reaction_difference_fp_size,
                           REACTION_DIFFERENCE_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.reaction_sss_fp_type",
                           "Type of structural chemical reaction fingerprint (1=Atompairs, 2=Topologcal torsion, 3=Morgan2, 4=RDKit, 5=Pattern)",
                           "Type of structural chemical reaction fingerprint (1=Atompairs, 2=Topologcal torsion, 3=Morgan2, 4=RDKit, 5=Pattern)",
                           &rdkit_reaction_sss_fp_type,
                           REACTION_SSS_FP_TYPE,
                           1,
                           5,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.reaction_difference_fp_type",
                           "Type of difference chemical reaction fingerprint (1=Atompairs, 2=Topologcal torsion, 3=Morgan2)",
                           "Type of difference chemical reaction fingerprint (1=Atompairs, 2=Topologcal torsion, 3=Morgan2)",
                           &rdkit_reaction_difference_fp_type,
                           REACTION_DIFFERENCE_FP_TYPE,
                           1,
                           3,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomBoolVariable(
                             "rdkit.ignore_reaction_agents",
                             "Should agents of a chemical reaction be taken into account ",
                             "If false (the default), agents (if provided) are taken into account for calculations/requests/similarity measures.",
                             &rdkit_ignore_reaction_agents,
                             false,
                             PGC_USERSET,
                             0,
  			   NULL,
  #if PG_VERSION_NUM >= 90000
                             NULL,
  #endif
                             NULL
                             );
  DefineCustomRealVariable(
                            "rdkit.agent_FP_bit_ratio",
                            "Allow to weight the impact of agents contained in a chemical reaction fingerprint",
                            "A scaling factor below 1.0 down weighted the agents of a chemical reaction",
                            &rdkit_agent_FP_bit_ratio,
                            0.2,
                            0.0,
                            3.0,
                            PGC_USERSET,
                            0,
 			   NULL,
 #if PG_VERSION_NUM >= 90000
                            NULL,
 #endif
                            NULL
                            );
  DefineCustomBoolVariable(
                             "rdkit.move_unmmapped_reactants_to_agents",
                             "Should unmapped reactant  agents of a chemical reaction be taken into account ",
                             "If true (the default), agents (if provided) are taken into account for calculations/requests/similarity measures.",
                             &rdkit_move_unmmapped_reactants_to_agents,
                             true,
                             PGC_USERSET,
                             0,
  			   NULL,
  #if PG_VERSION_NUM >= 90000
                             NULL,
  #endif
                             NULL
                             );
  DefineCustomRealVariable(
                             "rdkit.threshold_unmapped_reactant_atoms",
                             "Set the ratio of allowed unmapped reactant atoms",
                             "If the ratio of unmapped atoms in a reactant is larger than this the reactant is removed and added to the agents.",
                             &rdkit_threshold_unmapped_reactant_atoms,
                             0.2,
                             0.0,
                             1.0,
                             PGC_USERSET,
                             0,
  			   NULL,
  #if PG_VERSION_NUM >= 90000
                             NULL,
  #endif
                             NULL
                             );
  DefineCustomBoolVariable(
                             "rdkit.init_reaction",
                             "Checks if the reaction is ready for use",
                             "By default it is true, however reactions can be used uninitialized too but be aware of potential upcoming errors.",
                             &rdkit_init_reaction,
                             true,
                             PGC_USERSET,
                             0,
  			   NULL,
  #if PG_VERSION_NUM >= 90000
                             NULL,
  #endif
                             NULL
                             );
  DefineCustomIntVariable(
                           "rdkit.difference_FP_weight_agents",
                           "In reaction difference fingerprints weight factor for agents comapred to reactants and products",
                           "In reaction difference fingerprints weight factor for agents comapred to reactants and products",
                           &rdkit_difference_FP_weight_agents,
                           REACTION_DFP_WEIGHT_AGENTS,
                           -10,
                           10,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  DefineCustomIntVariable(
                           "rdkit.difference_FP_weight_nonagents",
                           "In reaction difference fingerprints weight factor for reactants and products comapred to agents",
                           "In reaction difference fingerprints weight factor for reactants and products comapred to agents",
                           &rdkit_difference_FP_weight_nonagents,
                           REACTION_DFP_WEIGHT_NONAGENTS,
                           1,
                           20,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );

  DefineCustomIntVariable(
                           "rdkit.avalon_fp_size",
                           "Size (in bits) of avalon fingerprints",
                           "Size (in bits) of avalon fingerprints",
                           &rdkit_avalon_fp_size,
                           AVALON_FP_SIZE,
                           64,
                           9192,
                           PGC_USERSET,
                           0,
			   NULL,
#if PG_VERSION_NUM >= 90000
                           NULL,
#endif
                           NULL
                           );
  rdkit_guc_inited = true;
}

double
getTanimotoLimit(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();

  return rdkit_tanimoto_smlar_limit;
}

double
getDiceLimit(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();

  return rdkit_dice_smlar_limit;
}

bool
getDoChiralSSS(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();

  return rdkit_do_chiral_sss;
}

int
getSubstructFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_sss_fp_size;
}
int
getMorganFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_morgan_fp_size;
}
int
getFeatMorganFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_featmorgan_fp_size;
}
int
getLayeredFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_layered_fp_size;
}
int
getRDKitFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_rdkit_fp_size;
}
int
getHashedTorsionFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_hashed_torsion_fp_size;
}
int
getHashedAtomPairFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_hashed_atompair_fp_size;
}
int
getAvalonFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_avalon_fp_size;
}
int
getReactionSubstructFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_reaction_sss_fp_size;
}
int
getReactionDifferenceFpSize(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_reaction_difference_fp_size;
}
int
getReactionSubstructFpType(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_reaction_sss_fp_type;
}
int
getReactionDifferenceFpType(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_reaction_difference_fp_type;
}
bool
getIgnoreReactionAgents(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_ignore_reaction_agents;
}
double
getReactionStructuralFPAgentBitRatio(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_agent_FP_bit_ratio;
}
bool
getMoveUnmappedReactantsToAgents(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_move_unmmapped_reactants_to_agents;
}
double
getThresholdUnmappedReactantAtoms(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_threshold_unmapped_reactant_atoms;
}
bool
getInitReaction(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_init_reaction;
}
int
getReactionDifferenceFPWeightAgents(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_difference_FP_weight_agents;
}
int
getReactionDifferenceFPWeightNonagents(void) {
  if (!rdkit_guc_inited)
    initRDKitGUC();
  return rdkit_difference_FP_weight_nonagents;
}

void _PG_init(void);
void
_PG_init(void) {
  initRDKitGUC();
}
