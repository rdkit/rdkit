// $Id$
//
//  Copyright (c) 2010-2015, Novartis Institutes for BioMedical Research Inc.
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
#ifndef _RDKIT_H_
#define _RDKIT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <postgres.h>

typedef bytea Mol;

#define DatumGetMolP(x) ((Mol *)PG_DETOAST_DATUM(x))
#define DatumGetMolPCopy(x) ((Mol *)PG_DETOAST_DATUM_COPY(x))
#define MolPGetDatum(x) (PointerGetDatum(x))

#define PG_GETARG_MOL_P(x) DatumGetMolP(PG_GETARG_DATUM(x))
#define PG_GETARG_MOL_P_COPY(x) DatumGetMolPCopy(PG_GETARG_DATUM(x))
#define PG_RETURN_MOL_P(x) PG_RETURN_DATUM(MolPGetDatum(x))

typedef bytea Bfp;

typedef struct {
  char vl_len_[4];
  uint16 weight;
  uint8 fp[FLEXIBLE_ARRAY_MEMBER];
} BfpSignature;

#define DatumGetBfpP(x) ((Bfp *)PG_DETOAST_DATUM(x))
#define DatumGetBfpPCopy(x) ((Bfp *)PG_DETOAST_DATUM_COPY(x))
#define BfpPGetDatum(x) (PointerGetDatum(x))

#define PG_GETARG_BFP_P(x) DatumGetBfpP(PG_GETARG_DATUM(x))
#define PG_GETARG_BFP_P_COPY(x) DatumGetBfpPCopy(PG_GETARG_DATUM(x))
#define PG_RETURN_BFP_P(x) PG_RETURN_DATUM(BfpPGetDatum(x))

#define BFP_SIGLEN(x) (VARSIZE(x) - VARHDRSZ)
  
typedef bytea Sfp;

#define DatumGetSfpP(x) ((Sfp *)PG_DETOAST_DATUM(x))
#define DatumGetSfpPCopy(x) ((Sfp *)PG_DETOAST_DATUM_COPY(x))
#define SfpPGetDatum(x) (PointerGetDatum(x))

#define PG_GETARG_SFP_P(x) DatumGetSfpP(PG_GETARG_DATUM(x))
#define PG_GETARG_SFP_P_COPY(x) DatumGetSfpPCopy(PG_GETARG_DATUM(x))
#define PG_RETURN_SFP_P(x) PG_RETURN_DATUM(SfpPGetDatum(x))

typedef bytea Reaction;

#define DatumGetReactionP(x) ((Reaction *)PG_DETOAST_DATUM(x))
#define DatumGetReactionPCopy(x) ((Reaction *)PG_DETOAST_DATUM_COPY(x))
#define ReactionPGetDatum(x) (PointerGetDatum(x))

#define PG_GETARG_REACTION_P(x) DatumGetReactionP(PG_GETARG_DATUM(x))
#define PG_GETARG_REACTION_P_COPY(x) DatumGetReactionPCopy(PG_GETARG_DATUM(x))
#define PG_RETURN_REACTION_P(x) PG_RETURN_DATUM(ReactionPGetDatum(x))

/*
 * From/to C/C++
 */

/* RDKit::ROMol */
typedef void *CROMol;
void freeCROMol(CROMol data);

CROMol constructROMol(Mol *data);
Mol *deconstructROMol(CROMol data);

CROMol parseMolBlob(char *data, int len);
char *makeMolBlob(CROMol data, int *len);
CROMol parseMolText(char *data, bool asSmarts, bool warnOnFail, bool asQuery);
CROMol parseMolCTAB(char *data, bool keepConformer, bool warnOnFail,
                    bool asQuery);
char *makeMolText(CROMol data, int *len, bool asSmarts);
char *makeCtabText(CROMol data, int *len, bool createDepictionIfMissing);
bool isValidSmiles(char *data);
bool isValidSmarts(char *data);
bool isValidCTAB(char *data);
bool isValidMolBlob(char *data, int len);

int molcmp(CROMol i, CROMol a);

int MolSubstruct(CROMol i, CROMol a);
int MolSubstructCount(CROMol i, CROMol a, bool uniquify);

bytea *makeMolSignature(CROMol data);

double MolAMW(CROMol i);
double MolLogP(CROMol i);
int MolHBA(CROMol i);
int MolHBD(CROMol i);
int MolNumAtoms(CROMol i);
int MolNumHeavyAtoms(CROMol i);
int MolNumRotatableBonds(CROMol i);
int MolNumHeteroatoms(CROMol i);
int MolNumRings(CROMol i);
int MolNumAromaticRings(CROMol i);
int MolNumAliphaticRings(CROMol i);
int MolNumSaturatedRings(CROMol i);
int MolNumAromaticHeterocycles(CROMol i);
int MolNumAliphaticHeterocycles(CROMol i);
int MolNumSaturatedHeterocycles(CROMol i);
int MolNumAromaticCarbocycles(CROMol i);
int MolNumAliphaticCarbocycles(CROMol i);
int MolNumSaturatedCarbocycles(CROMol i);
int MolNumHeterocycles(CROMol i);

double MolFractionCSP3(CROMol i);
double MolTPSA(CROMol i);
double MolChi0v(CROMol i);
double MolChi1v(CROMol i);
double MolChi2v(CROMol i);
double MolChi3v(CROMol i);
double MolChi4v(CROMol i);
double MolChi0n(CROMol i);
double MolChi1n(CROMol i);
double MolChi2n(CROMol i);
double MolChi3n(CROMol i);
double MolChi4n(CROMol i);
double MolKappa1(CROMol i);
double MolKappa2(CROMol i);
double MolKappa3(CROMol i);

int MolNumSpiroAtoms(CROMol i);
int MolNumBridgeheadAtoms(CROMol i);

char *makeMolFormulaText(CROMol data, int *len, bool separateIsotopes,
                         bool abbreviateHIsotopes);

const char *MolInchi(CROMol i, const char *opts);
const char *MolInchiKey(CROMol i, const char *opts);
CROMol MolMurckoScaffold(CROMol i);

CROMol MolAdjustQueryProperties(CROMol m, const char *params);

/* ExplicitBitVect */
typedef void *CBfp;
void freeCBfp(CBfp data);

CBfp constructCBfp(Bfp *data);
Bfp *deconstructCBfp(CBfp data);
BfpSignature *makeBfpSignature(CBfp data);

int CBfpSize(CBfp a);

double calcBitmapTanimotoSml(CBfp a, CBfp b);
double calcBitmapDiceSml(CBfp a, CBfp b);
double calcBitmapTverskySml(CBfp a, CBfp b, float ca, float cb);

/* SparseIntVect<boost::int32_t> */
typedef void *CSfp;
void freeCSfp(CSfp data);

CSfp constructCSfp(Sfp *data);
Sfp *deconstructCSfp(CSfp data);
bytea *makeSfpSignature(CSfp data, int numBits);
bytea *makeLowSparseFingerPrint(CSfp data, int numInts);

double calcSparseTanimotoSml(CSfp a, CSfp b);
double calcSparseDiceSml(CSfp a, CSfp b);
double calcSparseStringDiceSml(const char *a, unsigned int sza, const char *b,
                               unsigned int szb);
bool calcSparseStringAllValsGT(const char *a, unsigned int sza, int tgt);
bool calcSparseStringAllValsLT(const char *a, unsigned int sza, int tgt);
CSfp addSFP(CSfp a, CSfp b);
CSfp subtractSFP(CSfp a, CSfp b);

void countOverlapValues(bytea *sign, CSfp data, int numBits,
                        int *sum, int *overlapSum, int *overlapN);
void countLowOverlapValues(bytea *sign, CSfp data, int numInts,
                           int *querySum, int *keySum, int *overlapUp,
                           int *overlapDown);
/*
 * Various mol -> fp transformation
 */

CBfp makeLayeredBFP(CROMol data);
CBfp makeRDKitBFP(CROMol data);
CBfp makeMorganBFP(CROMol data, int radius);
CSfp makeMorganSFP(CROMol data, int radius);
CBfp makeFeatMorganBFP(CROMol data, int radius);
CSfp makeFeatMorganSFP(CROMol data, int radius);
CSfp makeAtomPairSFP(CROMol data);
CSfp makeTopologicalTorsionSFP(CROMol data);
CBfp makeAtomPairBFP(CROMol data);
CBfp makeTopologicalTorsionBFP(CROMol data);
CBfp makeMACCSBFP(CROMol data);
CBfp makeAvalonBFP(CROMol data, bool isQuery, unsigned int bitFlags);

/*
 * Indexes
 */

#define NUMBITS (2048)
#define NUMRANGE (120)

#define INTRANGEMAX (0xff)
typedef struct IntRange {
  uint8 low;
  uint8 high;
} IntRange;

#define RDKitTanimotoStrategy (1)
#define RDKitDiceStrategy (2)
#define RDKitOrderByTanimotoStrategy (3)
#define RDKitOrderByDiceStrategy (4)
#define RDKitContains (3)
#define RDKitContained (4)
#define RDKitEquals (6)
#define RDKitSmaller (7)
#define RDKitGreater (8)

bool calcConsistency(bool isLeaf, uint16 strategy, double nCommonUp,
                     double nCommonDown, double nKey, double nQuery);

/* Chemical Reactions
 * RDKit::ChemicalReaction */
typedef void *CChemicalReaction;

void freeChemReaction(CChemicalReaction data);

CChemicalReaction constructChemReact(Reaction *data);
Reaction *deconstructChemReact(CChemicalReaction data);

CChemicalReaction parseChemReactBlob(char *data, int len);
CChemicalReaction parseChemReactText(char *data, bool asSmarts,
                                     bool warnOnFail);
CChemicalReaction parseChemReactCTAB(char *data, bool warnOnFail);
char *makeChemReactBlob(CChemicalReaction data, int *len);
char *makeChemReactText(CChemicalReaction data, int *len, bool asSmarts);
char *makeCTABChemReact(CChemicalReaction data, int *len);

int ChemReactNumReactants(CChemicalReaction rxn);
int ChemReactNumProducts(CChemicalReaction rxn);
int ChemReactNumAgents(CChemicalReaction rxn);

/* Reaction substructure search */
bytea *makeReactionSign(CChemicalReaction data);
int ReactionSubstruct(CChemicalReaction rxn, CChemicalReaction rxn2);
int reactioncmp(CChemicalReaction rxn, CChemicalReaction rxn2);
CBfp makeReactionBFP(CChemicalReaction data, int size, int fpType);

/* Reaction difference fingerprint */
CSfp makeReactionDifferenceSFP(CChemicalReaction data, int size, int fpType);

char *computeMolHash(CROMol data, int *len);

char *findMCSsmiles(char *smiles, char *params);
void *addMol2list(void *lst, Mol *mol);
char *findMCS(void *lst, char *params);

#ifdef __cplusplus
}
#endif
#endif
