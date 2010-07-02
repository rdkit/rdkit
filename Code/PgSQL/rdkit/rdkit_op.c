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
#include "rdkit.h"
#include "fmgr.h"

/***************** Mol operations ***********************/
static int 
molcmp(Mol *a, Mol *b) {
	int res;

	res = memcmp(VARDATA(a), VARDATA(b), Min(VARSIZE(a), VARSIZE(b)) - VARHDRSZ);
	if ( res )
		return res;

	if (VARSIZE(a) == VARSIZE(b))
		return 0;
	return (VARSIZE(a) > VARSIZE(b)) ? 1 : -1; 
}


#define MOLCMPFUNC( type, action, ret )                 	\
PG_FUNCTION_INFO_V1(mol_##type);                        	\
Datum		mol_##type(PG_FUNCTION_ARGS);           		\
Datum                                                   	\
mol_##type(PG_FUNCTION_ARGS)                       			\
{                                                       	\
    Mol    *a, *b;											\
	int		res;											\
															\
	fcinfo->flinfo->fn_extra = SearchMolCache(				\
								fcinfo->flinfo->fn_extra,	\
								fcinfo->flinfo->fn_mcxt,	\
								PG_GETARG_DATUM(0), 		\
								&a, NULL, NULL, NULL);		\
	fcinfo->flinfo->fn_extra = SearchMolCache(				\
								fcinfo->flinfo->fn_extra,	\
								fcinfo->flinfo->fn_mcxt,	\
								PG_GETARG_DATUM(1), 		\
								&b, NULL, NULL, NULL);		\
    res = molcmp(a, b);                         			\
    PG_RETURN_##ret( res action 0 );                    	\
}   \
/* keep compiler quiet - no extra ; */                  	\
extern int no_such_variable

MOLCMPFUNC(lt, <, BOOL);
MOLCMPFUNC(le, <=, BOOL);
MOLCMPFUNC(eq, ==, BOOL);
MOLCMPFUNC(ge, >=, BOOL);
MOLCMPFUNC(gt, >, BOOL);
MOLCMPFUNC(ne, !=, BOOL);
MOLCMPFUNC(cmp, +, INT32);




PG_FUNCTION_INFO_V1(mol_tanimoto_sml);
Datum		mol_tanimoto_sml(PG_FUNCTION_ARGS);
Datum
mol_tanimoto_sml(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, NULL, &bfp, NULL);

	res = calcTanimotoSml(afp, bfp); 

	PG_RETURN_FLOAT8(res);		
}

PG_FUNCTION_INFO_V1(mol_tanimoto_sml_op);
Datum		mol_tanimoto_sml_op(PG_FUNCTION_ARGS);
Datum
mol_tanimoto_sml_op(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, NULL, &bfp, NULL);

	res = calcTanimotoSml(afp, bfp); 
	PG_RETURN_BOOL( res >= getTanimotoLimit() );
}

PG_FUNCTION_INFO_V1(mol_dice_sml);
Datum		mol_dice_sml(PG_FUNCTION_ARGS);
Datum
mol_dice_sml(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, NULL, &bfp, NULL);

	res = calcDiceSml(afp, bfp); 

	PG_RETURN_FLOAT8(res);		
}

PG_FUNCTION_INFO_V1(mol_dice_sml_op);
Datum		mol_dice_sml_op(PG_FUNCTION_ARGS);
Datum
mol_dice_sml_op(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, NULL, &bfp, NULL);

	res = calcDiceSml(afp, bfp); 
	PG_RETURN_BOOL( res >= getDiceLimit() );
}

/***************** Fingerprint operations ***********************/

static int 
fpcmp(FingerPrint *a, FingerPrint *b) {
	int res;

	res = memcmp(VARDATA(a), VARDATA(b), Min(VARSIZE(a), VARSIZE(b)) - VARHDRSZ);
	if ( res )
		return res;

	if (VARSIZE(a) == VARSIZE(b))
		return 0;
	return (VARSIZE(a) > VARSIZE(b)) ? 1 : -1; 
}


#define FPCMPFUNC( type, action, ret )         	        	\
PG_FUNCTION_INFO_V1(fp_##type);                	        	\
Datum		fp_##type(PG_FUNCTION_ARGS); 	   		    	\
Datum                                                   	\
fp_##type(PG_FUNCTION_ARGS)                       			\
{                                                       	\
    FingerPrint    *a, *b;									\
	int		res;											\
															\
	fcinfo->flinfo->fn_extra = SearchFPCache(				\
								fcinfo->flinfo->fn_extra,	\
								fcinfo->flinfo->fn_mcxt,	\
								PG_GETARG_DATUM(0), 		\
								&a, NULL, NULL);			\
	fcinfo->flinfo->fn_extra = SearchFPCache(				\
								fcinfo->flinfo->fn_extra,	\
								fcinfo->flinfo->fn_mcxt,	\
								PG_GETARG_DATUM(1), 		\
								&b, NULL, NULL);			\
    res = fpcmp(a, b); 		                        		\
    PG_RETURN_##ret( res action 0 );                    	\
}   \
/* keep compiler quiet - no extra ; */                  	\
extern int no_such_variable

FPCMPFUNC(lt, <, BOOL);
FPCMPFUNC(le, <=, BOOL);
FPCMPFUNC(eq, ==, BOOL);
FPCMPFUNC(ge, >=, BOOL);
FPCMPFUNC(gt, >, BOOL);
FPCMPFUNC(ne, !=, BOOL);
FPCMPFUNC(cmp, +, INT32);

PG_FUNCTION_INFO_V1(fp_tanimoto_sml);
Datum		fp_tanimoto_sml(PG_FUNCTION_ARGS);
Datum
fp_tanimoto_sml(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &bfp, NULL);

	res = calcTanimotoSml(afp, bfp); 

	PG_RETURN_FLOAT8(res);		
}

PG_FUNCTION_INFO_V1(fp_tanimoto_sml_op);
Datum		fp_tanimoto_sml_op(PG_FUNCTION_ARGS);
Datum
fp_tanimoto_sml_op(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &bfp, NULL);

	res = calcTanimotoSml(afp, bfp); 
	PG_RETURN_BOOL( res >= getTanimotoLimit() );
}

PG_FUNCTION_INFO_V1(fp_dice_sml);
Datum		fp_dice_sml(PG_FUNCTION_ARGS);
Datum
fp_dice_sml(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &bfp, NULL);

	res = calcDiceSml(afp, bfp); 

	PG_RETURN_FLOAT8(res);		
}

PG_FUNCTION_INFO_V1(fp_dice_sml_op);
Datum		fp_dice_sml_op(PG_FUNCTION_ARGS);
Datum
fp_dice_sml_op(PG_FUNCTION_ARGS) {
	MolFingerPrint 	afp,
					bfp;
	double 			res;

	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &afp, NULL);
	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &bfp, NULL);

	res = calcDiceSml(afp, bfp); 
	PG_RETURN_BOOL( res >= getDiceLimit() );
}

PG_FUNCTION_INFO_V1(fp_size);
Datum		fp_size(PG_FUNCTION_ARGS);
Datum
fp_size(PG_FUNCTION_ARGS) {
	MolFingerPrint	fp;

	fcinfo->flinfo->fn_extra = SearchFPCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &fp, NULL);

	PG_RETURN_INT32(MolFingerPrintSize(fp));
}

PG_FUNCTION_INFO_V1(mol_substruct);
Datum		mol_substruct(PG_FUNCTION_ARGS);
Datum
mol_substruct(PG_FUNCTION_ARGS) {
	CROMol 	i,
			a;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &i, NULL, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &a, NULL, NULL);

	PG_RETURN_BOOL(MolSubstruct(i, a));		
}

PG_FUNCTION_INFO_V1(mol_rsubstruct);
Datum		mol_rsubstruct(PG_FUNCTION_ARGS);
Datum
mol_rsubstruct(PG_FUNCTION_ARGS) {
	CROMol 	i,
			a;

	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(0), 
								NULL, &a, NULL, NULL);
	fcinfo->flinfo->fn_extra = SearchMolCache(
								fcinfo->flinfo->fn_extra,
								fcinfo->flinfo->fn_mcxt,
								PG_GETARG_DATUM(1), 
								NULL, &i, NULL, NULL);

	PG_RETURN_BOOL(MolSubstruct(i, a));		
}

