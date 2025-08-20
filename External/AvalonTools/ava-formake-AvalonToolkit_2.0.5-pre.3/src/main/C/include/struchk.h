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

/************************************************************************/
/*                                                                      */
/*  File:     struchk.h                                                 */
/*                                                                      */
/*  Purpose:  Header file for public interface of struchk module.       */
/*                                                                      */
/************************************************************************/

#define BAD_MOLECULE		1
#define ALIAS_CONVERSION_FAILED	2
#define TRANSFORMED		4
#define FRAGMENTS_FOUND		8
#define EITHER_WARNING		16
#define STEREO_ERROR		32
#define DUBIOUS_STEREO_REMOVED	64
#define ATOM_CLASH		128
#define ATOM_CHECK_FAILED	256
#define SIZE_CHECK_FAILED	512
#define RECHARGED		1024
#define STEREO_FORCED_BAD	2048
#define STEREO_TRANSFORMED	4096
#define TEMPLATE_TRANSFORMED	8192
#define TAUTOMER_TRANSFORMED	16384

#define BAD_SET	(BAD_MOLECULE 			| \
                 ALIAS_CONVERSION_FAILED 	| \
		 STEREO_ERROR			| \
		 STEREO_FORCED_BAD		| \
		 ATOM_CLASH			| \
		 ATOM_CHECK_FAILED		| \
		 SIZE_CHECK_FAILED)

#define TRANSFORMED_SET (TRANSFORMED		| \
			 FRAGMENTS_FOUND	| \
			 EITHER_WARNING		| \
			 DUBIOUS_STEREO_REMOVED	| \
			 STEREO_TRANSFORMED	| \
			 TEMPLATE_TRANSFORMED	| \
			 TAUTOMER_TRANSFORMED	| \
			 RECHARGED)

/**
 * Wrapper around RunStruchk() to be called from Oracle call-outs.
 * Note: InitCheckMol() may have been called before to initialize
 * options.
 */
extern int CheckMol(struct reaccs_molecule_t *mol);

/**
 * This function converts parameters from single, linefeed-separated
 * string representation to array of strings. This is an entry point
 * that is supposed to be called from the Oracle call-out glue code.
 */
extern int InitCheckMol(char *opt);


/**************************************************************
 * CloseOpenFiles closes the open global files
 **************************************************************/
 void CloseOpenFiles(void);


// AG */
/* -ft */
extern char *tmp_dir_name;          /* name for a tmp directory */
