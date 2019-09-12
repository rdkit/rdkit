/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

#ifndef RDL_TESTING_H
#define RDL_TESTING_H

#include "RingDecomposerLib.h"

#ifdef __cplusplus
extern "C" {
#endif

RDL_API
int RDL_no_stop_fun(const void* nothing);

RDL_API
int RDL_timeout_stop_fun(const void* p_timeout);

/*
 * validate the URFs of given ringsystem,
 * StopFunc specifies function to interrupt the process
 * return 0 if successful, 1 on difference, -1 on timeout
 */
RDL_API
int RDL_validateRingFamilies(const char** cycle_array,
                             unsigned* urf_numbers,
                             unsigned nof_rings,
                             unsigned nof_bonds,
                             int (*stop) (const void*));

/*
 * check internal consistency of results
 */
RDL_API
int RDL_checkConsistency(const RDL_data* data);

/*
 * compare two cycles (terminated by a 2)!
 */
RDL_API
int RDL_cmp_cycles(const void* vc1, const void* vc2);

/*
 * validate that a basis is indeed a basis
 */
RDL_API
int RDL_validateBasis(const char** relevant_cycle_array,
                      unsigned nof_rc,
                      const char** basis_cycle_arry,
                      unsigned nof_basis,
                      unsigned nof_bonds,
                      unsigned* bcc,
                      int (*stop) (const void*));


#ifdef __cplusplus
}
#endif

#endif
