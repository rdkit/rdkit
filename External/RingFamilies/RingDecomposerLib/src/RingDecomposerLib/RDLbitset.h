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

#ifndef RDL_BITSET_H
#define RDL_BITSET_H

void RDL_bitset_xor_inplace(unsigned char* dst, const unsigned char* src, unsigned size);
void RDL_bitset_or_inplace(unsigned char* dst, const unsigned char* src, unsigned size);
int RDL_bitset_empty(unsigned const char* cycle, unsigned const char* empty_cycle, unsigned size);

void RDL_bitset_set(unsigned char* cycle, unsigned pos);

void RDL_swap_columns(unsigned char** rows, unsigned nof_rows, unsigned col1, unsigned col2);
unsigned RDL_bitset_init(unsigned char** compressed, unsigned size);
unsigned RDL_bitset_compressed(unsigned char** compressed, const char* edges, unsigned size);

unsigned char RDL_bitset_test(unsigned char* bitset, unsigned pos);


#endif
