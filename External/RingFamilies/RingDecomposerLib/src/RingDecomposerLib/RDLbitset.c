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

#include <string.h>
#include <stdlib.h>

#include "RDLbitset.h"

/* add two compressed cycles inplace */
void RDL_bitset_xor_inplace(unsigned char* dst, unsigned const char* src, unsigned size)
{
  unsigned i;

  for (i = 0; i < size; ++i) {
    dst[i] ^= src[i];
  }
}

void RDL_bitset_or_inplace(unsigned char* dst, unsigned const char* src, unsigned size)
{
  unsigned i;

  for (i = 0; i < size; ++i) {
    dst[i] |= src[i];
  }
}

/* check if compressed cycle is empty (uses an empty template for efficiency */
int RDL_bitset_empty(const unsigned char* cycle, const unsigned char* empty_cycle, unsigned size)
{
  return !memcmp(cycle, empty_cycle, size * sizeof(*cycle));
}

/* define block size (should be 8 bit) */
static const unsigned BlockSize = 8*sizeof(unsigned char);

/* set position of compressed bitset */
void RDL_bitset_set(unsigned char* bitset, unsigned pos)
{
  bitset[pos / BlockSize] |= (unsigned char)(1 << (pos % BlockSize));
}

unsigned RDL_bitset_init(unsigned char** compressed, unsigned size)
{
  unsigned new_size;

  new_size = size / BlockSize;
  if (size % BlockSize) {
    ++new_size;
  }

  *compressed = malloc(new_size * sizeof(**compressed));
  memset(*compressed, 0, new_size * sizeof(**compressed));

  return new_size;
}

/* make a compressed bitset */
unsigned RDL_bitset_compressed(unsigned char** compressed, const char* edges, unsigned size)
{
  unsigned new_size, i;

  new_size = RDL_bitset_init(compressed, size);

  for (i = 0; i < size; ++i) {
    if (edges[i]) {
      RDL_bitset_set(*compressed, i);
    }
  }

  return new_size;
}

/* flip position of compressed bitset */
static void RDL_bitset_flip(unsigned char* bitset, unsigned pos)
{
  bitset[pos / BlockSize] ^= (unsigned char)(1 << (pos % BlockSize));
}

/* swap columns of compressed bitset */
void RDL_swap_columns(unsigned char** rows, unsigned nof_rows, unsigned col1, unsigned col2)
{
  unsigned char val1, val2;
  unsigned i;

  for (i = 0; i < nof_rows; ++i) {
    val1 = RDL_bitset_test(rows[i], col1);
    val2 = RDL_bitset_test(rows[i], col2);

    if (!val1 != !val2) {
      RDL_bitset_flip(rows[i], col1);
      RDL_bitset_flip(rows[i], col2);
    }
  }
}

/* test position of compressed bitset */
unsigned char RDL_bitset_test(unsigned char* bitset, unsigned pos)
{
  return (bitset[pos / BlockSize] >> (pos % BlockSize)) & 1;
}
