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

#include "RDLutility.h"

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

RDL_outputFunction RDL_outputFunc = RDL_writeNothing;

unsigned **RDL_alloc2DUIntArray(unsigned n, unsigned m)
{
  unsigned i;
  unsigned **arr = malloc(n * sizeof(*arr));
  for(i=0; i<n; ++i)
  {
    arr[i] = malloc(m*sizeof(**arr));
  }
  return arr;
}

char **RDL_alloc2DCharArray(unsigned n, unsigned m)
{
  unsigned i;
  char **arr = malloc(n * sizeof(*arr));
  for(i=0; i<n; ++i)
  {
    arr[i] = malloc(m*sizeof(**arr));
  }
  return arr;
}

void RDL_delete2DUIntArray(unsigned **arr, unsigned n)
{
  unsigned i;
  for (i = 0; i < n; ++i) {
    free(arr[i]);
  }
  free(arr);
}

void RDL_delete2DCharArray(char **arr, unsigned n)
{
  unsigned i;
  for (i = 0; i < n; ++i) {
    free(arr[i]);
  }
  free(arr);
}

void RDL_print2DUIntArray(unsigned **arr, unsigned n, unsigned m)
{
  unsigned i,j;
  for(i=0; i<n; ++i)
  {
    for(j=0; j<m; ++j)
    {
      printf("%d ",arr[i][j]);
    }
    printf("\n");
  }
}

void RDL_print2DCharArray(char **arr, unsigned n, unsigned m)
{
  unsigned i,j;
  for(i=0; i<n; ++i)
  {
    for(j=0; j<m; ++j)
    {
      printf("%d ",arr[i][j]);
    }
    printf("\n");
  }
}

void RDL_writeNothing(RDL_ERROR_LEVEL level, const char* fmt, ...)
{

}

void RDL_writeToStderr(RDL_ERROR_LEVEL level, const char* fmt, ...)
{
  unsigned len;
  char* new_fmt;
  char* specifier;

  len = strlen(fmt);
  len += 1 + 14;

  /* we don't need any output initialization, so skip */
  if (level == RDL_INITIALIZE_OUTPUT) {
    return;
  }

  switch(level) {
    case RDL_DEBUG:
      specifier = "RDL_DEBUG";
      break;
    case RDL_WARNING:
      specifier = "RDL_WARNING";
      break;
    case RDL_ERROR:
      specifier = "RDL_ERROR";
      break;
    default:
      specifier = "RDL_????";
  }

  new_fmt = malloc(len * sizeof(*new_fmt));
  sprintf(new_fmt, "%12s: %s", specifier, fmt);

  va_list va;
  va_start(va, fmt);
  vfprintf(stderr, new_fmt, va);
  va_end(va);

  free(new_fmt);
}
