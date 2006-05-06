
////////////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or              //
// modify it under the terms of the GNU General Public License                //
// version 2 as published by the Free Software Foundation.                    //
//                                                                            //
// This program is distributed in the hope that it will be useful, but        //
// WITHOUT ANY WARRANTY; without even the implied warranty of                 //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          //
// General Public License for more details.                                   //
//                                                                            //
// Written by François Fleuret                                                //
// Contact <francois.fleuret@epfl.ch> for comments & bug reports              //
// Copyright (C) 2004 EPFL                                                    //
////////////////////////////////////////////////////////////////////////////////

// $Id: fastentropy.h,v 1.1 2005/03/03 15:52:35 fleuret Exp $

#ifndef FASTENTROPY_H
#define FASTENTROPY_H

#ifdef _MSC_VER
#include <windows.h>
typedef UINT32 uint32_t;
typedef UINT16 uint16_t;
#else
   #include <stdint.h>
#endif


// Lookup tables to speed up the training

extern int fe_nb_bits[65536];
extern double fe_logn[65536], fe_nlogn[65536];

void fe_init_tables();

inline void fe_set_bit(int k, uint32_t *a, int v) {
  if(v) a[k/32] = a[k/32] |  (1 << (k%32));
  else  a[k/32] = a[k/32] & ~(1 << (k%32));
}

inline int fe_get_bit(int k, uint32_t *a) {
  return a[k/32] & (1 << (k%32));
}

inline int fe_count_and(int n, uint32_t *a, uint32_t *b) {
  uint16_t *aa = (uint16_t *) a, *bb = (uint16_t *) b;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[*aa++ & *bb++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & *aa & *bb];
  return t;
}

inline int fe_count_and_not(int n, uint32_t *a, uint32_t *b) {
  uint16_t *aa = (uint16_t *) a,  *bb = (uint16_t *) b;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[*aa++ & ~*bb++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & *aa & ~*bb];
  return t;
}

inline int fe_count_and_not_not(int n, uint32_t *a, uint32_t *b) {
  uint16_t *aa = (uint16_t *) a, *bb = (uint16_t *) b;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[65535 & ~*aa++ & ~*bb++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & ~*aa & ~*bb];
  return t;
}

// This selection maximises the conditional mutual information between
// the features and the class to predict, given any feature already
// picked

void fe_selection_cmim(int nb_samples,
                       int nb_total_features, uint32_t **x, uint32_t *y,
                       int nb_selected, int *selected);

// This selection maximises the mutual information between the
// features and the class to predict, without taking care of the
// redundancy

void fe_selection_mim(int nb_samples,
                      int nb_total_features, uint32_t **x, uint32_t *y,
                      int nb_selected, int *selected);

#endif
