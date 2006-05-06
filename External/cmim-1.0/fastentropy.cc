
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

// $Id: fastentropy.cc,v 1.2 2005/03/03 20:16:15 fleuret Exp $


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#include "misc.h"
#include "fastentropy.h"

const double FLOAT_TOL=1e-8;

int fe_nb_bits[65536];
double fe_logn[65536], fe_nlogn[65536];

#ifdef DEBUG
int fe_was_initialized = 0;
#endif

inline void fe_and(int n, uint32_t *a, uint32_t *b, uint32_t *ab) {
  for(int k = 0; k < (n+31)/32; k++) *ab++ = *a++ & *b++;
}

inline void fe_and_not(int n, uint32_t *a, uint32_t *b, uint32_t *ab) {
  for(int k = 0; k < (n+31)/32; k++) *ab++ = *a++ & ~*b++;
}

inline void fe_and_not_not(int n, uint32_t *a, uint32_t *b, uint32_t *ab) {
  for(int k = 0; k < (n+31)/32; k++) *ab++ = ~*a++ & ~*b++;
}

inline int fe_count(int n, uint32_t *a) {
  uint16_t *aa = (uint16_t *) a;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[*aa++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & *aa];
  return t;
}

inline double fe_entropy(int n, uint32_t *a) {
  int n1 = fe_count(n, a);
  return fe_logn[n] - (fe_nlogn[n1] + fe_nlogn[n-n1])/double(n);
}

inline double fe_entropy_couple(int n, uint32_t *a, uint32_t *b) {
  int n11 = fe_count_and(n, a, b);
  int n10 = fe_count_and_not(n, a, b);
  int n01 = fe_count_and_not(n, b, a);
  int n00 = fe_count_and_not_not(n, a, b);

  return fe_logn[n] - (fe_nlogn[n00] + fe_nlogn[n01] + fe_nlogn[n10] + fe_nlogn[n11])/double(n);
}

void fe_init_tables() {
#ifdef DEBUG
  fe_was_initialized = 1;
#endif
  for(int i = 0; i < 65536; i++) {
    if(i == 0) { fe_logn[i] = 0.0; fe_nlogn[i] = 0.0; }
    else { fe_logn[i] = log(double(i)); fe_nlogn[i] = double(i) * log(double(i)); }
    int n = 0;
    for(int j = 0; j < 16; j++) if(i & (1 << j)) n++;
    fe_nb_bits[i] = n;
  }
}

void fe_selection_cmim(int nb_samples,
                       int nb_tests, uint32_t **x, uint32_t *y,
                       int nb_selected, int *selected) {

#ifdef DEBUG
  if(!fe_was_initialized) {
    cerr << "fe_init_tables() was not called!\n";
    abort();
  }
#endif

  if(nb_samples > 65535) {
    cerr << "Too many pictures, the nlogn table is too small.\n";
    exit(1);
  }

  double *s = new double[nb_tests];
  double *ch = new double[nb_tests];
  int *m = new int[nb_tests];

  double h = fe_entropy(nb_samples, y);

  for(int i = 0; i < nb_tests; i++) {
    ch[i] = fe_entropy_couple(nb_samples, y, x[i]) - fe_entropy(nb_samples, x[i]);
    s[i] = h - ch[i];
    m[i] = 0;
  }

#if 0
  cerr <<"\n";cerr.flush();
  for(int i=0;i<nb_tests;i++){
    for(int j=0;j<nb_samples;j++){
      if(fe_get_bit(j,x[i])) cerr<<"1";
      else cerr << "0";
    }
    cerr <<"\n";cerr.flush();
  }
#endif

  int nb_uint32 = (nb_samples+31)/32;
  uint32_t *z00=new uint32_t[nb_uint32];
  uint32_t *z01=new uint32_t[nb_uint32];
  uint32_t *z10=new uint32_t[nb_uint32];
  uint32_t *z11=new uint32_t[nb_uint32];
  for(int n = 0; n < nb_selected; n++) {
    selected[n] = -1;
    double best_s = 0;
    for(int i = 0; i < nb_tests; i++) {
      //cerr << "S[i]: " << n << " " << i << " " << s[i] << " " << best_s << endl;cerr.flush();
      if(s[i]-best_s>FLOAT_TOL) {
	//cerr << "  m[i]: " << m[i] << " " << n << endl;cerr.flush();
        if(m[i] < n) {
          fe_and(nb_samples, x[i], y, z11);
          fe_and_not(nb_samples, x[i], y, z10);
          fe_and_not(nb_samples, y, x[i], z01);
          fe_and_not_not(nb_samples, x[i], y, z00);
          while((s[i]-best_s)>FLOAT_TOL && m[i] < n) {
            double h_y_xi_xmi = fe_logn[nb_samples] -
              (  fe_nlogn[fe_count_and    (nb_samples, z11, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and    (nb_samples, z10, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and    (nb_samples, z01, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and    (nb_samples, z00, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and_not(nb_samples, z11, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and_not(nb_samples, z10, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and_not(nb_samples, z01, x[selected[m[i]]])]
               + fe_nlogn[fe_count_and_not(nb_samples, z00, x[selected[m[i]]])])/double(nb_samples);
            double h_xi_xmi = fe_entropy_couple(nb_samples, x[i], x[selected[m[i]]]);
            double ss = ch[selected[m[i]]] - (h_y_xi_xmi - h_xi_xmi);
	    //cerr << "     S[i]: " << s[i] << " " << ss << endl;cerr.flush();
	    //cerr << "           " << h_y_xi_xmi << " " << h_xi_xmi << endl;cerr.flush();
            if(ss < s[i]) s[i] = ss;
            m[i]++;
          }
        }
        if(s[i]-best_s>FLOAT_TOL) {
          best_s = s[i];
          selected[n] = i;
        }
      }
    }
  }
  delete [] z00;
  delete [] z01;
  delete [] z10;
  delete [] z11;
  delete [] s;
  delete [] ch;
  delete [] m;
}

void fe_selection_mim(int nb_samples,
                      int nb_tests, uint32_t **x, uint32_t *y,
                      int nb_selected, int *selected) {

#ifdef DEBUG
  if(!fe_was_initialized) {
    cerr << "fe_init_tables() was not called!\n";
    abort();
  }
#endif

  if(nb_samples > 65535) {
    cerr << "Too many pictures, the nlogn table is too small.\n";
    exit(1);
  }

  Couple *tmp=new Couple[nb_tests];
  double h = fe_entropy(nb_samples, y);

  for(int i = 0; i < nb_tests; i++) {
    tmp[i].index = i;
    tmp[i].value = fe_entropy_couple(nb_samples, y, x[i]) - h - fe_entropy(nb_samples, x[i]);
  }

  qsort(tmp, nb_tests, sizeof(Couple), compare_couple);

  // Here we have the features sorted according to their mutual
  // information with the class to predict

  for(int n = 0; n < nb_selected; n++) selected[n] = tmp[n].index;
  delete [] tmp;
}

inline int fe_count_and_and(int n, uint32_t *a, uint32_t *b, uint32_t *c) {
  uint16_t *aa = (uint16_t *) a, *bb = (uint16_t *) b, *cc = (uint16_t *) c;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[*aa++ & *bb++ & *cc++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & *aa & *bb & *cc];
  return t;
}

inline int fe_count_and_not_and(int n, uint32_t *a, uint32_t *b, uint32_t *c) {
  uint16_t *aa = (uint16_t *) a,  *bb = (uint16_t *) b, *cc = (uint16_t *) c;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[*aa++ & ~*bb++ & *cc++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & *aa & ~*bb & *cc];
  return t;
}

inline int fe_count_and_not_not_and(int n, uint32_t *a, uint32_t *b, uint32_t *c) {
  uint16_t *aa = (uint16_t *) a, *bb = (uint16_t *) b, *cc = (uint16_t *) c;
  int t = 0;
  for(int k = 0; k < n/16; k++) t += fe_nb_bits[65535 & ~*aa++ & ~*bb++ & *cc++];
  if(n%16 > 0) t += fe_nb_bits[(65535 >> (16-(n%16))) & ~*aa & ~*bb & *cc];
  return t;
}

inline double fe_conditional_entropy(int n, uint32_t *a, uint32_t *c) {
  int nc = fe_count(n, c);
  int nac = fe_count_and(n, a, c);
  return fe_logn[nc] - (fe_nlogn[nac] + fe_nlogn[nc-nac])/double(nc);
}

inline double fe_conditional_entropy_couple(int n, uint32_t *a, uint32_t *b, uint32_t *c) {
  int nc = fe_count(n, c);
  int n11 = fe_count_and_and(n, a, b, c);
  int n10 = fe_count_and_not_and(n, a, b, c);
  int n01 = fe_count_and_not_and(n, b, a, c);
  int n00 = fe_count_and_not_not_and(n, a, b, c);

  return fe_logn[nc] - (fe_nlogn[n00] + fe_nlogn[n01] + fe_nlogn[n10] + fe_nlogn[n11])/double(nc);
}
