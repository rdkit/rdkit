#ifndef BOOST_NUMERIC_BINDINGS_MUMPS_464_CMUMPS_C_HPP
#define BOOST_NUMERIC_BINDINGS_MUMPS_464_CMUMPS_C_HPP

/*

   THIS FILE IS PART OF MUMPS VERSION 4.6.4
   This Version was built on Thu Jan 11 13:32:35 2007


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet, Hybrid scheduling for the parallel solution of linear
   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).

*/
/* $Id: cmumps_c.hpp 39206 2007-09-12 06:58:12Z karlmeerbergen $ */
/* Mostly written in march 2002 (JYL) */

// This file is modified by Karl Meerbergen for C++ users

/* Complex datatypes */
typedef struct {float r,i;} mumps_complex;

/* Next line defines CMUMPS_INT, CMUMPS_DOUBLE and CMUMPS_DOUBLE2 */
#include "cmumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    CMUMPS_INT sym, par, job;
    CMUMPS_INT comm_fortran;    /* Fortran communicator */
    CMUMPS_INT icntl[40];
    CMUMPS_DOUBLE2 cntl[5];
    CMUMPS_INT n;
   
    CMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    CMUMPS_INT nz; CMUMPS_INT *irn; CMUMPS_INT *jcn; CMUMPS_DOUBLE *a;
    /* Distributed entry */
    CMUMPS_INT nz_loc; CMUMPS_INT *irn_loc; CMUMPS_INT *jcn_loc; CMUMPS_DOUBLE *a_loc;
    /* Element entry */
    CMUMPS_INT nelt; CMUMPS_INT *eltptr; CMUMPS_INT *eltvar; CMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    CMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    CMUMPS_INT *sym_perm;
    /* column permutation */
    CMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    CMUMPS_DOUBLE *colsca; CMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    CMUMPS_DOUBLE *rhs, *rhs_sparse, *sol_loc;
    CMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    CMUMPS_INT nrhs, lrhs, nz_rhs, lsol_loc;
  CMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  CMUMPS_INT mblock, nblock, nprow, npcol;
    CMUMPS_INT info[40],infog[40];
    CMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    CMUMPS_INT deficiency; CMUMPS_DOUBLE * nullspace; CMUMPS_INT * mapping;
    /* Schur */
    CMUMPS_INT size_schur; CMUMPS_INT *listvar_schur; CMUMPS_DOUBLE *schur;
    /* Internal parameters */
    CMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
  } CMUMPS_STRUC_C;


extern "C" {
  void cmumps_c(CMUMPS_STRUC_C * cmumps_par);
}

#endif

