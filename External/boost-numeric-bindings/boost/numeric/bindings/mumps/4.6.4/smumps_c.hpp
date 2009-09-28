#ifndef BOOST_NUMERIC_BINDINGS_MUMPS_464_SMUMPS_C_HPP
#define BOOST_NUMERIC_BINDINGS_MUMPS_464_SMUMPS_C_HPP

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
/* $Id: smumps_c.hpp 39206 2007-09-12 06:58:12Z karlmeerbergen $ */
/* Mostly written in march 2002 (JYL) */

// This file is modified by Karl Meerbergen for C++ users

/* Complex datatypes */

/* Next line defines SMUMPS_INT, SMUMPS_DOUBLE and SMUMPS_DOUBLE2 */
#include <smumps_prec.h>
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    SMUMPS_INT sym, par, job;
    SMUMPS_INT comm_fortran;    /* Fortran communicator */
    SMUMPS_INT icntl[40];
    SMUMPS_DOUBLE2 cntl[5];
    SMUMPS_INT n;
   
    SMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    SMUMPS_INT nz; SMUMPS_INT *irn; SMUMPS_INT *jcn; SMUMPS_DOUBLE *a;
    /* Distributed entry */
    SMUMPS_INT nz_loc; SMUMPS_INT *irn_loc; SMUMPS_INT *jcn_loc; SMUMPS_DOUBLE *a_loc;
    /* Element entry */
    SMUMPS_INT nelt; SMUMPS_INT *eltptr; SMUMPS_INT *eltvar; SMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    SMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    SMUMPS_INT *sym_perm;
    /* column permutation */
    SMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    SMUMPS_DOUBLE *colsca; SMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    SMUMPS_DOUBLE *rhs, *rhs_sparse, *sol_loc;
    SMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    SMUMPS_INT nrhs, lrhs, nz_rhs, lsol_loc;
  SMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  SMUMPS_INT mblock, nblock, nprow, npcol;
    SMUMPS_INT info[40],infog[40];
    SMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    SMUMPS_INT deficiency; SMUMPS_DOUBLE * nullspace; SMUMPS_INT * mapping;
    /* Schur */
    SMUMPS_INT size_schur; SMUMPS_INT *listvar_schur; SMUMPS_DOUBLE *schur;
    /* Internal parameters */
    SMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
  } SMUMPS_STRUC_C;


extern "C" {
void smumps_c(SMUMPS_STRUC_C * smumps_par);
}

#endif

