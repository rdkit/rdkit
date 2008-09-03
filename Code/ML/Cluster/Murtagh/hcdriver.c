/* hcdriver.f -- translated by f2c (version 20020208).

   hand-modified by gL, but it's still ugly as sin

*/
#include <stdlib.h>

#include "f2c.h"


int distdriver_(n, len, d__, iopt, ia, ib, crit)
integer *n, *len;
doublereal *d__;
integer *iopt, *ia, *ib;
doublereal *crit;
{
  static logical *flag__;
  static doublereal *membr,*disnn;
  extern /* Subroutine */ int hc_();
  static integer *nn;
  flag__ = (logical *)malloc(*n*sizeof(logical));
  membr = (doublereal *)malloc(*n*sizeof(doublereal));
  disnn = (doublereal *)malloc(*n*sizeof(doublereal));
  nn = (integer *)malloc(*n*sizeof(integer));

  /* Parameter adjustments */
  --crit;
  --ib;
  --ia;
  --d__;

  /* Function Body */
  hc_(n, len, iopt, &ia[1], &ib[1], &crit[1], membr, nn, disnn, flag__, &
      d__[1]);

  free(membr);
  free(nn);
  free(disnn);
  free(flag__);
  return 0;
} /* distdriver_ */

