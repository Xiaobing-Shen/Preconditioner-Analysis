/*---------------------------------------------------------------------------*
 * Preconditioned CGNR                                                       *
 *---------------------------------------------------------------------------*
 * Na Li, Mar, 2005                                                          *
 *                                                                           *
 * Report bugs / send comments to: nli@cs.umn.edu                            *
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heads.h"
#include "matfun.h"

int pcgnr( matptr mat, mqrsptr mqrs, double *rhs, double *sol,
           double tol, int *itmax, int flag )
{
/*----------------------------------------------------------------------
|                 *** Preconditioned CGNR ***
+-----------------------------------------------------------------------
| N. L. Jul 2004
+-----------------------------------------------------------------------
| on entry:
|========== 
|
|(mat)    = matrix in SparMat format.
|
|(mqrs)   = MQR factorization of input matrix from IMQR
|                                                                      
| rhs     = real vector of length m containing the right hand side.
|           
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| (itmax) = max number of iterations allowed. 
| flag    = 0: normal case
|           1: implements the no preconditioning option here. useful
|              for test purposes.
|
| on return:
|==========
|          int =  0 --> successful return.
|          int =  1 --> convergence not achieved in itmax iterations.
|
| sol   == contains an approximate solution (upon successful return).
| itmax == has changed. It now contains the number of steps required
|          to converge -- 
+---------------------------------------------------------------------*/
  int retval = 1;
  int m = mat->m, n = mat->n, maxits = *itmax, its = 0;
  double *r = NULL, *rbar = NULL, *z = NULL, *p = NULL, *w = NULL;
  double alpha, beta, scale, tolnrm2;
  double *wk = NULL; /* work array */

  r = (double *)Malloc( m * sizeof(double), "pcgnr" );
  rbar = (double *)Malloc( n * sizeof(double), "pcgnr" );
  z = (double *)Malloc( n * sizeof(double), "pcgnr" );
  p = (double *)Malloc( n * sizeof(double), "pcgnr" );
  w = (double *)Malloc( m * sizeof(double), "pcgnr" );
  wk = (double *)Malloc( n * sizeof(double), "pcgnr" );

  matxvec( mat, sol, r );
  vec_sub( m, rhs, r, r );  /* r = b - Ax */
  matTxvec( mat, r, rbar );       /* rbar = A^T * r */
  tolnrm2 = tol * tol * vec_dot( n, rbar, rbar );

  if( flag != 0 ) {
    vec_cpy( n, rbar, z );
  } else {
    /* z = inv( (MQR)^T * (MQR) ) * rbar */
    mqrsol( mqrs, rbar, z, 1, wk );   /* (MQR)^T * z = rbar */
    mqrsol( mqrs, z, z, 0, wk );      /* (MQR)   * z = z    */
  }
  vec_cpy( n, z, p );             /* p = z */
  scale = vec_dot( n, z, rbar );
  alpha = 1.0e+6;
  its++;
  
  while( its < maxits ) {
    matxvec( mat, p, w );
    alpha = scale / vec_dot( m, w, w );
    vec_xadd( n, sol, p, alpha, sol );
    vec_xadd( m, r, w, -alpha, r );
    matTxvec( mat, r, rbar );
    if( vec_dot( n, rbar, rbar ) < tolnrm2 ) {
      retval = 0;
      break;
    }
    if( flag != 0 ) {
      vec_cpy( n, rbar, z );
    } else {
      /* z = inv( (MQR)^T * (MQR) ) * rbar */
      mqrsol( mqrs, rbar, z, 1, wk );   /* (MQR)^T * z = rbar */
      mqrsol( mqrs, z, z, 0, wk );      /* (MQR)   * z = z    */
    }
    beta = 1.0 / scale;
    scale = vec_dot( n, z, rbar );
    beta = beta * scale;
    vec_xadd( n, z, p, beta, p );
    its++;
  }

  *itmax = its; 

  free( r );
  free( rbar );
  free( z );
  free( p );
  free( w );
  free( wk );

  return retval;
}
/*-------------------- end of pcgnr ------------------------------------
----------------------------------------------------------------------*/
