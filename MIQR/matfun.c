/*---------------------------------------------------------------------------*
 * matfun.c of MQR
 *---------------------------------------------------------------------------*
 * Na Li, Jul 12, 2004                                                       *
 *                                                                           *
 * Report bugs / send comments to: saad@cs.umn.edu, nli@cs.umn.edu           *
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heads.h"
#include "defs.h"
#include "matfun.h"

int readhbc_c( matptr mat, io_t *pio, int manu_rhs )
{
    int job, ncol, nrow, nrhs, ierr, i, k, flag;
    int *ia = NULL, *ja = NULL;
    double *a = NULL, *rhs = NULL;
    char guesol[3];
    int nnz, nmax, nzmax, rhs_sz;

/* find out size of Harwell-Boeing matrix ---------------------------*/
    job = 0;
    nmax = nzmax = 1;
    nrhs = 0;
    readmtc( &nmax, &nzmax, &job, pio->Fname, a, ja, ia, rhs, &nrhs,
             guesol, &nrow, &ncol, &nnz,
             pio->title, pio->key, pio->type, &ierr );
    if( ierr != 0 ) {
        fprintf( stderr, "readhb: err in read matrix header = %d\n", ierr );
        return ierr;
    }
/* some consistency checks ------------------------------------------*/
    pio->mdim = nrow;
    pio->ndim = ncol;
    if( nrow < ncol ) {
        fprintf( stderr, "readhb: matrix size m x n must satisfy m >= n\n" );
        fprintf( stderr, "        for this version.\n" );
        return -1;
    }
/* allocate space ---------------------------------------------------*/
    job = 2 + nrhs;
    if( nrhs > 0 ) {
      rhs_sz = nrow * nrhs;
      if( guesol[0] == 'G' || guesol[0] == 'g' ) rhs_sz += ncol * nrhs;
      if( guesol[1] == 'X' || guesol[1] == 'x' ) rhs_sz += ncol * nrhs;
      rhs = (double *)Malloc( sizeof(double)*rhs_sz, "readhb" );
      nrhs = rhs_sz;
    }
    ia     = (int *)Malloc( sizeof(int)*(ncol+1), "readhb" );
    ja     = (int *)Malloc( sizeof(int)*nnz, "readhb" );
    a      = (double *)Malloc( sizeof(double)*nnz, "readhb" );
/* read matrix ------------------------------------------------------*/
    nmax = ncol+1;
    nzmax = nnz;
    readmtc( &nmax, &nzmax, &job, pio->Fname, a, ja, ia, rhs, &nrhs,
             guesol, &nrow, &ncol, &nnz,
             pio->title, pio->key, pio->type, &ierr );
    if( ierr != 0 ) {
      fprintf( stderr, "readhb: err in read matrix data = %d\n", ierr );
      return ierr;
    }
    if( job <= 1 ) {
      fprintf( stderr, "readhb: matrix values not available...\n" );
      return -1;
    }
    if( job <= 2 ) {
      nrhs = 0;
    }

/*---------------------------------------------------------------------- 
|	convert (fortran sytle) csc matrix into (c style) SparMat format
-----------------------------------------------------------------------*/
    flag = MAT_COLBASE_INDEX|MAT_COLBASE_VALUE|
           MAT_ROWBASE_INDEX|MAT_ROWBASE_VALUE;
    if( ( ierr = CSC2Mat( nrow, ncol, a, ja, ia, mat, flag ) ) != 0 ) {
        fprintf( stderr, "readhb: CSC2mat error\n" );
        return ierr;
    }
    pio->mdim = nrow;
    pio->ndim = ncol;
    pio->nnz  = nnz;
    pio->title[72] = '\0';
    pio->key[8] = '\0';
    pio->type[3] = '\0';
    pio->nrhs = nrhs;

    if( !manu_rhs && nrhs > 0 ) {
      int offset = 0;
      pio->rhs = (double *)Malloc( nrhs*nrow*sizeof(double), "readhb" );
      pio->guess = (double *)Malloc( nrhs*ncol*sizeof(double), "readhb" );

      memcpy( pio->rhs, rhs, nrhs*nrow*sizeof(double) );
      offset += nrhs*nrow;
      if( guesol[0] == 'G' || guesol[0] == 'g' ) {
        memcpy( pio->guess, rhs + offset, nrhs*ncol*sizeof(double) );
        offset += nrhs * ncol;
      } else {
        memset( pio->guess, 0, nrhs*ncol*sizeof(double) );
      }
      if( guesol[1] == 'X' || guesol[1] == 'x' ) {
        pio->sol = (double *)Malloc( nrhs*ncol*sizeof(double), "readhb" );
        memcpy( pio->sol, rhs + offset, nrhs*ncol*sizeof(double) );
      }
    } else {
      pio->nrhs = nrhs = 1;
      pio->sol = (double *)Malloc( nrhs*ncol*sizeof(double), "readhb" );
      pio->guess = (double *)Malloc( nrhs*ncol*sizeof(double), "readhb" );
      pio->rhs = (double *)Malloc( nrhs*nrow*sizeof(double), "readhb" );
      for( k = 0; k < nrhs; k++ ) {
        double *sol = pio->sol + ncol*k;
        double *guess = pio->guess + ncol*k;
        double *rhs = pio->rhs + nrow*k;
        for( i = 0; i < ncol; i++ ) {
          sol[i] = 1.0;
          guess[i] = 0.0; /* initial guess = 0 */
        }
        matxvec( mat, sol, rhs );
      }
    }

    free( a );    
    free( ja );
    free( ia );
    if( rhs ) free( rhs );
    
    return 0;
}

int setupMat( matptr mat, int m, int n, int flag, int bufsz, int incsz )
{
/*----------------------------------------------------------------------
| Initialize SparCol structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( mat )   =  Pointer to a SparMat struct.
|       m   =  the # of rows
|       n   =  the # of columns
| CSR_INDEX(flag): initialize row-based format
| CSC_INDEX(flag): initialize column-based format
|
| On return:
|===========
|
|  mat->m, n
|      ->rvec if CSR_INDEX(flag)
|      ->cvec if CSC_INDEX(flag)
|
| integer value returned:
|             0   --> successful return.
|--------------------------------------------------------------------*/
   mat->m = m;
   mat->n = n;
   mat->nnz = 0;
   mat->flag = flag;
   mat->bufsz = bufsz;
   mat->incsz = incsz;
   if( CSR_INDEX(flag) ) {
     mat->rvec = (int *)Malloc( (m+1)*sizeof(int), "setupMat" );
     mat->ridx = (int *)Malloc( bufsz*sizeof(int), "setupMat" );
     if( CSR_VALUE(flag) ) {
       mat->rval = (double *)Malloc( bufsz*sizeof(double), "setupMat" );
     } else {
       mat->rval = NULL;
     }
   } else {
     mat->rvec = NULL;
     mat->ridx = NULL;
     mat->rval = NULL;
   }
   if( CSC_INDEX(flag) ) {
     mat->cvec = (int *)Malloc( (n+1)*sizeof(int), "setupMat" );
     mat->cidx = (int *)Malloc( bufsz*sizeof(int), "setupMat" );
     if( CSC_VALUE(flag) ) {
       mat->cval = (double *)Malloc( bufsz*sizeof(double), "setupMat" );
     } else {
       mat->cval = NULL;
     }
   } else {
     mat->cvec = NULL;
     mat->cidx = NULL;
     mat->cval = NULL;
   }
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCSC
|--------------------------------------------------------------------*/

int cleanMat( matptr mat )
{
/*----------------------------------------------------------------------
| Free up memory allocated for SparMat structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( mat )  =  Pointer to a SparMat struct.
|--------------------------------------------------------------------*/
  if (mat == NULL) return 0;
  if( mat->rvec != NULL ) free( mat->rvec );
  if( mat->ridx != NULL ) free( mat->ridx );
  if( mat->rval != NULL ) free( mat->rval );
  if( mat->cvec != NULL ) free( mat->cvec );
  if( mat->cidx != NULL ) free( mat->cidx );
  if( mat->cval != NULL ) free( mat->cval );
  free(mat);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanMat
|--------------------------------------------------------------------*/

int checkMatBuffer( matptr mat, int nadd )
{
  if( mat->bufsz < mat->nnz + nadd ) {
    int sz = mat->bufsz;
    while( sz < mat->nnz + nadd ) {
      sz += mat->incsz;
    }
    if( CSR_INDEX(mat->flag) ) {
      mat->ridx = (int *)Realloc( mat->ridx, sz*sizeof(int), "chkBuf" );
    }
    if( CSR_VALUE(mat->flag) ) {
      mat->rval = (double *)Realloc( mat->rval, sz*sizeof(double), "chkBuf" );
    }
    if( CSC_INDEX(mat->flag) ) {
      mat->cidx = (int *)Realloc( mat->cidx, sz*sizeof(int), "chkBuf" );
    }
    if( CSC_VALUE(mat->flag) ) {
      mat->cval = (double *)Realloc( mat->cval, sz*sizeof(double), "chkBuf" );
    }
    mat->bufsz = sz;
  }
  return 0;
}

int loadVectorR( vecptr v, matptr mat, int index )
{
  int start = mat->rvec[index], end = mat->rvec[index+1];
  v->nnz = end - start;
  v->idx = mat->ridx + start;
  v->val = mat->rval + start;
  return 0;
} 

int loadVectorC( vecptr v, matptr mat, int index )
{
  int start = mat->cvec[index], end = mat->cvec[index+1];
  v->nnz = end - start;
  v->idx = mat->cidx + start;
  v->val = mat->cval + start;
  return 0;
} 

int printMat( matptr mat, char *name, FILE *fp )
{
  int m = mat->m, n = mat->n, i, j;
  Vector v;
  fprintf( fp, "********** %s *************\n", name );
  if( CSR_VALUE(mat->flag) ) {
    for( i = 0; i < m; i++ ) {
      loadVectorR( &v, mat, i );
      fprintf( fp, "Row %d: nnz = %d\n  ", i, v.nnz );
      for( j = 0; j < v.nnz; j++ ) {
        fprintf( fp, "(%d,%g) ", v.idx[j], v.val[j] );
        if( (j+1) % 6 == 0 ) fprintf( fp, "\n  " );
      }
      if( j % 6 != 0 ) fprintf( fp, "\n" );
    }
  } else if( CSC_VALUE(mat->flag) ) {
    for( i = 0; i < n; i++ ) {
      loadVectorC( &v, mat, i );
      fprintf( fp, "Column %d: nnz = %d\n  ", i, v.nnz );
      for( j = 0; j < v.nnz; j++ ) {
        fprintf( fp, "(%d,%g) ", v.idx[j], v.val[j] );
        if( (j+1) % 6 == 0 ) fprintf( fp, "\n  " );
      }
      if( j % 6 != 0 ) fprintf( fp, "\n" );
    }
  } else if( CSR_INDEX(mat->flag) ) {
    for( i = 0; i < m; i++ ) {
      loadVectorR( &v, mat, i );
      fprintf( fp, "Row %d: nnz = %d\n  ", i, v.nnz );
      for( j = 0; j < v.nnz; j++ ) {
        fprintf( fp, "(%d) ", v.idx[j] );
        if( (j+1) % 20 == 0 ) fprintf( fp, "\n  " );
      }
      if( j % 20 != 0 ) fprintf( fp, "\n" );
    }
  } else if( CSC_INDEX(mat->flag) ) {
    for( i = 0; i < n; i++ ) {
      loadVectorC( &v, mat, i );
      fprintf( fp, "Column %d: nnz = %d\n  ", i, v.nnz );
      for( j = 0; j < v.nnz; j++ ) {
        fprintf( fp, "(%d) ", v.idx[j] );
        if( (j+1) % 20 == 0 ) fprintf( fp, "\n  " );
      }
      if( j % 20 != 0 ) fprintf( fp, "\n" );
    }
  }
  return 0;
}

int nnzMat( matptr mat )
{
  if( mat == NULL ) return 0;
  return mat->nnz;
}

int setupMQR( mqrptr mqr, int n )
{
/*----------------------------------------------------------------------
| Initialize ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|  ( mqr )  =  Pointer to a SparMQR struct
|
| On return:
|===========
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    mqr->n = n;
    mqr->s = 0;
    mqr->perm = (int *)Malloc( n*sizeof(int), "setupMQR" );
    mqr->iperm = (int *)Malloc( n*sizeof(int), "setupMQR" );
    mqr->D_1 = NULL;
    mqr->E = NULL;
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupMQR
|--------------------------------------------------------------------*/

int cleanMQR( mqrptr mqr )
{
/*----------------------------------------------------------------------
| Free up memory allocated for SparMQR structs.
|----------------------------------------------------------------------
| on entry:
|==========
|  ( mqr )  =  Pointer to a ILUSpar struct.
|--------------------------------------------------------------------*/
  if( mqr == NULL ) return 0;
  if( mqr->perm ) free( mqr->perm );
  if( mqr->iperm ) free( mqr->iperm );
  if( mqr->D_1 ) free( mqr->D_1 );
  cleanMat( mqr->E );
  free( mqr );
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanMQR
|--------------------------------------------------------------------*/

int nnzMQR( mqrptr mqr )
{
  if( mqr == NULL ) return 0;
  return mqr->s + nnzMat(mqr->E);
}

int setupMQRS( mqrsptr mqrs )
{
  int i;
  mqrs->n = 0;
  mqrs->lvl = 0;
  for( i = 0; i < MAX_MQR_LEVEL; i++ ) {
    mqrs->amqrs[i] = NULL;
  }
  mqrs->R = (matptr)Malloc( sizeof(SparMat), "setupMQRS" );
  memset( mqrs->R, 0, sizeof(SparMat) );

  return 0;
}

int cleanMQRS( mqrsptr mqrs )
{
  int i;
  if( mqrs == NULL ) return 0;

  for( i = 0; i < MAX_MQR_LEVEL; i++ ) {
    cleanMQR( mqrs->amqrs[i] );
  }

  cleanMat( mqrs->R );
  free( mqrs );

  return 0;
}

int nnzMQRS( mqrsptr mqrs )
{
  int i, nz = 0;
  if( mqrs == NULL ) return 0;
  for( i = 0; i < mqrs->lvl; i++ ) {
    nz += nnzMQR( mqrs->amqrs[i] );
  }
  nz += nnzMat( mqrs->R );
  return nz;
}

int printMQRS( mqrsptr mqrs, char *name, FILE *fp )
{
  int i;
  mqrptr mqr;
  fprintf( fp, "=============== %s ==============\n", name );
  for( i = 0; i < mqrs->lvl; i++ ) {
    mqr = mqrs->amqrs[i];
    fprintf( fp, "%s, lvl %d:\n", name, i );
    vec_output( mqr->s, mqr->D_1, "inv(D)", fp );
    printMat( mqr->E, "E", fp );
  }
  printMat( mqrs->R, "R", fp );
  fprintf( fp, "========== end of %s ============\n", name );
  return 0;
}

int CSC2Mat(int m, int n, double *a, int *ja, int *ia, matptr mat, int flag)
{
/*----------------------------------------------------------------------
| Convert CSC matrix to SparMat struct
|----------------------------------------------------------------------
| on entry:
|==========
|         m  = the # of rows
|         n  = the # of columns
| a, ja, ia  = Matrix stored in CSC format (with FORTRAN indexing).
| (flag & MAT_ROWBASE_INDEX) != 0: setup row-based structure
| (flag & MAT_ROWBASE_VALUE) != 0: setup row-based data
| (flag & MAT_COLBASE_INDEX) != 0: setup col-based structure
| (flag & MAT_COLBASE_VALUE) != 0: setup col-based data
|
| On return:
|===========
|
| ( mat )  =  Matrix stored as SparMat struct.
| mat->ridx is in increasing order in each row
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int row, col, i, j, j1, j2, id, nnz;

  j1 = ia[0] - 1;
  j2 = ia[n] - 1;
  nnz = j2 - j1;
  setupMat( mat, m, n, flag, nnz, nnz );
  mat->nnz = nnz;

  if( CSR_INDEX(flag) ) {
    for( i = 0; i <= m; i++ ) {
      mat->rvec[i] = 0;
    }
    for( j = j1; j < j2; j++ ) {
      row = ja[j] - 1;
      mat->rvec[row+1]++;
    }
    for( i = 0; i < m; i++ ) {
      mat->rvec[i+1] += mat->rvec[i];
    }
  }
  if( CSC_INDEX(flag) ) {
    mat->cvec[0] = 0;
    for( i = 1; i <= n; i++ ) {
      mat->cvec[i] = ia[i] - ia[0];
    }
  }

  id = 0;
  for( col = 0; col < n; col++ ) {
    j1 = ia[col] - 1;
    j2 = ia[col+1] - 1;
    for( j = j1; j < j2; j++ ) {
      row = ja[j] - 1;
      if( CSR_INDEX(flag) ) {
        mat->ridx[mat->rvec[row]] = col;
        if( CSR_VALUE(flag) ) {
          mat->rval[mat->rvec[row]] = a[j];
        }
        mat->rvec[row]++;
      }
      if( CSC_INDEX(flag) ) {
        mat->cidx[id] = row;
        if( CSC_VALUE(flag) ) {
          mat->cval[id] = a[j];
        }
        id++;
      }
    }
  }
  if( CSR_INDEX(flag) ) {
    for( i = m; i > 0; i-- ) {
      mat->rvec[i] = mat->rvec[i-1];
    }
    mat->rvec[0] = 0;
  }

  return 0;
}

int transMat( matptr mat )
{
/*---------------------------------------------------------------------
| This function does the matrix tranpose
|----------------------------------------------------------------------
| on entry:
| mat  = a matrix (in SparMat format) stored column wise
|
| on return
| mat  = the matrix stored row wise as well
| integer value returned:
|             0   --> successful return.
|--------------------------------------------------------------------*/
  int m = mat->m, n = mat->n;
  int i, j, j1, j2, row, col;
  int *cvec = mat->cvec, *rvec = mat->rvec;
  int *cidx = mat->cidx, *ridx = mat->ridx;
  double *cval = mat->cval, *rval = mat->rval;
  j1 = cvec[0];
  j2 = cvec[n];
  if( CSR_INDEX(mat->flag) ) {
    for( i = 0; i <= m; i++ ) {
      rvec[i] = 0;
    }
    for( j = j1; j < j2; j++ ) {
      row = cidx[j];
      rvec[row+1]++;
    }
    for( i = 0; i < m; i++ ) {
      rvec[i+1] += rvec[i];
    }
    if( CSR_VALUE(mat->flag) ) {
      for( col = 0; col < n; col++ ) {
        j1 = cvec[col];
        j2 = cvec[col+1];
        for( j = j1; j < j2; j++ ) {
          row = cidx[j];
          ridx[rvec[row]] = col;
          rval[rvec[row]] = cval[j];
          rvec[row]++;
        }
      }
    } else {
      for( col = 0; col < n; col++ ) {
        j1 = cvec[col];
        j2 = cvec[col+1];
        for( j = j1; j < j2; j++ ) {
          row = cidx[j];
          ridx[rvec[row]] = col;
          rval[rvec[row]] = cval[j];
          rvec[row]++;
        }
      }
    }
    for( i = m; i > 0; i-- ) {
      rvec[i] = rvec[i-1];
    }
    rvec[0] = 0;
  }
  return 0;
}

/*---------------------------------------------------------------------
|     end of CSC2Mat
|--------------------------------------------------------------------*/

int matxvec( matptr mat, double *x, double *y )
{
/*---------------------------------------------------------------------
| This function does the matrix vector product y = A x.
|----------------------------------------------------------------------
| on entry:
| mat  = the matrix (in SparMat format)
| x    = a vector, size = mat->n
|
| on return
| y     = the product of A * x, size = mat->m
| integer value returned:
|             0   --> successful return.
|--------------------------------------------------------------------*/
  int m = mat->m, n = mat->n, i, j, j1, j2;
  int *idx;
  double *val, xval, yval;

  for( i = 0; i < m; i++ ) y[i] = 0.0;
  if( CSR_VALUE(mat->flag) ) {
    idx = mat->ridx;
    val = mat->rval;
    for( i = 0; i < m; i++ ) {
      j1 = mat->rvec[i];
      j2 = mat->rvec[i+1];
      yval = y[i];
      for( j = j1; j < j2; j++ ) {
        yval += val[j] * x[idx[j]];
      }
      y[i] = yval;
    }
  } else if( CSC_VALUE(mat->flag) ) {
    idx = mat->cidx;
    val = mat->cval;
    for( i = 0; i < n; i++ ) {
      j1 = mat->cvec[i];
      j2 = mat->cvec[i+1];
      xval = x[i];
      for( j = j1; j < j2; j++ ) {
        y[idx[j]] += val[j] * xval;
      }
    }
  } else {
    return -1;
  }
  return 0;
}

int matTxvec( matptr mat, double *x, double *y )
{
/*---------------------------------------------------------------------
| This function does the matrix vector product y = A^T * x.
|----------------------------------------------------------------------
| on entry:
| mat  = the matrix (in SparMat format)
| x    = a vector, size = mat->m
|
| on return
| y     = the product of A^T * x, size = mat->n
| integer value returned:
|             0   --> successful return.
|--------------------------------------------------------------------*/
  int m = mat->m, n = mat->n, i, j, j1, j2;
  int *idx;
  double *val, xval, yval;

  for( i = 0; i < n; i++ ) y[i] = 0.0;
  if( CSR_VALUE(mat->flag) ) {
    idx = mat->ridx;
    val = mat->rval;
    for( i = 0; i < m; i++ ) {
      j1 = mat->rvec[i];
      j2 = mat->rvec[i+1];
      xval = x[i];
      for( j = j1; j < j2; j++ ) {
        y[idx[j]] += val[j] * xval;
      }
    }
  } else if( CSC_VALUE(mat->flag) ) {
    idx = mat->cidx;
    val = mat->cval;
    for( i = 0; i < n; i++ ) {
      j1 = mat->cvec[i];
      j2 = mat->cvec[i+1];
      yval = y[i];
      for( j = j1; j < j2; j++ ) {
        yval += val[j] * x[idx[j]];
      }
      y[i] = yval;
    }
  } else {
    return -1;
  }
  return 0;
}

int mqrsol( mqrsptr mqrs, double *y, double *x, int transp, double *wk )
{
  /* transp = 0: Solve (MQR) * x = y */
  /* transp = 1: Solve (MQR)^T * x = y */
  /* the first n elements are inv of diags */
  int i, j, n = mqrs->n, lvl = mqrs->lvl;
  matptr R = mqrs->R;
  mqrptr mqr = NULL;
  double xval, *diag = R->cval, *xPosR, *xPosI, *xPosI1; 
  int *pos = mqrs->pos;
  Vector v;

  if( x != y ) {
    for( i = 0; i < n; i++ ) x[i] = y[i];
  }

  if( transp ) {
    for( i = 0; i < lvl; i++ ) {
      xPosI = x + pos[i];
      xPosI1 = x + pos[i+1];
      mqr = mqrs->amqrs[i];
      vec_perm( mqr->n, xPosI, mqr->perm, wk );
      for( j = 0; j < mqr->s; j++ ) {
        xPosI[j] = wk[j] * mqr->D_1[j];
      }
      matTxvec( mqr->E, xPosI, xPosI1 );
      vec_sub( mqr->n-mqr->s, wk+mqr->s, xPosI1, xPosI1 );
    }
    if( R->n > 0 ) {
      xPosR = x + pos[lvl];
      for( i = 0; i < R->n; i++ ) {
        loadVectorC( &v, R, i );
        xval = xPosR[i];
        /* v->nnz must be greater 1 and the last element is the diagonal */
        for( j = 0; j < v.nnz; j++ ) {
          xval -= v.val[j] * xPosR[v.idx[j]];
        }
        xPosR[i] = xval * diag[i]; /* diag = inv(Rii) */
      }
    }
  } else {
    if( R->n > 0 ) {
      xPosR = x + pos[lvl];
      for( i = R->n-1; i >= 0; i-- ) {
        loadVectorC( &v, R, i );
        xPosR[i] *= diag[i]; /* diag = inv(Rii) */
        xval = xPosR[i];
        for( j = 0; j < v.nnz; j++ ) {
          xPosR[v.idx[j]] -= v.val[j] * xval;
        }
      }
    }
    for( i = lvl-1; i >= 0; i-- ) {
      xPosI = x + pos[i];
      xPosI1 = x + pos[i+1];
      mqr = mqrs->amqrs[i];
      matxvec( mqr->E, xPosI1, wk );
      vec_sub( mqr->s, xPosI, wk, wk );
      for( j = 0; j < mqr->s; j++ ) {
        wk[j] *= mqr->D_1[j];
      }
      vec_cpy( mqr->n-mqr->s, xPosI1, wk+mqr->s );
      vec_perm( mqr->n, wk, mqr->iperm, xPosI );
    }
  }
  return 0;
}

double vec_nrminf( int n, double *v )
{
  int i;
  double nrm = 0.0;
  for( i = 0; i < n; i++ ) {
    if( fabs( v[i] ) > nrm ) nrm = fabs( v[i] );
  }
  return nrm;
}

double vec_nrm2( int n, double *v )
{
  int i;
  double val, nrm = 0.0;
  for( i = 0; i < n; i++ ) {
    val = v[i];
    nrm += val * val;
  }
  return sqrt(nrm);
}

double vec_nrm2_diff( int n, double *v1, double *v2 )
{
  /* Calculate || v1 - v2 || */
  int i;
  double diff, nrm = 0.0;
  for ( i = 0; i < n; i++ ) {
    diff = v1[i] - v2[i];
    nrm += diff * diff;
  }
  return sqrt(nrm);
}

int vec_add( int n, double *v1, double *v2, double *x )
{
  int i;
  for( i = 0; i < n; i++ ) x[i] = v1[i] + v2[i];
  return 0;
}

int vec_sub( int n, double *v1, double *v2, double *x )
{
  int i;
  for( i = 0; i < n; i++ ) x[i] = v1[i] - v2[i];
  return 0;
}

int vec_xadd( int n, double *v1, double *v2, double alpha, double *x )
{
  /* x = v1 + v2 * alpha, x could be either v1 or v2 */
  int i;
  for( i = 0; i < n; i++ ) x[i] = v1[i] + alpha * v2[i];
  return 0;
}

int vec_cpy( int n, double *v1, double *v2 )
{
  /* v2 = v1 */
  int i;
  for( i = 0; i < n; i++ ) v2[i] = v1[i];
  return 0;
}

int vec_perm( int n, double *v1, int *perm, double *v2 )
{
  /* v2 = perm(v1) */
  int i;
  for( i = 0; i < n; i++ ) {
    v2[i] = v1[perm[i]];
  }
  return 0;
}

double vec_dot( int n, double *v1, double *v2 )
{
  /* return <v1,v2> */
  int i;
  double ret = 0.0;
  for( i = 0; i < n; i++ ) ret += v1[i] * v2[i];
  return ret;
}

int vec_output( int n, double *v, char *name, FILE *fp )
{
  /* output vector v to file pointer fp */
  int i;
  fprintf( fp, "============ %s ==========\n", name );
  for( i = 0; i < n; i++ ) fprintf( fp, "%s[%d] = %g\n", name, i, v[i] );
  return 0;
}

double sparvec_nrm1( vecptr v )
{
  int i;
  double  nrm = 0.0;
  for( i = 0; i < v->nnz; i++ ) {
    nrm += fabs(v->val[i]);
  }
  return nrm;
}

double sparvec_nrm2( vecptr v )
{
  int i;
  double val, nrm = 0.0;
  for( i = 0; i < v->nnz; i++ ) {
    val = v->val[i];
    nrm += val * val;
  }
  return sqrt( nrm );
}

int fsparvec_add( fvecptr v1, vecptr v2, double alpha )
{
  /* v1 = v1 + v2 * alpha */
  int i, idx;
  double val;
  for( i = 0; i < v2->nnz; i++ ) {
    idx = v2->idx[i];
    val = v2->val[i];
    if( v1->pat[idx] == 1 ) {
      v1->val[idx] += val * alpha;
    } else {
      /* fillin */
      v1->val[idx] = alpha * val;
      v1->pat[idx] = 1;
      v1->idx[v1->nnz] = idx;
      v1->nnz++;
    }
  }
  return 0;
}

int fsparvec_iadd( fvecptr v1, vecptr v2, double alpha, double tol )
{
  /* incomplete version of fsparvec_add */
  int i, idx;
  double val, valadd;
  for( i = 0; i < v2->nnz; i++ ) {
    idx = v2->idx[i];
    val = v2->val[i];
    valadd = alpha * val;
    if( v1->pat[idx] == 1 ) {
      v1->val[idx] += valadd;
    } else if( fabs(valadd) > tol ) {
      /* fillin */
      v1->val[idx] = alpha * val;
      v1->pat[idx] = 1;
      v1->idx[v1->nnz] = idx;
      v1->nnz++;
    }
  }
  return 0;
}

double fsparvec_nrm1( fvecptr v )
{
  int i, nnz = v->nnz;
  double nrm = 0.0;
  for( i = 0; i < nnz; i++ ) {
    nrm += fabs( v->val[v->idx[i]] );
  }
  return sqrt( nrm );
}

double fsparvec_nrm2( fvecptr v )
{
  int i, nnz = v->nnz;
  double nrm = 0.0, val;
  for( i = 0; i < nnz; i++ ) {
    val = v->val[v->idx[i]];
    nrm += val * val;
  }
  return sqrt( nrm );
}

double fsparvec_dot( fvecptr v1, vecptr v2 )
{
  int i, idx;
  double retval = 0.0;
  for( i = 0; i < v2->nnz; i++ ) {
    idx = v2->idx[i];
    if( v1->pat[idx] == 1 ) {
      retval += v1->val[idx] * v2->val[i];
    }
  }
  return retval;
}

int fsparvec_loadC(  fvecptr v1, matptr mat, int index )
{
  /* load v2 to v1 */
  int j1, j2, j, idx;
  j1 = mat->cvec[index];
  j2 = mat->cvec[index+1];
  v1->dim = mat->m;
  v1->nnz = 0;
  for( j = j1; j < j2; j++ ) {
    idx = mat->cidx[j];
    v1->idx[v1->nnz] = idx;
    v1->pat[idx] = 1;
    v1->val[idx] = mat->cval[j];
    v1->nnz++;
  }
  return 0;
}

int fsparvec_loadR(  fvecptr v1, matptr mat, int index )
{
  /* load v2 to v1 */
  int j1, j2, j, idx;
  j1 = mat->rvec[index];
  j2 = mat->rvec[index+1];
  v1->dim = mat->n;
  v1->nnz = 0;
  for( j = j1; j < j2; j++ ) {
    idx = mat->ridx[j];
    v1->idx[v1->nnz] = idx;
    v1->pat[idx] = 1;
    v1->val[idx] = mat->rval[j];
    v1->nnz++;
  }
  return 0;
}

