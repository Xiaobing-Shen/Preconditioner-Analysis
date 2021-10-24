/*-----------------------------------------------------------------------*
 * Multilevel Incomplete QR Factorization                                *
 *-----------------------------------------------------------------------*
 * Na Li, Mar, 2005                                                      *
 *                                                                       *
 * Report bugs / send comments to: nli@cs.umn.edu                        *
 *-----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "heads.h"
#include "matfun.h"

#define MIN_NUM (1e-14)

/* working arrays: */
double *vecnval = NULL, *vecmval = NULL;
int *vecnpat = NULL, *vecnidx = NULL;
int *vecmpat = NULL, *vecmidx = NULL;
grphnd nodes = NULL;
double *Anrms = NULL;

int imqr_step( matptr A1, mqrptr *pmqr, matptr *pA2,
               double tol_cos, double tol_drop,
               int lfil, int nosort );
int iqr( matptr mat, matptr R, double tol_drop, int lfil );

int comp_nodes( const void *buf1, const void *buf2 )
{
  const grphnd nd1 = (grphnd)buf1;
  const grphnd nd2 = (grphnd)buf2;
  /*if( nd1->deg < nd2->deg ) return -1;
  else if( nd1->deg > nd2->deg ) return 1;
  return 0;*/
  return nd1->deg - nd2->deg;
}

int imqr( matptr mat, mqrsptr mqrs, int nomulti, int nosort,
          double tol_level, int maxlvl, double tol_cos,
          double tol_cos_adaptive, double tol_drop,
          int lfil, double *wrk, int *iwrk, FILE *fp )
{
/*---------------------------------------------------------------------
 * Incomplete multi-level QR factorization
 *---------------------------------------------------------------------
 * Parameters
 *---------------------------------------------------------------------
 * on entry:
 * =========
 * mat      : matrix stored in SparMat format -- see heads.h for details
 * mqrs     : pointer to a SparMQRS struct -- see heads.h for details
 * nomulti  : 0: IMQR, 1: IQR
 * tol_level: tolerance for stopping MQR. Stop if
 *            reduced size < tol_level * previous size
 * tol_cos  :   orthogonal tolerance: if (v1,v2) < tol_cos*||v1||*||v2||
 *            vectors v1 and v2 are considered as orthogoanl
 * tol_drop :  dropping tolerance
 * lfil     : max # of elements allowed in each column
 * fp       : file pointer for error log (might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr != 0   --> error
 * mqrs
 *-----------------------------------------------------------------------
 * work array: wrk at least m + n, iwrk at least 2m + 2n, iwrk[i] = 0 
 *----------------------------------------------------------------------*/
  printf("begin imqr\n");
  int i, lvl = 0, m = mat->m, n = mat->n, min_n = (int)(0.01*n), *pos;
  matptr A1 = mat, A2 = NULL;
  clock_t clk;

  vecnval = wrk;
  vecmval = vecnval + n;
  vecnpat = iwrk;
  vecnidx = vecnpat + n;
  vecmpat = vecnidx + n;
  vecmidx = vecmpat + m;
  nodes = (grphnd)Malloc( n*sizeof(GraphNode), "imqr" );
  Anrms = (double *)Malloc( n*sizeof(double), "imqr" );

  setupMQRS( mqrs );
  clk = clock();
  if( nomulti == 0 ) {
    /* implement multilevel IQR */
    int done = 0;
    mqrptr mqr;
    while( !done ) {
      imqr_step( A1, &mqr, &A2, tol_cos, tol_drop, lfil, nosort );
      printf("END imqr_step OUT\n");
      mqrs->amqrs[lvl++] = mqr;
      if( (double)mqr->s <= tol_level * (double)mqr->n ) done = 1;
      if( A2->n < min_n ) done = 1;
      if( lvl >= maxlvl ) done = 1;
      if( lvl > 1 ) cleanMat( A1 );
      A1 = A2;
      tol_cos += tol_cos_adaptive; /* adaptive angle tolerance */
    }
  }
  printf("TEST--------0\n");
  clk = clock() - clk;
  mqrs->tm_lvl = clk/(double)CLOCKS_PER_SEC;

  clk = clock();
  printf("TEST--------2\n");
  iqr( A1, mqrs->R, tol_drop, lfil );
  clk = clock() - clk;
  mqrs->tm_iqr = clk/(double)CLOCKS_PER_SEC;
  if( lvl > 1 ) cleanMat( A1 );

  mqrs->n = mat->n;
  mqrs->lvl = lvl;
  pos = mqrs->pos;
  pos[0] = 0;
  if( lvl > 0 ) {
    for( i = 1; i <= lvl; i++ ) {
      pos[i] = pos[i-1] + mqrs->amqrs[i-1]->s;
    }
  }

  free( nodes );
  free( Anrms );


  return 0;
}

int imqr_step( matptr A1, mqrptr *pmqr, matptr *pA2,
               double tol_cos, double tol_drop, int lfil, int nosort )
{
  printf("begin imqr_step\n");
  int m = A1->m, n = A1->n, i, j, j1, j2, s, row, col, id, flag;
  int *perm = NULL, *iperm = NULL;
  SparVec Ei, a2i;
  matptr A2 = NULL, graph = NULL, E = NULL;
  matptr max_ortho_columns( matptr, int *, int *, double, int, int * );
  mqrptr mqr;
  mqr = (mqrptr)Malloc( sizeof(SparMQR), "imqr_step" );
  mqr->E = (SparMat *)Malloc( sizeof(SparMat), "imqr_step" );
  E = mqr->E;
  Vector vec, vecE;
  double alpha, diag, nrm2;
  clock_t clk;
  setupMQR( mqr, n );
  perm = mqr->perm;
  iperm = mqr->iperm;

  /* A1 = ( A11 | A12 ), where columns of A11 are orthogonal to each other */
  clk = clock();
  /* calculate 2-norm for each column of A in order to calculate
   * the cosine of the angle between two columns later */
  for( i = 0; i < n; i++ ) {
    loadVectorC( &vec, A1, i );
    Anrms[i] = sparvec_nrm2( &vec );
  }
  graph = max_ortho_columns( A1, perm, iperm, tol_cos, nosort, &s );
  clk = clock() - clk;
  /*fprintf( stderr, "max_othor time: %g\n", clk/(double)CLOCKS_PER_SEC );*/

  /* D_1 = diag( 1/||A11(1)||, 1/||A11(2)||, ..., 1/||A11(s)|| ) */
  mqr->s = s;
  mqr->D_1 = (double *)Malloc( s*sizeof(double), "imqr_step" );
  for( i = 0; i < s; i++ ) {
    diag = Anrms[perm[i]];
    if( diag < MIN_NUM ) mqr->D_1[i] = 0.0;
    else mqr->D_1[i] = 1.0 / diag;
  }
  Ei.idx = vecnidx;
  Ei.pat = vecnpat;
  Ei.val = vecnval;

  setupMat( E, s, n-s, MAT_COLBASE_INDEX|MAT_COLBASE_VALUE,
            A1->bufsz/2, A1->incsz/2 );
  E->cvec[0] = E->nnz;
  /* calculate E = Q1^T * A12, where Q1 = A11 * inv(D) */
  /* the graph matrix A^T * A can be reused */
  clk = clock();
  for( i = 0; i < n-s; i++ ) {
    loadVectorR( &vec, graph, perm[s+i] );
    Ei.nnz = 0;
    for( j = 0; j < vec.nnz; j++ ) {
      id = iperm[vec.idx[j]];
      if( id < s ) {
        Ei.idx[Ei.nnz] = id;
        Ei.val[Ei.nnz] = vec.val[j] * mqr->D_1[id];
        Ei.nnz++;
      }
    }

    /* copy Ei to the i-th column of E */
    checkMatBuffer( E, Ei.nnz );
    for( j = 0; j < Ei.nnz; j++ ) {
      E->cidx[E->nnz] = Ei.idx[j];
      E->cval[E->nnz] = Ei.val[j];
      E->nnz++;
    }
    E->cvec[i+1] = E->nnz;
  }
  cleanMat( graph );
  clk = clock() - clk;
  /*fprintf( stderr, "Calculating E: %g\n", clk/(double)CLOCKS_PER_SEC );*/
  a2i.val = vecmval;
  a2i.pat = vecmpat;
  a2i.idx = vecmidx;
  A2 = (matptr)Malloc( sizeof(SparMat), "imqr_step" );
  flag = MAT_ROWBASE_INDEX|MAT_ROWBASE_VALUE|
         MAT_COLBASE_INDEX|MAT_COLBASE_VALUE;
  setupMat( A2, m, n-s, flag, A1->bufsz, A1->incsz );
  A2->cvec[0] = A2->nnz;
  /* Calculate A2_new = A2 - Q1 * E */
  clk = clock();
  for( i = 0; i < n-s; i++ ) {
    fsparvec_loadC( &a2i, A1, perm[s+i] );
    loadVectorC( &vecE, E, i );
    for( j = 0; j < vecE.nnz; j++ ) {
      col = vecE.idx[j];
      /* A2i = A2i - <A2i,qj>*qj */
      alpha = - vecE.val[j] * mqr->D_1[col];
      loadVectorC( &vec, A1, perm[col] );
      fsparvec_add( &a2i, &vec, alpha );
    }

    /* copy a2i to the i-th column of A2 */
    checkMatBuffer( A2, a2i.nnz );
    for( j = 0; j < a2i.nnz; j++ ) {
      row = a2i.idx[j];
      A2->cidx[A2->nnz] = row;
      A2->cval[A2->nnz] = a2i.val[row];
      A2->nnz++;
      a2i.pat[row] = 0; /* reset pat to 0 array */
    }
    A2->cvec[i+1] = A2->nnz;
  }
  transMat( A2 ); /* store A2 row-wise */
  clk = clock() - clk;
  /*fprintf( stderr, "Calculating A2: %g\n", clk/(double)CLOCKS_PER_SEC );*/

  *pmqr = mqr;
  *pA2 = A2;
  printf("END imqr_step\n");
  return 0;
}

matptr max_ortho_columns( matptr A, int *perm, int *iperm,
                          double tol_cos, int nosort, int *ps )
{
  int n = A->n, i, j, k, row, col, id, pm, ipm, s, nnz, j1, j2;
  double cs, cosnrm, aTval;
  Vector rvec;
  SparVec ci;
  matptr C;

  /* initial ci, used as the i-th row of C */
  ci.val = vecnval;
  ci.pat = vecnpat;
  ci.idx = vecnidx;

  C = (matptr)Malloc( sizeof(SparMat), "max_ortho_columns" );
  setupMat( C, n, n, MAT_ROWBASE_INDEX|MAT_ROWBASE_VALUE,
            A->bufsz*2, A->incsz );
  C->rvec[0] = C->nnz;

  /* nodes will be used to store the degree of each vertex, so that it */
  /* can be sorted by qsort                                            */
  for( i = 0; i < n; i++ ) {
    nodes[i].deg = 0;
    nodes[i].idx = i;
  }

  for( i = 0; i < n; i++ ) {
    /* calculate the i-th row of C: ci = ai^T * a(:,i+1:n) */
    ci.nnz = 0;
    j1 = A->cvec[i];
    j2 = A->cvec[i+1];
    for( j = j1; j < j2; j++ ) {
      row = A->cidx[j];
      aTval = A->cval[j];
      loadVectorR( &rvec, A, row );
      for( k = rvec.nnz-1; k >= 0; k-- ) {
        /* row indices are in increasing order, stop when idx <= i */
        col = rvec.idx[k];
        if( col <= i ) break;
        if( ci.pat[col] == 0 ) {
          /* fill-in */
          ci.idx[ci.nnz] = col;
          ci.pat[col] = 1;
          ci.val[col] = aTval * rvec.val[k];
          ci.nnz++;
        } else {
          ci.val[col] += aTval * rvec.val[k];
        }
      }
    }

    /* setup the upper part of the graph of normalized A^T * A */
    /* copy ci to the i-th row of C */
    nnz = nodes[i].deg + ci.nnz;
    checkMatBuffer( C, nnz );
    C->nnz += nodes[i].deg;
    if( ci.nnz > 0 ) {
      /* if (ai,aj) < tol_cos * |ai| * |aj|, assume ai is orthogonal to aj */
      cosnrm = tol_cos * Anrms[i];
      for( j = 0; j < ci.nnz; j++ ) {
        col = ci.idx[j];
        cs = ci.val[col];
        if( fabs(cs) > cosnrm * Anrms[col] ) {
          nodes[col].deg++;
          C->ridx[C->nnz] = col;
          C->rval[C->nnz] = cs;
          C->nnz++;
        }
        ci.pat[col] = 0; /* reset pattern to 0 */
      }
    }
    C->rvec[i+1] = C->nnz;
  }

  for( i = 0; i < n; i++ ) {
    ci.idx[i] = C->rvec[i] + nodes[i].deg; /* pointing to boundary index */
    nodes[i].deg = C->rvec[i]; /* pointing to the next available position */
  }
  /* convert to full graph from the upper part */
  for( i = 0; i < n; i++ ) {
    j1 = ci.idx[i];
    j2 = C->rvec[i+1];
    for( j = j1; j < j2; j++ ) {
      row = C->ridx[j];
      C->ridx[nodes[row].deg] = i;
      C->rval[nodes[row].deg] = C->rval[j];
      nodes[row].deg++;
    }
  }
  for( i = 0; i < n; i++ ) {
    nodes[i].deg = C->rvec[i+1] - C->rvec[i];
  }

  if( !nosort ) {
    qsort( nodes, n, sizeof(GraphNode), comp_nodes );
  }
  s = 0;
  for( i = 0; i < n; i++ ) {
    perm[i] = i;
    iperm[i] = i;
    ci.idx[i] = 0;
  }
  for( i = 0; i < n; i++ ) {
    id = nodes[i].idx;
    if( ci.idx[id] == 1 ) continue;

    pm = perm[s];
    ipm = iperm[id];
    perm[s] = id;
    perm[ipm] = pm;
    iperm[id] = s;
    iperm[pm] = ipm;
    s++;

    j1 = C->rvec[id];
    j2 = C->rvec[id+1];
    for( j = j1; j < j2; j++ ) {
      ci.idx[C->ridx[j]] = 1;
    }
  }
  
  *ps = s;
  return C;
}

typedef struct tagElement
{
    double  val;
    int     col;
    int     next;
} Element;

int iqr( matptr mat, matptr R, double tol_drop, int lfil )
{
  /*==================================================================*/
  /* IQR using Saad's ILQ method                                      */
  /*==================================================================*/
  /* Input:                                                           */
  /*                                                                  */
  /*==================================================================*/
  /* work array: wk at least m + n, iwk at least 2m + n               */
  /* iwk initially is 0 array                                         */
  /*==================================================================*/
  /*                                                                  */
  /*==================================================================*/
  int m = mat->m, n = mat->n, i, j, k, nnz;
  int Apos, Qnext, Ridx, QTpos;
  double Aval, Rval, invRval, value, tnorm;
  int szBase, szQlist, szQbuf, szQinc, Qnnz, index;
  const int *idx, *Qidx;
  const double *val, *Qval;


  matptr Q = (matptr)Malloc( sizeof(SparMat), "iqr" );
  SparVec qi, ri;
  Element *Qlist;
  Vector cvec, vec;


  /* initialize Q and R (column wise) */
  setupMat( Q, m, n, MAT_COLBASE_INDEX|MAT_COLBASE_VALUE,
            mat->bufsz, mat->incsz );
  setupMat( R, n, n, MAT_COLBASE_INDEX|MAT_COLBASE_VALUE,
            n+mat->bufsz, mat->incsz );
  Q->cvec[0] = Q->nnz;
  R->nnz = n;
  R->cvec[0] = R->nnz; /* the first n elements are diagonals */

  /* initialize buffers for qi and ri (the i-th column of Q and R resp.) */
  qi.val = vecmval;
  qi.pat = vecmpat;
  qi.idx = vecmidx;
  ri.val = vecnval;
  ri.pat = vecnpat;
  ri.idx = vecnidx;

  /* initialize the linked list of Q (row wise) using a linear array */
  /* Qlist (cf. Saad's ILQ). Qlist[i] points to the i-th row of Q    */
  /* szQbuf: the current buffer size                                 */
  /* szQinc: the size to increase when the buffer is full            */
  /* QTpos : the next available position in Qlist                    */
  szBase  =   nnzMat( mat );
  szQlist =   m;
  szQbuf  =   m + mat->bufsz;
  szQinc  =   mat->incsz;
  Qlist   =   (Element *)malloc( szQbuf * sizeof(Element) );
  for( i = 0; i < m; i++ )
  {
      Qlist[i].next = -1;
  }
  QTpos = m;

  for( i = 0; i < n; i++ )
  {
    loadVectorC( &cvec, mat, i );
    idx = cvec.idx;
    val = cvec.val;
    nnz = cvec.nnz;
        
    /* reset the i-th column of Q and R */
    qi.dim  =   m;
    ri.dim  =   n;
    qi.nnz  =   0;
    ri.nnz  =   0;

    /* load the i-th column of A and use Saad's ILQ method to calculate */
    /* entries in the i-th column of R                                  */
    tnorm = 0.0;
    for( j = 0; j < nnz; j++ )
    {
      Apos = idx[j];
      Aval = val[j];

      tnorm += fabs( Aval );

      qi.idx[qi.nnz] = Apos;
      qi.val[Apos] = Aval;
      qi.pat[Apos] = 1;
      qi.nnz++;
      Qnext = Qlist[Apos].next;
      while( Qnext != -1 )
      {
        index = Qlist[Qnext].col;
        value = Qlist[Qnext].val;

        if( ri.pat[index] == 0 )
        {
          /* fill-in */
          ri.idx[ri.nnz]  =   index;
          ri.pat[index]   =   1;
          ri.val[index]   =   value * Aval;
          ri.nnz++;
        }
        else
        {
          /* not a fill-in */
          ri.val[index]  +=   value * Aval;
        }
        Qnext = Qlist[Qnext].next;
      }
    }
    tnorm *= tol_drop;

    for( j = 0; j < ri.nnz; j++ )
    {
      ri.pat[ ri.idx[j] ] = 0;    /* reset pattern to 0 */
    }
    /* keep only the lfil largest elements (in absolute value) in R(:,i) */
    if( ri.nnz > lfil ) {
      qsplit2( ri.val, &n, ri.idx, &ri.nnz, &lfil );
      ri.nnz = lfil;
    }

    /* calculate Ai - Sum <Ai,Qj>Qj, i.e., Ai - Sum(rj*Qj) */
    /* store Ri = <Ai,Qj> */
    checkMatBuffer( R, ri.nnz );
    for( j = 0; j < ri.nnz; j++ )
    {
      Ridx = ri.idx[j];
      Rval = ri.val[Ridx];

      loadVectorC( &vec, Q, Ridx );
      fsparvec_add( &qi, &vec, -Rval );

      if( fabs(Rval) >= tnorm ) { /* drop very small elements */

        R->cidx[R->nnz] = ri.idx[j];
        R->cval[R->nnz] = Rval;
        R->nnz++;
      }
    }
    R->cvec[i+1] = R->nnz;

    for( j = 0; j < qi.nnz; j++ )
    {
      qi.pat[ qi.idx[j] ] = 0; /* reset pattern to 0 */
    }
    /* keep only the lfil largest elements (in absolute value) in Q(:,i) */
    if( qi.nnz > lfil ) {
      qsplit2( qi.val, &m, qi.idx, &qi.nnz, &lfil );
      qi.nnz = lfil;
    }

    Rval = fsparvec_nrm2( &qi );
    if( Rval == 0.0 ) {
      fprintf( stderr, "Zero column encountered.\n" );
      Rval = MIN_NUM;
    }
    invRval = 1.0 / Rval;
    R->cval[i] = invRval; /* the first n elements store inv of R diag */

    /* update Q and Qlist */
    if( szQlist + qi.nnz > szQbuf )
    {
        /* if the current buffer is not large enough, enlarge it */
        szQbuf += szQinc;
        while( szQlist + qi.nnz > szQbuf )
        {
            szQbuf += szQinc;
        }
        Qlist = (Element *)realloc( Qlist, szQbuf * sizeof(Element) );
    }
    checkMatBuffer( Q, qi.nnz );
    for( j = 0; j < qi.nnz; j++ )
    {
      index = qi.idx[j];
      value = qi.val[index] * invRval;

      /* ignore very small entries */
      /* if( fabs(value) < tnorm ) continue; */

      Q->cidx[Q->nnz] = index;
      Q->cval[Q->nnz] = value;
      Q->nnz++;

      /* add to Qlist's row Qidx */
      Qlist[QTpos].col = i;
      Qlist[QTpos].val = value;
      Qlist[QTpos].next = Qlist[index].next;
      Qlist[index].next = QTpos;
      szQlist++;
      QTpos++;
    }
    Q->cvec[i+1] = Q->nnz;
  }

  if( Qlist ) free( Qlist );
  cleanMat( Q );
  return 0;
}

/*-----------------------------------------------------------------------*/
int cond_est( matptr mat, mqrsptr mqrs, double *y, double *x, FILE *fp )
{
  int n = mqrs->n, i;
  double norm;
  for(i = 0; i < n; i++ ) {
    y[i] = 1.0;
  }
  printf("begin cond_est------1\n");
  mqrsol( mqrs, y, x, 1, y );
  printf("begin cond_est------2\n");
  mqrsol( mqrs, x, x, 0, y );
  norm = vec_nrminf( n, x );
  printf("begin cond_est------3\n");
  fprintf(fp, "  MQR inf-norm lower bound : %16.2f\n", norm );
  if( norm > 1e30 ) {
    return -1;
  }
  return 0;
}
/*-----------------------------------------------------------------------*/
