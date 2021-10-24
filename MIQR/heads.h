/*---------------------------------------------------------------------------*
 * Sparse Matrix/Vector Data Structures                                      *
 *---------------------------------------------------------------------------*
 * Na Li, Jul 12, 2004                                                       *
 *                                                                           *
 * Report bugs / send comments to: nli@cs.umn.edu                            *
 *---------------------------------------------------------------------------*/

/* sparse vector that allows fill-ins to be inserted */
#define MAT_ROWBASE_INDEX    (0x00000001)
#define MAT_ROWBASE_VALUE    (0x00000002)
#define MAT_COLBASE_INDEX    (0x00000004)
#define MAT_COLBASE_VALUE    (0x00000008)
#define CSR_INDEX(flag)      ((flag) & MAT_ROWBASE_INDEX)
#define CSR_VALUE(flag)      ((flag) & MAT_ROWBASE_VALUE)
#define CSC_INDEX(flag)      ((flag) & MAT_COLBASE_INDEX)
#define CSC_VALUE(flag)      ((flag) & MAT_COLBASE_VALUE)
typedef struct _SparVector {
  int dim;      /* the size of the vector */
  int nnz;      /* the # of nonzeros */
  double *val;  /* full vector array */
  int *pat;     /* pat[i] = 1 if i-th element is nonzero, otherwise 0 */
  int *idx;     /* idx[1..nnz] are the indices of nonzeros */
} SparVec, * fvecptr;

/* sparse vector */
typedef struct _Vector {
    int     nnz;       /* the # of nonzeros */
    const int *idx;    /* the row/column indices */
    const double *val; /* the values */
} Vector, * vecptr;

typedef struct _SparMat {
    int m, n;          /* size of the matrix (m x n) */

    int *rvec;         /* matrix stroed in row-based. the i-th row indices
                        * are from rvec[i] to rvec[i+1]-1 */
    int *ridx;          /* index array row wise */
    double *rval;       /* value array row wise */

    int *cvec;         /* matrix stored in column-based. the i-th column
                        * indices are from cvec[i] to cvec[i+1]-1 */
    int *cidx;          /* index array column wise */
    double *cval;       /* value array column wise */

    int nnz;           /* current # of entries */
    int bufsz;         /* current buffer size */
    int incsz;         /* the size to increase when needd */
    int flag;
/*=======================================================================*/
/*=========================================================================
 * if (flag & MAT_ROWBASE_INDEX) ==> rvec[row].idx is available
 * if (flag & MAT_ROWBASE_VALUE) ==> rvec[row].val is available
 * if (flag & MAT_COLBASE_INDEX) ==> cvec[col].idx is available
 * if (flag & MAT_COLBASE_VALUE) ==> cvec[col].val is available
**=======================================================================*/
} SparMat, * matptr;

typedef struct _SparMQR {
  int n;
  int s;
  int *perm;
  int *iperm;
  double *D_1;  /* inv(D) */
  matptr E;
} SparMQR, * mqrptr;

typedef struct _SparMQRS {
    int n;
    int lvl;
#define MAX_MQR_LEVEL 100
    mqrptr amqrs[MAX_MQR_LEVEL];
    int pos[MAX_MQR_LEVEL+1];
    matptr R;
    double tm_lvl;
    double tm_iqr;
} SparMQRS, * mqrsptr;

typedef struct _GraphNode {
  int deg;     /* degree */
  int idx;     /* index  */
} GraphNode, * grphnd;

