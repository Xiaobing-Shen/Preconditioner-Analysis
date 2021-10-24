void *Malloc( int, char * );
void *Realloc( void *, int, char * );
double vec_nrminf( int, double * );
double vec_nrm2( int, double * );
double vec_nrm2_diff( int, double *, double * );
double vec_dot( int, double *, double * );
int loadVectorR( vecptr, matptr, int );
int loadVectorC( vecptr, matptr, int );
double sparvec_nrm1( vecptr );
double sparvec_nrm2( vecptr );
double fsparvec_nrm1( fvecptr );
double fsparvec_nrm2( fvecptr );
double fsparvec_dot( fvecptr, vecptr );
int matxvec( matptr mat, double *x, double *y );
int matTxvec( matptr mat, double *x, double *y );
int vec_add( int n, double *v1, double *v2, double *x );
int vec_sub( int n, double *v1, double *v2, double *x );
int vec_xadd( int n, double *v1, double *v2, double alpha, double *x );
int vec_cpy( int n, double *v1, double *v2 );
int vec_perm( int n, double *v1, int *perm, double *v2 );
int vec_output( int n, double *v, char *name, FILE *fp );
int fsparvec_loadR(  fvecptr, matptr, int );
int fsparvec_loadC(  fvecptr, matptr, int );
int fsparvec_add( fvecptr v1, vecptr v2, double alpha );
int fsparvec_iadd( fvecptr v1, vecptr v2, double alpha, double tol );
int CSC2Mat(int m, int n, double *a, int *ja, int *ia, matptr mat, int flag);
int setupMQRS( mqrsptr mqrs );
int cleanMQRS( mqrsptr mqrs );
int nnzMQRS( mqrsptr mqrs );
int setupMat( matptr mat, int m, int n, int flag, int bufsz, int incsz );
int cleanMat( matptr mat );
int checkMatBuffer( matptr mat, int nadd );
int nnzMat( matptr mat );
int printMat( matptr mat, char *name, FILE *fp );
int transMat( matptr mat );
int setupMQR( mqrptr mqr, int n );
int cleanMQR( mqrptr mqr );
int cond_est( matptr mat, mqrsptr mqrs, double *y, double *x, FILE *fp );

int mqrsol( mqrsptr mqrs, double *y, double *x, int transp, double *wk );
int imqr( matptr mat, mqrsptr mqrs, int nomulti, int nosort,
          double tol_level, int maxlvl, double tol_cos,
          double tol_cos_adaptive, double tol_drop,
          int lfil, double *wrk, int *iwrk, FILE *fp );
int pcgnr( matptr mat, mqrsptr mqrs, double *rhs, double *sol,
           double tol, int *itmax, int flag );

#ifdef _SGI
#define readmtc readmtc_
#define qsplit qsplit_
#define qsplit2 qsplit2_
#else
#ifdef _LINUX
#define readmtc readmtc_
#define qsplit qsplit_
#define qsplit2 qsplit2_
#else
#ifdef _IBM
#define readmtc readmtc
#define qsplit qsplit
#define qsplit2 qsplit2
#else
#define readmtc readmtc_
#define qsplit qsplit_
#define qsplit2 qsplit2_
#endif
#endif
#endif

void readmtc( int *nmax, int *nzmax, int *job, char *fname,
               double *a, int *ja, int *ia, double *rhs, int *nrhs,
               char *guesol, int *nrow, int *ncol, int *nnz, char *title,
               char *key, char *type, int *ierr );


void qsplit2( double *a, int *m, int *ind, int *nnz, int *ncut );
void qsplit( double *a, int *ind, int *n, int *ncut );

