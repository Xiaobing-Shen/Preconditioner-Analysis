/*---------------------------------------------------------------------------*
 * Definitions of MIQR                                                       *
 *---------------------------------------------------------------------------*
 * Na Li, Jul 12, 2004                                                       *
 *                                                                           *
 * Report bugs / send comments to: nli@cs.umn.edu                            *
 *---------------------------------------------------------------------------*/

#define MAX_PATH 256
#define MAX_LINE MAX_PATH

typedef struct _io_t {
    FILE *fout;                 /* output file handle              */
    char outfile[MAX_PATH];     /* output filename                 */
    char Fname[MAX_PATH];       /* matrix filename                 */
    char HBnameF[MAX_PATH];     /* HB name                         */
    char title[73];             /* Matrix title                    */
    char key[9];                /* Matrix key                      */
    char type[4];               /* HB type                         */
    int  nrhs;                  /* the # of rhs                    */
    double *rhs;                /* the nrhs right hand side        */
    double *guess;              /* the nrhs initial guess          */
    double *sol;                /* the nrhs exact solution         */
    int mdim, ndim;             /* matrix size m x n               */
    int nnz;                    /* number of nonzero               */

/* parameters from inputs -----------------------------------------*/
    int nparam;                 /* number of tests for each matrix */
    double tol_level;           /* tolerance (between 0 and 1) to stop MQR */
    int maxlvl;                 /* maximal levels allowed */
    double tol_cos;             /* tolerance for orthogonal */
    double tol_cos_adaptive;    /* increment for each level */
    double tol_delta;           /* increment for orthogonal tolerance */
    double tol_drop;            /* tolerance for dropping elements */
    double tol_lfil;            /* only tol_lfil*(nnz/n) fill-ins allowed */
    int lfil;                   /* = tol_lfil * (nnz/n) */
    int maxits;                 /* maximum number of CG iters  */
    double tol_its;             /* tolerance for stopping CG       */

/* result for output ----------------------------------------------*/
    double tm_p;                /* time for preconditioner (s)     */
    double tm_i;                /* time for iteration (s)          */
    double tm_total;            /* totoal time tm_p + Sum(tm_i)_nrhs */
    int fillin;                 /* memory used during precondition */
    int its;                    /* number of iterations            */
    double enorm;               /* error norm:          || x- x0|| */
    double rnorm;               /* final residual norm: ||Ax-Ax0|| */
}   io_t;

void *Malloc( int, char * );
void *Realloc( void *, int, char * );
int read_inputs( char *in_file, io_t *pio );
int get_matrix_info( FILE *fmat, io_t *pio );
int output_header( io_t *pio );
int output_prec( io_t *pio, mqrsptr mqrs, int output_levels, int nparam );
int output_result( io_t *pio, int iparam, int nrhs );
int readhbc_c( matptr mat, io_t *pio, int manu_rhs );

#define max(a,b)  (((a)>(b))?(a):(b))

