from libc.string cimport memset, strlen
from libc.stdlib cimport atoi, free, calloc
from libc.stdio cimport stdout, fopen, sprintf, fprintf, fclose, printf, fgets
cimport cython
from scipy import sparse
from time import time
import os
cdef extern from "stdio.h":
    ctypedef struct _iobuf:
        char *_ptr;
        int   _cnt;
        char *_base;
        int   _flag;
        int   _file;
        int   _charbuf;
        int   _bufsiz;
        char *_tmpfname;
    ctypedef _iobuf FILE;
cdef extern from "heads.h":
    cdef enum:
        MAT_ROWBASE_INDEX = (0x00000001)
        MAT_ROWBASE_VALUE = (0x00000002)
        MAT_COLBASE_INDEX = (0x00000004)
        MAT_COLBASE_VALUE = (0x00000008)
        MAX_MQR_LEVEL = 100
    ctypedef struct _SparMat:
        int m, n;
        int *rvec;
        int *ridx;
        double *rval;
        int *cvec;
        int *cidx;
        double *cval;
        int nnz;
        int bufsz;
        int incsz;
        int flag;

    ctypedef _SparMat SparMat
    ctypedef _SparMat * matptr
    ctypedef struct _SparMQR:
        int n;
        int s;
        int *perm;
        int *iperm;
        double *D_1;
        matptr E;
    ctypedef _SparMQR * mqrptr
    ctypedef _SparMQR SparMQR

    ctypedef struct _SparMQRS:
        int n;
        int lvl;
        mqrptr amqrs[MAX_MQR_LEVEL];
        int pos[MAX_MQR_LEVEL+1];
        matptr R;
        double tm_lvl;
        double tm_iqr;
    ctypedef _SparMQRS * mqrsptr
    ctypedef _SparMQRS SparMQRS
cdef extern from "defs.h":
    cdef enum:
        MAX_PATH = 256
        MAX_LINE = 256

    void *Malloc( int, char * );
    ctypedef struct _io_t:
        FILE *fout;                 #/* output file handle              */
        char outfile[MAX_PATH];     #/* output filename                 */
        char Fname[MAX_PATH];       #/* matrix filename                 */
        char HBnameF[MAX_PATH];     #/* HB name                         */
        char title[73];             #/* Matrix title                    */
        char key[9];                #/* Matrix key                      */
        char type[4];               #/* HB type                         */
        int  nrhs;                  #/* the # of rhs                    */
        double *rhs;                #/* the nrhs right hand side        */
        double *guess;              #/* the nrhs initial guess          */
        double *sol;                #* the nrhs exact solution         */
        int mdim, ndim;             #* matrix size m x n               */
        int nnz;                    #/* number of nonzero               */
        int nparam;                 #/* number of tests for each matrix */
        double tol_level;           #/* tolerance (between 0 and 1) to stop MQR */
        int maxlvl;                 #/* maximal levels allowed */
        double tol_cos;             #/* tolerance for orthogonal */
        double tol_cos_adaptive;    #/* increment for each level */
        double tol_delta;           #/* increment for orthogonal tolerance */
        double tol_drop;            #/* tolerance for dropping elements */
        double tol_lfil;            #/* only tol_lfil*(nnz/n) fill-ins allowed */
        int lfil;                   #/* = tol_lfil * (nnz/n) */
        int maxits;                 #/* maximum number of CG iters  */
        double tol_its;             #/* tolerance for stopping CG       */
        double tm_p;                #/* time for preconditioner (s)     */
        double tm_i;                #/* time for iteration (s)          */
        double tm_total;            #/* totoal time tm_p + Sum(tm_i)_nrhs */
        int fillin;                 #/* memory used during precondition */
        int its;                    #/* number of iterations            */
        double enorm;               #/* error norm:          || x- x0|| */
        double rnorm;               #/* final residual norm: ||Ax-Ax0|| */
    ctypedef _io_t io_t;
    int read_inputs( char *in_file, io_t *pio );
    int get_matrix_info( FILE *fmat, io_t *pio );
    int readhbc_c( matptr mat, io_t *pio, int manu_rhs );
    int output_header( io_t *pio );
    int imqr( matptr mat, mqrsptr mqrs, int nomulti, int nosort,
          double tol_level, int maxlvl, double tol_cos,
          double tol_cos_adaptive, double tol_drop,
          int lfil, double *wrk, int *iwrk, FILE *fp );
    int output_prec( io_t *pio, mqrsptr mqrs, int output_levels, int nparam );
    int output_result( io_t *pio, int iparam, int nrhs );


cdef extern from "matfun.h":
    void * Malloc( int, char * );
    int nnzMQRS( mqrsptr mqrs );
    int cond_est( matptr mat, mqrsptr mqrs, double *y, double *x, FILE *fp );
    int pcgnr( matptr mat, mqrsptr mqrs, double *rhs, double *sol,
           double tol, int *itmax, int flag );
    double vec_nrm2_diff( int, double *, double * );
    int matxvec( matptr mat, double *x, double *y );
    int matTxvec( matptr mat, double *x, double *y );
    int vec_sub( int n, double *v1, double *v2, double *x );
    double vec_nrm2( int, double * );
    int vec_output( int n, double *v, char *name, FILE *fp );
    int cleanMQRS( mqrsptr mqrs );
    int cleanMat( matptr mat );


# def  CSR_INDEX(flag):
#     return ((flag) & MAT_ROWBASE_INDEX)
# def CSR_VALUE(flag):
#     return ((flag) & MAT_ROWBASE_VALUE)
# def CSC_INDEX(flag):
#     return ((flag) & MAT_COLBASE_INDEX)
# def CSC_VALUE(flag):
#     return ((flag) & MAT_COLBASE_VALUE)

def help():
    print( "======================================================\n" )
    print( "Usage: parameters [np] [nm] [ns] [nr] [l] [k] [o] [h]\n" )
    print( "======================================================\n" )
    print( "Options:\n" )
    print( "np: 1/true to disable preconditioner when solving the system\n" )
    print( "nm: 1/true to disable multi-level preconditioner (IQR instead)\n" )
    print( "ns: 1/true to do NOT sort the graph A'*A according the degrees\n" )
    print( "nr: 1/true to do NOT read rhs from file. use artificial rhs\n" )
    print( "l:  1/true to output the reduced size for each level\n" )
    print( "k:  1/true to skip iterations\n" )
    print( "o:  1/true to output the solutions\n" )
    print( "h:  1/true to display this help\n" )
    print( "A:  matrix A x = b\n" )
    print( "b:  matrix A x = b\n" )
    print( "======================================================\n" )

def changeMat(matpath):
    if os.path.exists(matpath):
        filenames = os.listdir(matpath)
        Nfile = str(len(filenames)) + "\n"
        filePath = [os.path.join(matpath, x) for x in filenames]
        with open("matfile", "w") as f:
            f.write(Nfile)
            for i in range(len(filenames)):
                f.write( "'" + filePath[i] + "' '" + filenames[i] + "' \n")
        return 0
    else:
        return -1
class Inputs(object):
    def __init__(self, nparam = 11, tol_level = 0.3, maxlvl = 4, tol_cos = 0.00, tol_cos_adaptive = 0.01,
                 tol_delta = 0.02, tol_drop = 0.01, lfil = 6.0, maxits = 2000, tol_its = 1.0e-8):

        self.nparam = nparam
        self.tol_level = tol_level
        self.maxlvl = maxlvl
        self.tol_cos = tol_cos
        self.tol_cos_adaptive = tol_cos_adaptive
        self.tol_delta = tol_delta
        self.tol_drop = tol_drop
        self.tol_lfil = lfil
        self.maxits = maxits
        self.tol_its = tol_its

    def usage(self):
        print("=================================================================")
        print("                      parameters of Inputs                       ")
        print("=================================================================")
        print("nparam    = number of tests for each matrix")
        print("tol_level = tolerance for stopping MQR, typical 0.20")
        print("maxlvl    = maximal levels allowed, typical 4")
        print("tol_cos   = tolerance for orthogonal, typical 0.03")
        print("tol_cos_adaptive = increment for each additional level")
        print("tol_delta = increment for the orthogonal threshold after each test")
        print("tol_drop  = tolerance for dropping small terms")
        print("lfil      = up to lfil*(nnz/n) of fills allowed for each column")
        print("maxits    = maxits in outer fgmres.")
        print("tol_its   = tolerance for stopping iteration")
        print("=================================================================")


def preconIMQR(*args, **kwargs):
    # for key in kwargs:
    #     print(key)
    cdef int ierr = 0;


    cdef int skip_its = 0, output_levels = 0, output_sol = 0;
    cdef int noprec = 0, nomulti = 0, nosort = 0, manu_rhs = 0;
    cdef double tol_cos;
    cdef int nrhs;
    cdef double *sol = NULL, *guess = NULL, *rhs = NULL;
    cdef double *x = NULL, *y = NULL;
    cdef double *wk = NULL;
    cdef int *iwk = NULL;
    cdef int nmat, numat, iparam, i;
    cdef matptr mat = NULL;

    cdef mqrsptr mqrs = NULL;
    cdef FILE *fmat = NULL;
    cdef io_t io;
    cdef FILE *flog = stdout;
    print("beginning setting choice......................")
    for key in kwargs:
        if ("A" not in kwargs.keys()) or ("b" not in kwargs.keys()):
            print("not A or b data..........")
            return 0
        else:
            A = kwargs["A"]
            b = kwargs["b"]
        if (key != "A") and (key != "b") and (kwargs[key] != 0):
            if key == "np":
                noprec = 1
            if key == "nm":
                nomulti = 1
            if key == "ns":
                nosort = 1
            if key == "nr":
                manu_rhs = 1
            if key == "k":
                skip_its = 1
                break
            if key == "l":
                output_levels = 1
                break
            if key == "o":
                output_sol = 1
                break
            if key == "h":
                help()
                return 0

    memset( &io, 0, sizeof(io))
    print("get inputsfile param......................")
    inputs = Inputs()
    io.tol_its = inputs.tol_its
    io.maxits = inputs.maxits
    io.tol_lfil = inputs.tol_lfil
    io.tol_delta = inputs.tol_delta
    io.tol_cos_adaptive = inputs.tol_cos_adaptive
    io.nparam = inputs.nparam
    io.tol_level = inputs.tol_level
    io.maxlvl = inputs.maxlvl
    io.tol_drop = inputs.tol_drop
    io.tol_cos = inputs.tol_cos
    io.nrhs = 1
    # io.sol = <double *>Malloc( io.nrhs * ncol * sizeof(double), "preconIMQR" )


    # if( read_inputs( INPUT, &io ) != 0 ):
    #     fprintf( flog, "Invalid inputs file...\n" )
    #     return -1
    # print("beginning read file......................1")
    # fmat = fopen( MAT, "r" )
    # if( NULL == fmat  ):
    #     fprintf( flog, "Can't open matfile...\n" )
    #     return -1
    # print("beginning read file......................2")
    # memset( line, 0, MAX_LINE )
    # fgets( line, MAX_LINE, fmat )
    # print("beginning read file......................3")
    #
    # with open("inputs", "r") as f:
    #     line = f.readlines()
    # numat = atoi( line )
    # if( numat <= 0 ):
    #     fprintf(flog, "Invalid count of matrices...\n" )
    #     return -1

# for nmat in range(numat):
#     if( get_matrix_info( fmat, &io ) != 0 ):
#       print("Invalid format in matfile...\n" )
    print("\n========================================\n" )
    print("            Solving MATRIX\n")
    print("========================================\n" )

    mat = <matptr>Malloc( sizeof(SparMat), "preconIMQR" )
    AmatrixC = sparse.csc_matrix(A)
    AmatrixR = sparse.csr_matrix(A)

    mat.nnz = AmatrixC.nnz
    mat.m = AmatrixC.shape[0]
    mat.n = AmatrixC.shape[1]
    mat.bufsz = AmatrixC.nnz
    mat.incsz = AmatrixC.nnz

    mat.cvec = <int *>Malloc( (len(AmatrixC.indptr) ) * sizeof(int), "preconIMQR" )
    mat.cidx = <int *>Malloc( (len(AmatrixC.indices) ) * sizeof(int), "preconIMQR" )
    mat.cval = <double *>Malloc( (len(AmatrixC.data) ) * sizeof(double), "preconIMQR" )

    mat.rvec = <int *>Malloc( (len(AmatrixR.indptr) ) * sizeof(int), "preconIMQR" )
    mat.ridx = <int *>Malloc( (len(AmatrixR.indices) ) * sizeof(int), "preconIMQR" )
    mat.rval = <double *>Malloc( (len(AmatrixR.data) ) * sizeof(double), "preconIMQR" )

    for i in range(len(AmatrixC.indptr)):
        mat.cvec[i] = AmatrixC.indptr[i]

    for i in range(len(AmatrixC.indices)):
        mat.cidx[i] = AmatrixC.indices[i]

    for i in range(len(AmatrixC.data)):
        mat.cval[i] = AmatrixC.data[i]

    for i in range(len(AmatrixR.indptr)):
        mat.rvec[i] = AmatrixR.indptr[i]

    for i in range(len(AmatrixR.indices)):
        mat.ridx[i] = AmatrixR.indices[i]

    for i in range(len(AmatrixR.data)):
        mat.rval[i] = AmatrixR.data[i]

    io.mdim = mat.m
    io.ndim = mat.n
    io.nnz = mat.nnz
    io.nrhs = 1
    io.rhs = <double *>Malloc( io.nrhs * io.mdim * sizeof(double), "preconIMQR" )
    io.guess = <double *>Malloc( io.nrhs * io.ndim*sizeof(double), "preconIMQR" )
    io.sol = <double *>Malloc( io.nrhs*io.ndim*sizeof(double), "preconIMQR" );
    for i in range(len(b)):
        io.rhs[i] = b[i]
    for i in range(io.ndim):
          io.guess[i] = 0.0

    io.lfil = <int>(io.tol_lfil * <double>io.nnz / <double>io.ndim )
    if( io.lfil < 1 ):
        io.lfil = 1
    # io.HBnameF = "youlike"
    # io.type = "csc"
    # output_header( &io )
    tol_cos = io.tol_cos


    x = <double *>Malloc( io.ndim * sizeof(double), "preconIMQR" )
    y = <double *>Malloc( max(io.ndim,io.mdim)*sizeof(double), "preconIMQR" )
    wk = <double *>Malloc( (io.mdim+io.ndim)*sizeof(double), "preconIMQR" )
    iwk = <int *>Malloc( 2*(io.mdim+io.ndim)*sizeof(int), "preconIMQR" )
    memset( iwk, 0, 2*(io.mdim+io.ndim)*sizeof(int) )

    for iparam in range(io.nparam):
        print("**** test # %d ****\n"%(iparam+1) )
        mqrs = <mqrsptr>Malloc( sizeof(mqrsptr), "preconIMQR" )#SparMQRS
        print("begin imqr-------------------------")
        clk = time()

        ierr = imqr( mat, mqrs, nomulti, nosort, io.tol_level, io.maxlvl,
                   tol_cos, io.tol_cos_adaptive, io.tol_drop, io.lfil,
                   wk, iwk, flog )
        clk = time() - clk
        io.tm_p = <double>clk
        io.tm_total = io.tm_p
        if( ierr != 0 ):
            print("*** mqr error, ierr != 0 ***\n" )
            return 0

        io.fillin = nnzMQRS( mqrs )
        # print("  MQR: level = %d, time = %8.1e (s), mem used = %d\n"%
        #       (mqrs.lvl, io.tm_p, io.fillin ))
        # output_prec( &io, mqrs, output_levels, iparam )
        print("begin cond_est")
        cierr = cond_est( mat, mqrs, y, x, flog )
        print("end cond_est")
        if( skip_its ):

            io.its = -1
            io.tm_i = -1
            io.enorm = -1
            io.rnorm = -1
        elif( cierr != 0 ):
            print("  Not attempting iterative solution.\n" )
            print( "Not attempting iterative solution.\n" )
            io.its = -1
            io.tm_i = -1
            io.enorm = -1
            io.rnorm = -1
        else:
            print("begin for")
            for nrhs in range(io.nrhs):
                sol = io.sol + nrhs * io.ndim
                guess = io.guess + nrhs * io.ndim
                rhs = io.rhs + nrhs * io.mdim
                for i in range(io.ndim):
                  x[i] = guess[i]
                io.its = io.maxits
                clk = time()
                print("begin pcgnr")
                ierr = pcgnr( mat, mqrs, rhs, x, io.tol_its, &io.its, noprec )
                clk = time() - clk
                io.tm_i = <double>clk
                io.tm_total += io.tm_i
                print("  PCGNR: time = %8.1e (s), "%(io.tm_i) )
                if( ierr == 0 ):
                  print("Converged in %d iterations.\n" % (io.its ))
                else:
                  print("Not Converge in %d iterations.\n" % (io.maxits) )
                # /* calculate error norm */
                if( sol != NULL ):
                  io.enorm = vec_nrm2_diff( io.ndim, x, sol )
                else:
                  io.enorm = -1.0

                # /* calculate residual norm */
                matxvec( mat, x, y )
                vec_sub( io.mdim, y, rhs, y )
                matTxvec( mat, y, wk )
                io.rnorm = vec_nrm2( io.ndim, wk )

                # output_result( &io, iparam, nrhs )
                tempX = []
                for i in range(io.ndim):
                    tempX.append(x[i])

                if( output_sol ):
                    print("x is ", tempX)
                  # vec_output( io.ndim, x, "x", io.fout )
        tol_cos += io.tol_delta

        cleanMQRS( mqrs )
    _return__dic = {}
    if "sol" not in _return__dic.keys():
        _return__dic["sol"] = []
        for i in range(io.ndim):
            _return__dic["sol"].append(x[i])
    # fclose(io.fout)
    cleanMat( mat )
    if( x ):
        free( x )
    if( y ):
        free( y )
    if( wk ):
        free( wk )
    if( iwk ):
        free( iwk )
    if( io.rhs ):
        free( io.rhs )
        io.rhs = NULL
    if( io.guess ):
        free( io.guess )
        io.guess = NULL
    if( io.rhs ):
        free( io.rhs )
        io.rhs = NULL
    # if( flog is not stdout ):
    #     fclose ( flog )
    # fclose(fmat)
    return _return__dic["sol"]
