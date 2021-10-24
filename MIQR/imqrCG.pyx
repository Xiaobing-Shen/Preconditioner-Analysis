from libc.stdio cimport stdout
from libc.string cimport memset
from scipy import sparse
import numpy as np
cimport cython
# cdef extern from "heads.h":
#     ctypedef struct _GraphNode:
#         int deg;     #/* degree */
#         int idx;
#     ctypedef _GraphNode GraphNode
#     ctypedef _GraphNode * grphnd
#
# cdef enum :
#     MIN_NUM = (1e-14)
# cdef double *vecnval = NULL, *vecmval = NULL;
# cdef int *vecnpat = NULL, *vecnidx = NULL;
# cdef int *vecmpat = NULL, *vecmidx = NULL;
# cdef grphnd nodes = NULL;
# cdef double *Anrms = NULL;
#
# cdef extern
cdef extern from "heads.h":
    cdef enum:
        MAX_MQR_LEVEL = 100
    ctypedef struct _SparMat:
        int m, n; # /* size of the matrix (m x n) */
        int *rvec; #/* matrix stroed in row-based. the i-th row indices
                        # * are from rvec[i] to rvec[i+1]-1 */
        int *ridx; # /* index array row wise */
        double *rval; #   /* value array row wise */
        int *cvec; # /* matrix stored in column-based. the i-th column
                        # * indices are from cvec[i] to cvec[i+1]-1 */
        int *cidx; # /* index array column wise */
        double *cval; # /* value array column wise */
        int nnz;
        int bufsz; # /* current buffer size */
        int incsz; # /* the size to increase when needd */
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

cdef extern from "matfun.h":
    void * Malloc( int, char * );


def help():
    print("===================================================================")
    print("Incomplete multi-level QR factorization Parameters definition")
    print("===================================================================")
    print("mat: matrix stored in SparMat format -- see heads.h for details")
    print("mqrs     : pointer to a SparMQRS struct -- see heads.h for details")
    print("nomulti  : 0: IMQR, 1: IQR")
    print("nosort :1 if don't sort the graph A'*A according the degrees else 0")
    print("tol_level: tolerance for stopping MQR. Stop if:")
    print("       reduced size < tol_level * previous size")
    print("tol_cos  :   orthogonal tolerance: if:")
    print("    (v1,v2) < tol_cos*||v1||*||v2|| ")
    print("    vectors v1 and v2 are considered as orthogoanl")
    print("tol_drop :  dropping tolerance")
    print("lfil     : max # of elements allowed in each column")
    print("fp       : file pointer for error log (might be stdout )")
    print("on return:")
    print("ierr     = return value.")
    print("           ierr  = 0   --> successful return.")
    print("           ierr != 0   --> error")
    print("mqrs")
    print("===================================================================")
    print("work array: wrk at least m + n, iwrk at least 2m + 2n, iwrk[i] = 0")
    print("===================================================================")


cdef extern from "defs.h":
    int imqr( matptr mat, mqrsptr mqrs, int nomulti, int nosort,
          double tol_level, int maxlvl, double tol_cos,
          double tol_cos_adaptive, double tol_drop,
          int lfil, double *wrk, int *iwrk, FILE *fp );


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


cdef matDic(matptr F):
    print("beginging.......7")
    T = {}
    T["m"] = F.m
    T["n"] = F.n
    T["nnz"] = F.nnz
    T["rvec"], T["rval"], T["ridx"] = [], [], []
    T["cvec"], T["cval"], T["cidx"] = [], [], []
    # print("F.m  is ", F.m )
    # print("F.n  is ", F.n )
    # print("F.rvec size", sizeof(F.rvec) / sizeof(F.rvec[0])  )
    # print("F.rval size", sizeof(F.rval) / sizeof(F.rval[0])  )
    # print("F.ridx size", sizeof(F.ridx) / sizeof(F.ridx[0])  )
    # print("F.rvec", F.rvec[0])
    #
    # print("T[rval] is ", T["rval"])
    # for i in range(F.m):
    #     T["rvec"].append(F.rvec[i])
    # print("beginging.......8")
    # for i in range(F.nnz):
    #     T["rval"].append(F.rval[i])
    #
    # for i in range(F.nnz):
    #     T["ridx"].append(F.ridx[i])
    #
    #
    #
    # for i in range(F.n + 1):
    #     T["cvec"].append(F.cvec[i])
    #
    # for i in range(F.nnz):
    #     T["cval"].append(F.cval[i])
    #
    # for i in range(F.nnz):
    #     T["cidx"].append(F.cidx[i])
    return T


cdef mqrDic(mqrptr F):
    print("beginging.......6")
    T = {}
    T["n"] = F.n
    T["s"] = F.s
    T["perm"] = []
    T["iperm"] = []

    for i in range(F.n):
        T["perm"].append(F.perm[i])
    for i in range(F.n):
        T["iperm"].append(F.iperm[i])
    T["E"] = matDic(F.E)
    T["D_1"] = []
    for i in range(F.n):
        T["D_1"].append(F.D_1[i])
    return T

cdef mqrsDic(mqrsptr F):
    print("beginging.......5")
    T = {}
    T["n"] = F.n
    T["lvl"] = F.lvl
    amqrs = {}
    for i in range(F.lvl):
        amqrs[str(i)] = mqrDic(F.amqrs[i])
    T["amqrs"] = amqrs
    T["tm_lvl"] = F.tm_lvl
    T["tm_iqr"] = F.tm_iqr
    return T

def imqrPrecon(a_matrix, nomulti = 0, nosort = 0):
    inputs = Inputs()
    cdef matptr mat = NULL;
    cdef mqrsptr mqrs = NULL;
    cdef int ierr;
    print("beginging.......1")
    mqrs = <mqrsptr>Malloc( sizeof(SparMQRS), "imqrPrecon")
    mat = <matptr>Malloc( sizeof(SparMat), "imqrPrecon" )
    print("beginging.......2")


    AmatrixC = sparse.csc_matrix(a_matrix)
    AmatrixR = sparse.csr_matrix(a_matrix)

    mat.nnz = AmatrixC.nnz
    mat.m = AmatrixC.shape[0]
    mat.n = AmatrixC.shape[1]
    print("beginging.......3")
    print("mat.m is ", mat.m)
    print("mat.n is ", mat.n)

    mat.cvec = <int *>Malloc( (len(AmatrixC.indptr) ) * sizeof(int), "imqrPrecon" )
    mat.cidx = <int *>Malloc( (len(AmatrixC.indices) ) * sizeof(int), "imqrPrecon" )
    mat.cval = <double *>Malloc( (len(AmatrixC.data) ) * sizeof(double), "imqrPrecon" )

    mat.rvec = <int *>Malloc( (len(AmatrixR.indptr) ) * sizeof(int), "imqrPrecon" )
    mat.ridx = <int *>Malloc( (len(AmatrixR.indices) ) * sizeof(int), "imqrPrecon" )
    mat.rval = <double *>Malloc( (len(AmatrixR.data) ) * sizeof(double), "imqrPrecon" )
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

    inputs.tol_lfil = <int>(inputs.tol_lfil * <double>mat.nnz / <double>mat.n)
    if( inputs.tol_lfil < 1 ):
        inputs.tol_lfil  = 1
    cdef FILE *flog = stdout

    wk = <double *>Malloc( (mat.m + mat.n)*sizeof(double), "imqrPrecon")
    iwk = <int *>Malloc( 2*(mat.m + mat.n)*sizeof(int), "imqrPrecon")
    memset( iwk, 0, 2*(mat.m + mat.n)*sizeof(int) )
    tol_cos = inputs.tol_cos

    _to_return = {}
    print("beginging.......4")
    if 'amqrs' not in _to_return:
        # we want to calculate matrix inversion only once...
        # _to_return['inverted'] = \
        ierr = imqr( mat, mqrs, nomulti, nosort, inputs.tol_level, inputs.maxlvl,
                   inputs.tol_cos, inputs.tol_cos_adaptive, inputs.tol_drop, inputs.tol_lfil,
                   wk, iwk, flog )
        if ierr > 0:
            print("imqr error")
            return -1
        else:
            _to_return = mqrsDic(mqrs)
            # cval, cidx, cvec  = [], [], []
            # for i in range((mqrs.R.n + 1)):
            #     cvec.append(mqrs.R.cvec[i])
            # for i in range(mqrs.R.nnz):
            #     cval.append(mqrs.R.cval[i])
            #     cidx.append(mqrs.R.cidx[i])
            # _to_return["cval"] = cval
            # _to_return["cidx"] = cidx
            # _to_return["cvec"] = cvec
            # _to_return["n"] = mqrs.n
            # _to_return["lvl"] = mqrs.lvl
            # for i in range(mqrs.lvl):
            #
            # _to_return["amqrs_n"] = mqrs.lvl
            # R = sparse.csc_matrix((cval, cidx, cvec), shape = (mqrs.R.m, mqrs.R.n))
            # _to_return['Noninverted'] = mqrs.R.cval
            # R = sparse.csc_matrix((mqrs.R.cval, mqrs.R.cidx, mqrs.R.cvec), shape = (mqrs.R.m, mqrs.R.n))
            # _to_return['Noninverted'] = R.T * R
    # return np.matrix(_to_return['Noninverted'])
    # return _to_return['Noninverted']
    return _to_return



