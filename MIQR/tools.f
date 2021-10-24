c-----------------------------------------------------------------------
c     some routines extracted/ modified from SPARSKIT2 + one from blas
c-----------------------------------------------------------------------
      subroutine readmtc (nmax,nzmax,job,fname,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine.
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Jul 20, 1998 by Irene Moulitsas
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax   =  max column dimension  allowed for matrix. The array ia should
c           be of length at least ncol+1 (see below) if job.gt.0
c nzmax  = max number of nonzeros elements allowed. the arrays a,
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c
c job    = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c                     rhs may contain initial guesses and exact
c                     solutions appended to the actual right hand sides.
c                     this will be indicated by the output parameter
c                     guesol [see below].
c
c fname = name of the file where to read the matrix from.
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c on return:
c----------
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs
c          is provided then it is rest to job=2 or job=1 depending on
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a      = the a matrix in the a, ia, ja (column) storage format
c ja     = column number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c          if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c          if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol   = number of columns in matrix
c nnz    = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c
c title  = character*72 = title of matrix test ( character a*72).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages
c         * ierr  =  0 means that  the matrix has been read normally.
c         * ierr  =  1 means that  the array matrix could not be read
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read
c         because nnz .gt. nzmax
c         * ierr  =  3 means that  the array matrix could not be read
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...)
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of
c         rhs as specified by the input value of nrhs is not
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c
c Notes:
c-------
c 1) This routine can be interfaced with the C language, since only
c    the name of the  file needs to be passed and no iounti number.
c
c 2) Refer to the  documentation on  the Harwell-Boeing formats for
c    details on the format assumed by readmt.
c    We summarize the format here for convenience.
c
c    a) all  lines in  inout are  assumed to be 80  character long.
c    b) the file  consists of a header followed by the block of the
c       column start  pointers followed  by  the  block  of the row
c       indices,  followed  by  the  block  of the  real values and
c       finally the  numerical  values of  the right-hand-side if a
c       right hand side is supplied.
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first  line  contains  the  title  (72  characters  long)
c         followed  by  the  8-character  identifier (name  of  the
c         matrix, called key) [ A72,A8 ]
c       * second line  contains the number of lines for each of the
c         following data blocks (4 of them) and the total number of
c         lines excluding the header.  [5i4]
c       * the   third  line  contains  a  three   character  string
c         identifying the type of  matrices as they  are referenced
c         in  the Harwell  Boeing documentation [e.g., rua, rsa,..]
c         and the number of rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth  line contains the variable fortran format for
c         the following data blocks. [2A16,2A20]
c       * The fifth  line is  only present if  right-hand-sides are
c         supplied. It  consists  of  three  one  character-strings
c         containing the  storage  format for the  right-hand-sides
c         ('F'= full,'M'=sparse=same as matrix), an  initial  guess
c         indicator  ('G' for yes),  an  exact  solution  indicator
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices.  [A3,11X,2I14]
c     d) The three  following blocks follow the header as described
c        above.
c     e) In case the right hand-side are in sparse formats then the
c        fourth  block  uses the  same  storage  format as  for the
c        matrix to  describe  the NRHS right  hand  sides provided,
c        with a column being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     &     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     &     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax)
      real*8 a(nzmax), rhs(*)
      character fname*100
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c
      iounit=15
      open(iounit,file = fname)
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd,
     &     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt,
     &     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c
c     anything else to read ?
c
      if (job .le. 0) goto 12
c     ---- check whether matrix is readable ------
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) goto 12
c     ---- read pointer and row numbers ----------
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  goto 12
c     --- and if available -----------------------
      if (valcrd .le. 0) then
         job = 1
         goto 12
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ----------------
      if (job .le. 2)  goto 12
c     --- and if available -----------------------
      if ( rhscrd .le. 0) then
         job = 2
         goto 12
      endif
c
c     --- read right-hand-side.--------------------
c
      if (rhstyp(1:1) .eq. 'M') then
         ierr = 4
         goto 12
      endif
c
      guesol = rhstyp(2:3)
c
      nvec = 1
      size = nrhs * nrow
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        nvec=nvec+1
        size = size + nrhs * ncol
      endif
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        nvec=nvec+1
        size = size + nrhs * ncol
      endif
c
      len = nrhs*nrow
c
      if (size .gt. lenrhs) then
        ierr = 5
        goto 12
      endif
c
c     read right-hand-sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c
c     read initial guesses if available
c
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        next = iend + 1
        iend = iend+ nrhs*ncol
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
c     read exact solutions if available
c
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        next = iend + 1
        iend = iend+nrhs*ncol 
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
 12   close(iounit)
      return
c---------end-of-readmt_c-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
        subroutine qsplit2(a,m,ind,nnz,ncut)
        real*8 a(m)
        integer m, nnz, ind(nnz), ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:m). is a real array
c     ind is a zero-based index array
c     only a(ind(i)+1) .ne. 0, i = 1,2,...,nnz
c     on output ind(nnz) is permuted such that its elements satisfy:
c
c     abs(a(ind(i)+1)) .ge. abs(a(ind(ncut)+1)) for i .lt. ncut and
c     abs(a(ind(i)+1)) .le. abs(a(ind(ncut)+1)) for i .gt. ncut
c
c     ind(1:nnz) is an integer array which will be permuted
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = nnz
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(ind(mid)+1))
        do 2 j=first+1, last
           if (abs(a(ind(j)+1)) .gt. abskey) then
              mid = mid+1
c     interchange
              itmp = ind(mid)
              ind(mid) = ind(j)
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit2-----------------------------------------
c-----------------------------------------------------------------------
        end
c-----------------------------------------------------------------------
        subroutine qsplit(a,ind,n,ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end

