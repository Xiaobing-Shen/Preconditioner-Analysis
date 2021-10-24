/*---------------------------------------------------------------------------*
 * main test driver for MIQR                                                 *
 *---------------------------------------------------------------------------*
 * Na Li, Jul 17, 2004                                                       *
 *                                                                           *
 * Report bugs / send comments to: nli@cs.umn.edu                            *
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "heads.h"
#include "defs.h"
#include "matfun.h"

int usage()
{
  printf( "======================================================\n" );
  printf( "Usage: mqr.ex [-np] [-nm] [-ns] [-nr] [-l] [-k] [-o] [-h]\n" );
  printf( "======================================================\n" );
  printf( "Options:\n" );
  printf( "-np: disable preconditioner when solving the system\n" );
  printf( "-nm: disable multi-level preconditioner (IQR instead)\n" );
  printf( "-ns: do NOT sort the graph A'*A according the degrees\n" );
  printf( "-nr: do NOT read rhs from file. use artificial rhs\n" );
  printf( "-l:  output the reduced size for each level\n" );
  printf( "-k:  skip iterations\n" );
  printf( "-o:  output the solutions\n" );
  printf( "-h:  display this help\n" );
  printf( "======================================================\n" );
  printf( "make sure file 'inputs' and file 'matfile' are in the \n" );
  printf( "same directory as the executable mqr.ex \n" );
  printf( "======================================================\n" );
  return 0;
}

int main( int argc, char **argv )
{
  int ierr = 0;

/*-------------------------------------------------------------------
 * options
 *-----------------------------------------------------------------*/
  int argvlen, skip_its = 0, output_levels = 0, output_sol = 0;
  int noprec = 0, nomulti = 0, nosort = 0, manu_rhs = 0;
  
  matptr mat = NULL;
  mqrsptr mqrs = NULL;
  double tol_cos;

  int nrhs;
  double *sol = NULL, *guess = NULL, *rhs = NULL;
  double *x = NULL, *y = NULL;

  /* work arrays */
  double *wk = NULL;
  int *iwk = NULL;
  
  FILE *flog = stdout, *fmat = NULL;
  io_t io;
  clock_t clk;
  
  int nmat, numat, iparam, i;
  char line[MAX_LINE];
  
/*-------------------------------------------------------------------
 * process options
 *-----------------------------------------------------------------*/
  i = 1;
  while( i < argc ) {
    if( argv[i][0] == '-' ) {
      switch( argv[i][1] ) {
      case 'k':
	skip_its = 1;
	break;
      case 'l':
        output_levels = 1;
        break;
      case 'n':
        argvlen = strlen( argv[i] );
        if( argvlen >= 3 && argv[i][2] == 'p' ) noprec = 1;
        else if( argvlen >= 3 && argv[i][2] == 'm' ) nomulti = 1;
        else if( argvlen >= 3 && argv[i][2] == 's' ) nosort = 1;
        else if( argvlen >= 3 && argv[i][2] == 'r' ) manu_rhs = 1;
        else {
	  printf( "Invalid option...\n" );
	  usage();
	  return 0;
        }
        break;
      case 'o':
        output_sol = 1;
        break;
      case 'h':
	usage();
	return 0;
      default:
	printf( "Invalid option...\n" );
	usage();
	return 0;
      }
    }
    i++;
  }
  
/*-------------------------------------------------------------------
 * reads matrix A from Harwell Boeing file
 *
 * solves A^T * A * x = rhs using MQR preconditioned CGNR
 *
 *-----------------------------------------------------------------*/
  
  memset( &io, 0, sizeof(io) );
/*-----------------------------------------------------------------*/
  if( read_inputs( "inputs", &io ) != 0 ) {
    fprintf( flog, "Invalid inputs file...\n" );
    goto ERROR_HANDLE;
  }
/*-----------------------------------------------------------------*/
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    goto ERROR_HANDLE;
  }
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  if( ( numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    goto ERROR_HANDLE;
  }
/* LOOP THROUGH MATRICES ------------------------------------------*/
  for( nmat = 1; nmat <= numat; nmat++ ) {
    if( get_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile...\n" );
      goto ERROR_HANDLE;
    }
    fprintf( flog, "\n========================================\n" );
    fprintf( flog, "Solving MATRIX: %s...\n", io.HBnameF );
    fprintf( flog, "========================================\n" );
/* Read in matrix and allocate memory------------------------------*/
    mat = (matptr)Malloc( sizeof(SparMat), "main" );
    ierr = readhbc_c( mat, &io, manu_rhs );
    if( ierr != 0 ) {
      fprintf( flog, "readhb_c error = %d\n", ierr );
      goto ERROR_HANDLE;
    }

    sprintf( io.outfile, "OUTPUT/%s.out", io.HBnameF );
    if( NULL == ( io.fout = fopen( io.outfile, "w" ) ) ) {
      fprintf(flog,"Can't open output file %s...\n", io.outfile);
      fprintf(flog,"Make sure sub-directory 'OUTPUT' exists.\n" );
      goto ERROR_HANDLE;
    }

    io.lfil = (int)(io.tol_lfil * (double)io.nnz / (double)io.ndim);
    if( io.lfil < 1 ) io.lfil = 1;
    output_header( &io );
    
    tol_cos = io.tol_cos;

    x = (double *)Malloc( io.ndim * sizeof(double), "main" );
    y = (double *)Malloc( max(io.ndim,io.mdim)*sizeof(double), "main" );
    wk = (double *)Malloc( (io.mdim+io.ndim)*sizeof(double), "main" );
    iwk = (int *)Malloc( 2*(io.mdim+io.ndim)*sizeof(int), "main" );
    memset( iwk, 0, 2*(io.mdim+io.ndim)*sizeof(int) );

/* LOOP THROUGH PARAMETERS ---------------------------------------*/
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
        fprintf( flog, "**** test # %d ****\n", iparam );

         /* Preconditioning */
        mqrs = (mqrsptr)Malloc( sizeof(SparMQRS), "main" );
        clk = clock();
        ierr = imqr( mat, mqrs, nomulti, nosort, io.tol_level, io.maxlvl,
                   tol_cos, io.tol_cos_adaptive, io.tol_drop, io.lfil,
                   wk, iwk, flog );
        clk = clock() - clk;
        io.tm_p = (double)clk / (double)CLOCKS_PER_SEC;
        io.tm_total = io.tm_p;
        if( ierr != 0 ) {
            fprintf( flog, "*** mqr error, ierr != 0 ***\n" );
	        goto ERROR_HANDLE;
        }
        io.fillin = nnzMQRS( mqrs );
        fprintf( flog, "  MQR: level = %d, time = %8.1e (s), mem used = %d\n",
               mqrs->lvl, io.tm_p, io.fillin );
      
        output_prec( &io, mqrs, output_levels, iparam );

        if( skip_its ) {
	        io.its = -1;
	        io.tm_i = -1;
	        io.enorm = -1;
	        io.rnorm = -1;
	        goto NEXT_PARA;
        }

        if( cond_est( mat, mqrs, y, x, flog ) != 0 ) {
	        fprintf( flog, "  Not attempting iterative solution.\n" );
	        fprintf( io.fout, "Not attempting iterative solution.\n" );
	        io.its = -1;
	        io.tm_i = -1;
	        io.enorm = -1;
	        io.rnorm = -1;
	        goto NEXT_PARA;
        }

        /* Iterations */
        for( nrhs = 0; nrhs < io.nrhs; nrhs++ ) {
            sol = io.sol + nrhs * io.ndim;
            guess = io.guess + nrhs * io.ndim;
            rhs = io.rhs + nrhs * io.mdim;
            for( i = 0; i < io.ndim; i++ ) {
                x[i] = guess[i];
            }
            io.its = io.maxits;
            clk = clock();
            ierr = pcgnr( mat, mqrs, rhs, x, io.tol_its, &io.its, noprec );
            clk = clock() - clk;
            io.tm_i = (double)clk / (double)CLOCKS_PER_SEC;
            io.tm_total += io.tm_i;
            fprintf( flog, "  PCGNR: time = %8.1e (s), ", io.tm_i );
            if( ierr == 0 ) {
                fprintf( flog, "Converged in %d iterations.\n", io.its );
            } else {
                fprintf( flog, "Not Converge in %d iterations.\n", io.maxits );
            }

            /* calculate error norm */
            if( sol != NULL ) {
                io.enorm = vec_nrm2_diff( io.ndim, x, sol );
            } else {
                io.enorm = -1.0;
            }

            /* calculate residual norm */
            matxvec( mat, x, y );
            vec_sub( io.mdim, y, rhs, y );
            matTxvec( mat, y, wk );
            io.rnorm = vec_nrm2( io.ndim, wk );

            output_result( &io, iparam, nrhs );
            if( output_sol ) {
                vec_output( io.ndim, x, "x", io.fout );
            }
        }
      
NEXT_PARA:
      tol_cos += io.tol_delta;
      cleanMQRS( mqrs );
    }
    
    fclose( io.fout );
    
    cleanMat( mat );
    if( x ) free( x );
    if( y ) free( y );
    if( wk ) free( wk );
    if( iwk ) free( iwk );
    if( io.rhs ) { free( io.rhs ); io.rhs = NULL; }
    if( io.guess ) { free( io.guess ); io.guess = NULL; }
    if( io.rhs ) { free( io.rhs ); io.rhs = NULL; }
  }
  
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  
  return 0;
  
ERROR_HANDLE:
  exit( -1 );
}
