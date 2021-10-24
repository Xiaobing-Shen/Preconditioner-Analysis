/*---------------------------------------------------------------------------*
 * auxf.c of MQR
 *---------------------------------------------------------------------------*
 * Na Li, Jul 12, 2004                                                       *
 *                                                                           *
 * Report bugs / send comments to: nli@cs.umn.edu                            *
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "heads.h"
#include "defs.h"

void *Malloc( int nbytes, char *msg )
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL)
  {
    fprintf( stderr, "%s:Not enough mem for %d bytes\n", msg, nbytes );
    exit( -1 );
  }

  return ptr;
}

void *Realloc( void *buf, int nbytes, char *msg )
{
  void *ptr;

  ptr = (void *)realloc( buf, nbytes );
  if (ptr == NULL)
  {
    fprintf( stderr, "%s:Not enough mem for %d bytes\n", msg, nbytes );
    exit( -1 );
  }

  return ptr;
}

int read_integer( FILE *fp )
{
  char line[MAX_LINE] = "", *p1, *p2;
  fgets( line, MAX_LINE, fp );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  return atoi( p1 );
}

double read_double( FILE *fp )
{
  char line[MAX_LINE] = "", *p1, *p2;
  fgets( line, MAX_LINE, fp );
  for( p1 = line; ' ' == *p1; p1++ );
  for( p2 = p1; ' ' != *p2; p2++ );
  *p2 = '\0';
  return atof( p1 );
}

int read_inputs( char *in_file, io_t *pio )
{
  FILE *finputs;
  if( NULL == ( finputs = fopen( in_file, "r" ) ) )
  {
    return -1;
  }
  pio->nparam = read_integer( finputs );
  pio->tol_level = read_double( finputs );
  pio->maxlvl = read_integer( finputs );
  pio->tol_cos = read_double( finputs );
  pio->tol_cos_adaptive = read_double( finputs );
  pio->tol_delta = read_double( finputs );
  pio->tol_drop = read_double( finputs );
  pio->tol_lfil = read_double( finputs );
  pio->maxits = read_integer( finputs );
  pio->tol_its = read_double( finputs );
  
  fclose( finputs );
  
  return 0;
}

int get_matrix_info( FILE *fmat, io_t *pio )
{
  char line[MAX_LINE], *p1, *p2;

  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  for( p1 = line; '\'' != *p1; p1++ );
  p1++;
  for( p2 = p1; '\'' != *p2; p2++ );
  *p2 = '\0';
  strcpy( pio->Fname, p1 );
  
  for( p1 = p2+1; '\'' != *p1; p1++ );
  p1++;
  for( p2 = p1; '\'' != *p2; p2++ );
  *p2 = '\0';
  strcpy( pio->HBnameF, p1 );
  
  return 0;
}

int output_header( io_t *pio )
{
  FILE *f = pio->fout;

  fprintf( f, "========================================================\n" );
  fprintf( f, "|                  MATRIX INFORMATION                  |\n" );
  fprintf( f, "========================================================\n" );
  fprintf( f, "| Name | %-16s | Type    | %-16s |\n",
           pio->HBnameF, pio->type );
  fprintf( f, "--------------------------------------------------------\n" ); 
  fprintf( f, "| Rows |%11d | Cols |%11d | nRHS |%6d |\n",
           pio->mdim, pio->ndim, pio->nrhs );
  fprintf( f, "--------------------------------------------------------\n" ); 
  fprintf( f, "| nnz |%12d | drptol |%7.4f | lvltol |%6.2f |\n",
           pio->nnz, pio->tol_drop, pio->tol_level );
  fprintf( f, "========================================================\n" );
  fflush( f );
  return 0;
}

int output_prec( io_t *pio, mqrsptr mqrs, int output_levels, int nparam )
{
  FILE *f = pio->fout;
  double tol_cos = pio->tol_cos + pio->tol_delta * (nparam-1);
  double ratio = (double)pio->fillin / (double)pio->nnz;
  fprintf( f, "\nTest NO. (%d)\n", nparam );
  fprintf( f, "========================================================\n" );
  fprintf( f, "|  ortho-tol  |    lfil    |  rsdl R dim |  nnzMQR/nnz |\n" );
  fprintf( f, "--------------------------------------------------------\n" );
  fprintf( f, "|%12.4f |%11d |%12d |%12.3f |\n",
           tol_cos, pio->lfil, mqrs->R->n, ratio );
  fprintf( f, "--------------------------------------------------------\n" );
  fprintf( f, "|  # of lvls  | Mul-Lvl-T  |    IQR-T    |    Prec-T   |\n" );
  fprintf( f, "--------------------------------------------------------\n" );
  fprintf( f, "|%12d |%11.2f |%12.2f |%12.2f |\n", mqrs->lvl,
           mqrs->tm_lvl, mqrs->tm_iqr, pio->tm_p );
  if( output_levels ) {
    int i;
    fprintf( f, "========================================================\n" );
    fprintf( f, "| Size reduced at each level:                          |\n" );
    for( i = 0; i < mqrs->lvl; i++ ) {
      if( i % 6 == 0 ) fprintf( f, "|" );
      fprintf( f, "%8d ", mqrs->amqrs[i]->s );
      if( (i+1) % 6 == 0 ) fprintf( f, "|\n" );
    }
    fprintf( f, "\n" );
  }
  fprintf( f, "========================================================\n" );
  return 0;
}

int output_result( io_t *pio, int iparam, int nrhs )
{
  FILE *f = pio->fout;
  if( nrhs == 0 ) {
    fprintf(f, "| nrhs |  Its  |   Its-T   |   Err Norm  |  Rsdl Norm  |\n" );
    fprintf(f, "-------------------------------------------------------\n" );
  }
  fprintf( f, "|%5d |%6d |%10.2f |%12.4g |%12.4g |\n", nrhs + 1,
           pio->its, pio->tm_i, pio->enorm, pio->rnorm );
  if( nrhs == pio->nrhs - 1 ) {
    fprintf(f, "========================================================\n");
    fprintf(f, "| Total time: %-10.2f                               |\n",
            pio->tm_total );
    fprintf(f, "========================================================\n");
  }
  fflush( f );

  return 0;
}

