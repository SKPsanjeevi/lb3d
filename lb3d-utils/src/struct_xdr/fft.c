#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"
#include "nrutil.h"

extern double **dmatrix(long, long, long, long);
extern void   free_dmatrix(double **, long, long, long, long);
extern void   rlft3_d(double ***, double **, unsigned long, unsigned long,
		      unsigned long, int);
extern int    dump_second_field, dump_third_field;

double        calc_avg_field(double ***, unsigned long, unsigned long,
			     unsigned long, Stats *);
int           fft(double ***, Params *, Stats *);


double
calc_avg_field(double ***data, unsigned long nx, unsigned long ny,
	       unsigned long nz, Stats *stats)
{
  int         i, j, k;
  double      field, avg;
  /*  FILE        *fp; */

  field = 0.;

  /* fp = fopen("blowup.dat", "w"); */

  for(i = 1; i <= nx; i++)
    for(j = 1; j <= ny; j++)
      for(k = 1; k <= nz; k++)
	{
	  field += data[i][j][k];
	  /* fprintf(fp, "%d %d %d %e %e\n",  
	     i, j, k, data[i][j][k], field); fflush(fp); 
	  */
	}

  /* fclose(fp); */

  avg = field / (double)(nx * ny * nz);
  if(dump_second_field)
     stats->avg2 = avg;
  else if(dump_third_field)
     stats->avg3 = avg;
  else
     stats->avg = avg;

  fflush(stdout);

  return(avg);

}




int
fft(double ***data, Params *params, Stats *stats)
{
   double             avg_field;
   double             **speq;
   unsigned long int  i, j, k, nx, ny, nz;
   FILE               *fp;
   static int         counter = 0;
   int                cond;

   cond=0;
/*
   cond=((counter==11)&&(!dump_second_field));
*/
   if(cond)
      if((fp = fopen("/home3/nelido/tmp/fft_sur_t4000_DOUBLE_-O0.dat", "w"))==NULL) {
	 fprintf(stderr,"\nError opening fft dumpfile\n");
	 exit(1);
      }

   
   /* TAYLOR LATTICE SIZE PARAMETERS' TYPE FOR rlft3_d() */
   nx = (unsigned long int) params->nx;
   ny = (unsigned long int) params->ny;
   nz = (unsigned long int) params->nz;
   
   speq = dmatrix(1, nx, 1, 2 * ny);
   
   avg_field = calc_avg_field(data, nx, ny, nz, stats);
   
   for(i = 1; i <= nx; i++)
      for(j = 1; j <= ny; j++)
	 for(k = 1; k <= nz; k++)  
	 {
	    data[i][j][k] -= (double)avg_field;
	 }
   
   rlft3_d(data, speq, nx, ny, nz, +1);
   
   for(i = 1; i <= nx; i++)
      for(j = 1; j <= ny; j++)
	 for(k = 1; k <= nz; k++) {           /* DONE IN-PLACE!! */
	    if(k <= (nz / 2))
	       data[i][j][k] = (data[i][j][2 * k - 1] *
				data[i][j][2 * k - 1] + 
				data[i][j][2 * k] *
				data[i][j][2 * k]) / 
		  (double)(nx * ny * nz);
	    else 
	       data[i][j][k] = data[i][j][nz - k + 1];   
	    /* Symmetry
	       property */
	    
	    
	    if(cond) {
	       fprintf(fp, "%lu %lu %lu %15.7e\n", i, j, k,
		       data[i][j][k]);
	       fflush(fp);
	    }
	    
	 }
   
   if(cond) {
      printf("\nWrote fft to disc...");
      fclose(fp);
   }
   
   free_dmatrix(speq, 1, nx, 1, 2 * ny);
   
   counter++;
   
   return(1);
}














