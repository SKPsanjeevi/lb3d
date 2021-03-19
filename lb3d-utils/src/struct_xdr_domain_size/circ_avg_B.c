#include "main.h"


int 
circ_avg(double ***data, Params *params, CASF *casf)
{
   int               m;
   double            x, y, z, kdist;
   unsigned long     i, j, k, nx, ny, nz;

   nx = (unsigned long int) params->nx;
   ny = (unsigned long int) params->ny;
   nz = (unsigned long int) params->nz;

   for (m = 0 ; m < params->kradius ; m++) {
      casf->S[m]         = 0.0;
      casf->sum_sqS[m]   = 0.0;
      casf->errS[m]      = 0.0;
      casf->sum_shell[m] = 0;
   }

   /* Now loop through k-space, k = 2*PI*n/L , n: index, L: lattice
      side length */
   for (i = 1 ; i <= nx ; i++)                             
      for (j = 1 ; j <= ny ; j++) 
	 for (k = 1 ; k < nz ; k++) {  
/*JENS: Why k=1..nz-1 ? */

         /* BEWARE: Before, ONLY POSITIVE k_z
	    are in data[][][]. Now I use
	    a symmetry property in sf.c
	    to fill the rest. 
	 */
	    
	    /* For negative-frequency points, kdist is
	       calculated differently */
	    if(i > nx/2)  x = nx - i + 1; else x = i - 1;
	    if(j > ny/2)  y = ny - j + 1; else y = j - 1;
	    if(k > nz/2)  z = nz - k; else z = k - 1;
	    
	    kdist = sqrt(x * x + y * y + z * z);
	    
	    for (m = 0 ; m < params->kradius ; m++) 
	       if(kdist >= (double)(m - 0.5) && kdist < (double)(m +
								 0.5))
		  /* only if it's in the annulus */
	       {
		  casf->S[m] += data[i][j][k];       
	//	  casf->sum_sqS[m] += data[i][j][k] * data[i][j][k];
		  casf->sum_shell[m] += 1;
	       }
	 }   

   for (m = 0 ; m < params->kradius ; m++) 
      casf->S[m] /= (double)casf->sum_shell[m];

   for (i = 1 ; i <= nx ; i++)
      for (j = 1 ; j <= ny ; j++)
         for (k = 1 ; k < nz ; k++) 
	    for (m = 0 ; m < params->kradius ; m++) {
 	       if(i > nx/2)  x = nx - i + 1; else x = i - 1;
               if(j > ny/2)  y = ny - j + 1; else y = j - 1;
               if(k > nz/2)  z = nz - k; else z = k - 1;
	       kdist = sqrt(x * x + y * y + z * z);
	       if(kdist >= (double)(m - 0.5) &&
	          kdist < (double)(m + 0.5)) {
	          casf->errS[m] += (data[i][j][k]-casf->S[m])*(data[i][j][k]-casf->S[m]);
	       }
            }

   for (m = 0 ; m < params->kradius ; m++) 
      if(casf->sum_shell[m] > 1)
         casf->errS[m] = sqrt(casf->errS[m] /(double)(casf->sum_shell[m] - 1) / (double)casf->sum_shell[m]);
      else
         casf->errS[m] = 0.0;

 /*
   printf("\n");
   for (m = 0 ; m < params->kradius ; m++) 
      printf("\nAvg of S in shell %lu = %15.7g +- %15.7g  |  %lu points", 
             m, casf->S[m], casf->errS[m], casf->sum_shell[m]);
 */

   return(1);
}  
  

