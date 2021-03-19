/*
Compute the integral:   
$\sigma=\int_{-\Infty}^{+\Infty}(P_N(z)-P_T(z))dz$,
where the infinite limits are substituted by the lattice limits,         
on assuming that for the lattice size considered an 'asymptotic'
limit to zero is achieved. Integration formula is taken from 
Press et al: Numerical Recipes in C. 2nd ed. CUP, 1995. (p.134)
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"


/********************************************************
 ** ARGS:   ts                     : current timestep
 **         params                 
 **         profile[1..params->nz]
 **
 **
 ********************************************************
 */


int 
plinteg(unsigned long ts, Params *params, float *profile, float *sigma) 
{
   int         ix, iz, iy;
   float       factor, sigma_previous, rel_diff;
   FILE        *fp;
   static int  count=0;

   sigma_previous = *sigma;
   *sigma = 0.0;

   if(VERBOSE)
	printf("\nplinteg.c:params->nz=%d", params->nz);

   if(params->direction_z) {
      //for(iz = 0; iz < params->nz; iz++) {
	for(iz = 0; iz < params->nz; iz++) {	
	 if((iz == 0) || (iz == (params->nz-1))) factor = 3./8.;
	 else if((iz == 1) || (iz == (params->nz-2))) factor = 7./6.;
	 else if((iz == 2) || (iz == (params->nz-3))) factor = 23./24.;
	 else factor = 1.;
	 
	 *sigma += factor * profile[iz];
      }
   } else if(params->direction_x) {
      for(ix = 0; ix < params->nx; ix++) {
	 if((ix == 0) || (ix == (params->nx-1))) factor = 3./8.;
	 else if((ix == 1) || (ix == (params->nx-2))) factor = 7./6.;
	 else if((ix == 2) || (ix == (params->nx-3))) factor = 23./24.;
	 else factor = 1.;
	 
	 *sigma += factor * profile[ix];      
      }
   } else if(params->direction_y) {
      for(iy = 0; iy < params->ny; iy++) {
	 if((iy == 0) || (iy == (params->ny-1))) factor = 3./8.;
	 else if((iy == 1) || (iy == (params->ny-2))) factor = 7./6.;
	 else if((iy == 2) || (iy == (params->ny-3))) factor = 23./24.;
	 else factor = 1.;
	 
	 *sigma += factor * profile[iy];      
      }
   }
   
   *sigma /= 2.;  /* DIVIDE BY 2 BECAUSE THERE ARE 2 INTERFACES */
   
   rel_diff = (*sigma - sigma_previous)/sigma_previous;

   /* Write onto stdout */
   if(count == 0) {
      if(remove(params->outfile) == 0) {
	 printf("\nAlready existing %s removed...\n",
		params->outfile);
	 count++;
      }
   }
   if((fp = fopen(params->outfile, "a"))==NULL) {
      printf("\nFile %s can't be opened\n", params->outfile);
      exit(0);
   }
   fprintf(fp, "%lu\t%16.8f\t%16.8f\n", ts, *sigma, rel_diff * 100.);
   fclose(fp);

   return(1);
}





