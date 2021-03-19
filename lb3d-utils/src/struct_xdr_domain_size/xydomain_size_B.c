#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"


int 
xydomain_size(double ***data, Params *params, XYAVGSF *xyavgsf)
{
   double            x, y, z, norm, vol, npre;
   unsigned long     i, j, k, nx, ny, nz;
   
   nx            = (unsigned long) params->nx;
   ny            = (unsigned long) params->ny;
   nz            = (unsigned long) params->nz;
   vol           = (double)nx*ny*nz;
   xyavgsf->xxS  = 0.0;
   xyavgsf->yyS  = 0.0;
   xyavgsf->zzS  = 0.0;
   norm          = 0.0;

   /* Now loop through k-space, k = 2*PI*n/L , n: index, L: lattice
      side length */
   for (i = 1 ; i <= nx ; i++)                             
      for (j = 1 ; j <= ny ; j++) 
	 for (k = 1 ; k < nz ; k++) {  
	    /* 
	     *  BEWARE: Before, ONLY POSITIVE k_z are in
	     * data[][][]. Now I use a symmetry property in sf.c to
	     * fill the rest. 
	     * For negative-frequency points, kdist is calculated
	     * differently */


	    if(i > nx/2)  x = (double)nx - i + 1; else x = (double)i - 1;
	    if(j > ny/2)  y = (double)ny - j + 1; else y = (double)j - 1;
	    if(k > nz/2)  z = (double)nz - k;     else z = (double)k - 1;

	   /* x               = DPI / (2 * params->kradius) * x;
	    y               = DPI / (2 * params->kradius) * y;
	    z               = DPI / (2 * params->kradius) * z;
	    x               = DPI * x;
	    y               = DPI * y;
	    z               = DPI * z;
	    */
	    /* There should be a DPI here and further down, but they */
	    /* cancel out */
	    x	= x/(double)nx;
	    y	= y/(double)ny;
	    z	= z/(double)nz;

	    norm           += data[i][j][k];
	    xyavgsf->xxS   += x * x * data[i][j][k];       
	    xyavgsf->yyS   += y * y * data[i][j][k];
	    xyavgsf->zzS   += z * z * data[i][j][k];
	    
//
// Comment:
// The following two lines are computed for potential future 
// applications. I'm not implementing the next function called 
// by caller.c (dump_to_files.c) to write anything derived from 
// them. To use them, just uncomment them.      --NGS
//
//	    xyavgsf->absxS += fabs(x) * data[i][j][k];
//	    xyavgsf->absyS += fabs(y) * data[i][j][k];
//	    xyavgsf->abszS += fabs(z) * data[i][j][k];
	 }
   
   if(norm==0) 
   {
      fprintf(stderr, 
	      "\nxydomain_size_B.c: division by norm==0 !\n");
      exit(0);
   }
   xyavgsf->xxS /= norm;
   xyavgsf->yyS /= norm;
   xyavgsf->zzS /= norm;

//   xyavgsf->absxS /= norm;
//   xyavgsf->absyS /= norm;
//   xyavgsf->abszS /= norm;

   for (i = 1 ; i <= nx ; i++)
      for (j = 1 ; j <= ny ; j++)
         for (k = 1 ; k < nz ; k++) {

	    if(i > nx/2)  x = (double)nx - i + 1; else x = (double)i - 1;
	    if(j > ny/2)  y = (double)ny - j + 1; else y = (double)j - 1;
	    if(k > nz/2)  z = (double)nz - k;     else z = (double)k - 1;
	    xyavgsf->xxSerr +=
	       (data[i][j][k] - xyavgsf->xxS)*(data[i][j][k] - xyavgsf->xxS);
	    xyavgsf->yySerr +=
	       (data[i][j][k] - xyavgsf->yyS)*(data[i][j][k] - xyavgsf->yyS);
	    xyavgsf->zzSerr +=
	       (data[i][j][k] - xyavgsf->zzS)*(data[i][j][k] - xyavgsf->zzS);

//
// See comment above 
//	    xyavgsf->absxSerr +=
//	       (data[i][j][k] - xyavgsf->absxS)*
//	       (data[i][j][k] - xyavgsf->absxS);
//	    xyavgsf->absySerr +=
//	       (data[i][j][k] - xyavgsf->absyS)*
//	       (data[i][j][k] - xyavgsf->absyS);

	 }
 
   xyavgsf->xxSerr = sqrt(xyavgsf->xxSerr /(vol-1.0) /vol);
   xyavgsf->yySerr = sqrt(xyavgsf->yySerr /(vol-1.0) /vol);
   xyavgsf->zzSerr = sqrt(xyavgsf->zzSerr /(vol-1.0) /vol);

//
// See comment above  
//   xyavgsf->absxSerr = sqrt(xyavgsf->absxSerr /(vol-1.0) /vol);
//   xyavgsf->absySerr = sqrt(xyavgsf->absySerr /(vol-1.0) /vol);

// Prefactor to normalize domain sizes. We divide by sqrt(3) to be able to
   // compare to old cubic data.
 
   //npre = sqrt((double)(nx*nx+ny*ny+nz*nz))/sqrt(3) ;
   //npre = (double)(nx*nx+ny*ny+nz*nz)/3 ;
   //npre = (double)(nx+ny+nz)/3;
   //npre = 1/sqrt(3) ;
   xyavgsf->Rx    = 1/sqrt(xyavgsf->xxS);
   xyavgsf->Ry    = 1/sqrt(xyavgsf->yyS);
   xyavgsf->Rz    = 1/sqrt(xyavgsf->zzS);

   if((xyavgsf->xxS==0) || 
      (xyavgsf->yyS==0) ||
      (xyavgsf->zzS==0))
   {
      fprintf(stderr, 
	      "\nxydomain_size_B.c: division by xxS or yyS==0 or zzS==0 in computing the error !\n");
      exit(0);  
   }

   xyavgsf->Rxerr = xyavgsf->xxSerr * PI /pow(xyavgsf->xxS, 1.5);
   xyavgsf->Ryerr = xyavgsf->yySerr * PI /pow(xyavgsf->yyS, 1.5);
   xyavgsf->Rzerr = xyavgsf->zzSerr * PI /pow(xyavgsf->zzS, 1.5);

   return(1);
}  
  




