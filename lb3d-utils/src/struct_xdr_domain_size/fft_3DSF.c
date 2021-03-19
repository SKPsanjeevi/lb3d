#include "main.h"
#include "nrutil.h"
#include <rpc/rpc.h>

extern double **dmatrix(long, long, long, long);
extern void   free_dmatrix(double **, long, long, long, long);
extern void   rlft3_d(double ***, double **, unsigned long, unsigned long,
		      unsigned long, int);
extern int    dump_second_field;

double        calc_avg_field(double ***, unsigned long, unsigned long,
			     unsigned long, Stats *);
int           fft(double ***, Params *, Stats *, unsigned long);


double
calc_avg_field(double ***data, unsigned long nx, unsigned long ny,
	       unsigned long nz, Stats *stats)
{
  int         i, j, k;
  double      field, avg;
  float	      *ap;
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
  else
     stats->avg = avg;

  fflush(stdout);

  return(avg);

}




int
fft(double ***data, Params *params, Stats *stats, unsigned long timestep)
{
   double             avg_field;
   double             **speq;
   float	      ap,sumz,sumy,sumx;
   long int           x, y, z;
   unsigned long int  i, j, k, ii, jj, kk, nx, ny, nz,xi;
   char               aux[128] = "\0";
   char               sf3d_file[128] = "\0";
   char               sf3d_xdrfile[128] = "\0";
   FILE               *fp,*fp2,*fp3,*fp4,*outf,*outf2;
   static int         counter = 0;
   XDR		      xdrs,xdrs2;


   strcpy(aux, params->outfile);
   //sprintf(sf3d_file, "%s_t%0"TIMECHARSTR"lu.sf3d",aux,timestep);
   //if((fp = fopen(sf3d_file, "w"))==NULL) {
   //   fprintf(stderr,"\nError opening sf3d_file\n");
   //   exit(1);
   //}
   sprintf(sf3d_file, "%s_t%0"TIMECHARSTR"lu.sf3dzsum",aux,timestep);
   if((fp2 = fopen(sf3d_file, "w"))==NULL) {
      fprintf(stderr,"\nError opening sf3d_file zsum\n");
      exit(1);
   }
   sprintf(sf3d_file, "%s_t%0"TIMECHARSTR"lu.sf3dysum",aux,timestep);
   if((fp3 = fopen(sf3d_file, "w"))==NULL) {
      fprintf(stderr,"\nError opening sf3d_file ysum\n");
      exit(1);
   }
   sprintf(sf3d_file, "%s_t%0"TIMECHARSTR"lu.sf3dxsum",aux,timestep);
   if((fp4 = fopen(sf3d_file, "w"))==NULL) {
      fprintf(stderr,"\nError opening sf3d_file xsum\n");
      exit(1);
   }

       /* Open the XDR output file */
     //   sprintf(sf3d_xdrfile, "%s_t%0"TIMECHARSTR"lu.sf3d.xdr",aux,timestep);
     //  if (NULL==(outf=fopen(sf3d_xdrfile,"w"))) {
     //           perror("Could not open XDR output file");
     //           return -1;
     //   }

        /* Set up an XDR stream */
    //    xdrstdio_create(&xdrs,outf,XDR_ENCODE);

        //sprintf(sf3d_xdrfile, "%s_t%0"TIMECHARSTR"lu.sf3dzsum.xdr",aux,timestep);
        //if (NULL==(outf2=fopen(sf3d_xdrfile,"w"))) {
        //        perror("Could not open XDR output file");
        //        return -1;
        //}

        /* Set up an XDR stream */
        //xdrstdio_create(&xdrs2,outf2,XDR_ENCODE);

   
   /* TAYLOR LATTICE SIZE PARAMETERS' TYPE FOR rlft3_d() */
   nx = (unsigned long int) params->nx;
   ny = (unsigned long int) params->ny;
   nz = (unsigned long int) params->nz;
   
	// This is for xfarbe:
   //fprintf(fp2, "%ld %ld\n", nx,ny);
   //fprintf(fp3, "%ld %ld\n", nx,nz);
   //fprintf(fp4, "%ld %ld\n", ny,nz);
	// This is for VTK:
   // fprintf(fp2,"# vtk DataFile Version 2.0\n");
   //fprintf(fp3,"# vtk DataFile Version 2.0\n");
   //fprintf(fp4,"# vtk DataFile Version 2.0\n");
   //fprintf(fp2, "z-sum\n");
   //fprintf(fp3, "y-sum\n");
   //fprintf(fp4, "x-sum\n");
   //fprintf(fp2, "ASCII\n");
   //fprintf(fp3, "ASCII\n");
   //fprintf(fp4, "ASCII\n");
   //fprintf(fp2, "DATASET STRUCTURED_POINTS\n");
   //fprintf(fp3, "DATASET STRUCTURED_POINTS\n");
   //fprintf(fp4, "DATASET STRUCTURED_POINTS\n");
   //fprintf(fp2, "DIMENSIONS %d %d 1\n",nx,ny);
   //fprintf(fp3, "DIMENSIONS %d %d 1\n",nx,nz);
   //fprintf(fp4, "DIMENSIONS %d %d 1\n",ny,nz);
   //fprintf(fp2, "SPACING 1 1 1\n");
   //fprintf(fp3, "SPACING 1 1 1\n");
   //fprintf(fp4, "SPACING 1 1 1\n");
   //fprintf(fp2, "ORIGIN 0 0 0\n");
   //fprintf(fp3, "ORIGIN 0 0 0\n");
   //fprintf(fp4, "ORIGIN 0 0 0\n");
   //fprintf(fp2, "POINT_DATA %d\n",nx*ny);
   //fprintf(fp3, "POINT_DATA %d\n",nx*nz);
   //fprintf(fp4, "POINT_DATA %d\n",nz*ny);
   //fprintf(fp2, "SCALARS scalars float\n");
   //fprintf(fp3, "SCALARS scalars float\n");
   //fprintf(fp4, "SCALARS scalars float\n");
   //fprintf(fp2, "LOOKUP_TABLE default\n");
   //fprintf(fp3, "LOOKUP_TABLE default\n");
   //fprintf(fp4, "LOOKUP_TABLE default\n");
	//
   speq = dmatrix(1, nx, 1, 2 * ny);
   
   avg_field = calc_avg_field(data, nx, ny, nz, stats);
   
   for(i = 1; i <= nx; i++)
      for(j = 1; j <= ny; j++)
	 for(k = 1; k <= nz; k++)  
	 {
	    data[i][j][k] -= (double)avg_field;
	 }
   
   rlft3_d(data, speq, nx, ny, nz, +1);
   
   for(ii = 1; ii <= nx; ii++) {
      if(ii <= ny/2) {
	 i = ii + ny/2;
	 x = -(long int)(nx - i + 1);
      }
      else {
	 i = ii - ny/2;
	 x = (long int)i - 1;
      }
		fprintf(fp2, "\n");
      for(jj = 1; jj <= ny; jj++) {
	 if(jj <= ny/2) {
	    j = jj + ny/2;
	    y = -(long int)(ny - j + 1);
	 }
	 else {
	    j = jj - ny/2;
	    y = (long int)j - 1;
	 }
		//fprintf(fp2, "\n");
	 
	 for(k = 1; k <= nz/2; k++)            /* DONE IN-PLACE!! */
	    data[i][j][k] = (data[i][j][2 * k - 1] *
			     data[i][j][2 * k - 1] + 
			     data[i][j][2 * k] *
			     data[i][j][2 * k]) / 
	       (double)(nx * ny * nz);
	 sumz=0;
	 for(k = nz/2+1; k <= nz; k++) {
	    data[i][j][k] = data[i][j][nz - k + 1];   

	    z = -(long int)(nz - k);
	    //fprintf(fp, "%ld %ld %ld %5.4e\n",
	    //	    x, y, z,
	    //	    data[i][j][k]);
	//	ap = data[i][j][k];
	//	xdr_float(&xdrs,&ap);
		sumz+=(float)data[i][j][k];
	 }

	 for(k = 1; k <= nz/2; k++) {
	    z = (long int)k - 1;
	    //fprintf(fp, "%ld %ld %ld %5.4e\n",
	    //    x, y, z,
	    //   data[i][j][k]);
	//	ap = data[i][j][k];
	//	xdr_float(&xdrs,&ap);
		sumz+=(float)data[i][j][k];
	 }
	//	xdr_float(&xdrs2,&sumz);
	//	xfarbe/gnuplot:
	  fprintf(fp2, "%ld %ld %5.4e\n",x, y, sumz);
	//	    VTK:
	    // fprintf(fp2, "%5.4e\n", sumz);
      }
   }

/*
	IMPORTANT: the k-loop reaching up to nz-1 plus the
	de-wrapping around on z being different than that
	for x and y agrees with the Nyquist critical freq.
	on z, f_3=-f_c, not being in data[][][] but in
	speq[][] instead. We ignore it.
*/
/*
	    fprintf(fp, "%ld %ld %ld %5.4e\n", 
	            x, y, z,
	            data[i][j][k]);
*/

   // Sum up in y direction.
	 for(k = nz/2+1; k <= nz; k++) {
		 fprintf(fp3,"\n");
   	  for(i = nx/2; i <= nx; i++) {
		 sumy = 0;
	 // ???? xi = i - ny/2;
	 //xi = (long int)xi - 1;
	    for(j = 1; j <= ny; j++) {
	         sumy+=(float)data[i][j][k];
	    }
	// VTK:	 
	//       fprintf(fp3, "%5.4e\n", sumy);
	//	xfarbe/gnuplot:
	  fprintf(fp3, "%ld %ld %5.4e\n",-(long int)(nx - i), -(long int)(nz - k +1), sumy);
	  }
	 
	  for(i=1;i<nx/2;i++) {
	   sumy = 0;

	    for(j = 1; j <= ny; j++) {
	        sumy+=(float)data[i][j][k];
	    }
	  //  fprintf(fp3, "%5.4e\n", sumy);
	  fprintf(fp3, "%ld %ld %5.4e\n",i, -(long int)(nz - k + 1), sumy);
	  }
         }
	   
	 for(k = 1; k <= nz/2; k++){
		 fprintf(fp3,"\n");
   	  for(i = nx/2; i <= nx; i++) {
		 sumy = 0;
	    for(j = 1; j <= ny; j++) {
	         sumy+=(float)data[i][j][k];
	    }
	      // fprintf(fp3, "%5.4e\n", sumy);
	       fprintf(fp3, "%ld %ld %5.4e\n",-(long int)(nx - i), (long int)(k-1), sumy);
	  }
	  for(i=1;i<nx/2;i++) {
	   sumy = 0;
	    for(j = 1; j <= ny; j++) {
	        sumy+=(float)data[i][j][k];
	    }
	   // fprintf(fp3, "%5.4e\n", sumy);
	    fprintf(fp3, "%ld %ld %5.4e\n",i, (long int)(k-1), sumy);
	  }
        }


   // Sum up in x direction.
   
   	for(j = ny/2; j <= ny; j++) {
		 fprintf(fp4,"\n");
	 for(k = nz/2+1; k <= nz; k++) {
	       sumx = 0;
	    for(i = 1; i <= nx; i++) {
	       sumx+=(float)data[i][j][k];
	    }
	       //fprintf(fp4, "%5.4e\n", sumx);
	       fprintf(fp4, "%ld %ld %5.4e\n",(long int)(j-ny), (long int)(k-nz-1), sumx);
	  }
	 
	  for(k = 1; k <= nz/2; k++){
	        sumx = 0;
	    for(i = 1; i <= nx; i++) {
	        sumx+=(float)data[i][j][k];
	    }
	    //fprintf(fp4, "%5.4e\n", sumx);
	       fprintf(fp4, "%ld %ld %5.4e\n",(long int)(j-ny), (long int)(k-1), sumx);
	  }
         }

   	for(j = 1; j < ny/2; j++) {
		 fprintf(fp4,"\n");
	 for(k = nz/2+1; k <= nz; k++) {
	       sumx = 0;
	    for(i = 1; i <= nx; i++) {
	       sumx+=(float)data[i][j][k];
	    }
	       //fprintf(fp4, "%5.4e\n", sumx);
	       fprintf(fp4, "%ld %ld %5.4e\n",(long int)(j-1), (long int)(k-nz-1), sumx);
	  }
	 
	  for(k = 1; k <= nz/2; k++){
	        sumx = 0;
	    for(i = 1; i <= nx; i++) {
	        sumx+=(float)data[i][j][k];
	    }
	    //fprintf(fp4, "%5.4e\n", sumx);
	       fprintf(fp4, "%ld %ld %5.4e\n",(long int)(j-1), (long int)(k-1), sumx);
	  }
         }
	   
	   

   //fflush(fp);
   fflush(fp2);
   fflush(fp3);
   fflush(fp4);
   
   printf("\nWrote fft to disc...");
   //fclose(fp);
   fclose(fp2);
   fclose(fp3);
   fclose(fp4);
   //xdr_destroy(&xdrs);
   
   free_dmatrix(speq, 1, nx, 1, 2 * ny);
   
   counter++;
   
   return(1);
}


