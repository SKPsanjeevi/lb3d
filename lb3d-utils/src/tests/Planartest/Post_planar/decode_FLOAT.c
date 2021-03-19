 /* Read XDR from the specified file, write ASCII to stdout.
 *
 * Syntax: xdr2all nx ny nz ox oy oz type filename 
 *
 * where type is one of scalar, 2scalar, 3scalar, vector
 *
 */

#include<stdio.h>
#include<rpc/rpc.h>
#include<sys/types.h>
#include<sys/stat.h>
#include <hdf5.h>

#include "main.h"

extern off_t   file_size(char *);
extern int     extrema(Stats *, float *);







int decode(char *prefix, Params *params, unsigned long timestep,
	   float ***data, Stats *stats)
{                      
   int     i, j, k;
   /* NB: type "double" should be at least 64 bits wide.
    * if not, it's trivial to change this code to use long double
    * for the buffers.
    */
   float   a=0.0;
   char    *xdrfile;
   FILE	   *infile;
   FILE *READPIPE;
   XDR	   xdrs;

   char	   *xdrtestfile;
   char    *filetest;
   char	    *alt_filename;
   
   xdrfile = (char *)malloc(PATHNAMELEN * sizeof(char));
   xdrtestfile = (char *)malloc(PATHNAMELEN * sizeof(char));
   alt_filename = (char *)malloc(PATHNAMELEN * sizeof(char));
   filetest = (char *)malloc((PATHNAMELEN +20) * sizeof(char));

   /* Open the data file */
   
   sprintf(xdrfile, "%s/%s_%s_t%06lu.xdr",  params->datapath, prefix, params->gr_out_file, timestep);

   sprintf(xdrtestfile, "%s/%s_%s_t%06lu-*.xdr", params->datapath, prefix, params->gr_out_file, timestep);	
   sprintf(filetest, "ls %s", xdrtestfile);	


   if(file_size(xdrfile) == 0) {
	   if( NULL==(READPIPE=popen(filetest, "r"))){
		   fprintf(stderr, "Cannot open pipe with command ' %s '!\n", filetest);
		   return(-1);
	   }
	   else{
			fscanf(READPIPE, "%s",alt_filename );
			if(file_size(alt_filename) == 0) {
				fprintf(stderr, "\n\tError: File %s has zero length\n", alt_filename);
				fflush(stderr);
				free(xdrfile);
				/*  return(1) is decode() OK, 
				*  return(0) is zero filelength, 
				*  return(-1) is error */
				return(0);
				}
			else{
				sprintf(xdrfile, "%s", alt_filename);
			}	
	   }
   }
   
   if (NULL == (infile = fopen(xdrfile,"r"))) {
      fprintf(stderr, "\n\tError: Could not open data file %s\n",   xdrfile);
      fprintf(stderr, "\a\nError in reading size of file %s\n"
		      "\tIn function tools.c:file_size() called by decode()\n",
		      xdrfile);
      fflush(stderr);
      free(xdrfile);
      return(-1);
   }
   printf("\n\tReading file %s", xdrfile);
   
   /* Set up an XDR stream */
   
   xdrstdio_create(&xdrs, infile, XDR_DECODE);
   
   
   /* Now read and parse stdin.
      The order of loops matches 
      with how data was dumped onto 
      disk by fortran. */
   
   /*
   fp = fopen("data_xyz.dat", "w");
   */

   if(VERBOSE)
      printf("\ndecode_FLOAT.c: (JUST BEFORE READ DATA)"
	     " params->nx=%d params->ny=%d",
	     params->nx, params->ny);
   for(k = 1; k <= params->nz; k++) {
      for(j = 1; j <= params->ny; j++) {
	 for(i = 1; i <= params->nx; i++) {
	    xdr_float(&xdrs, &a);
	    data[i][j][k] =  a;    
	    extrema(stats, &a);
	 }		 
      }
   }
   /*
    *  This following part added on January 2001
    *  to deal with other surface tension mechanical
    *  definition.
    */
   /*
   if(strcmp(prefix,"pzz")==0) {
   */
     /*  Only for PERPENDICULAR component.
      *  Choose the value away from the interface,
      *  which is the isotropic pressure and we 
      *  assume it's at k==params->nz/4.
      *  There's symmetry on x and y so we choose
      *  a middle point on this plane.
      *  We also check that this value is the maximum
      *  of the profile.
      */
   /*
     value = data[params->nx/2][params->ny/2][params->nz/4];
     printf("\nMax Pzz value is %16.8f"
	    "\nPzz at (%d,%d,%d) is %16.8f\n", 
	    stats->max,
	    params->nx/2, params->ny/2, params->nz/4, 
	    value);

     for(k = 1; k <= params->nz; k++) {
	for(j = 1; j <= params->ny; j++) {
	for(i = 1; i <= params->nx; i++) {		  
	   data[i][j][k] = value;
	   }
	}
     }
   }
   */

   /*
     fclose(fp);
     
     fp = fopen("data_z.dat", "w");
     for(k = 1; k <= params->nz; k++) {
     fprintf(fp, "\n%d %16.8f", k, data[params->nx / 2][params->ny / 2][k]);
     fflush(fp); 
     }
   */
   
   /*
     printf("%16.8f < field < %16.8f\n", 
     stats->min, stats->max);  
     fflush(stdout);
   */
	  
   /* Close the XDR stream */
      
   xdr_destroy(&xdrs);
      
   /*
     fclose(fp);
   */
   fclose(infile);
   free(xdrfile);
   return(1);   /* And they all lived happily ever after. */
   
}








