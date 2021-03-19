 /* Read XDR from the specified file, write ASCII to stdout.
 *
 * Syntax: xdr2all nx ny nz ox oy oz type filename 
 *
 * where type is one of scalar, 2scalar, 3scalar, vector
 *
 */

#include<stdio.h>
#include<float.h>
#include<rpc/rpc.h>
#include<sys/types.h>
#include<sys/stat.h>
#include "main.h"

extern off_t   file_size(char *);




int decode(Params *params, unsigned long timestep, Dat *dat,
	   Stats *stats)
{                      
   int     i, j, k;
   /* NB: type "double" should be at least 64 bits wide.
    * if not, it's trivial to change this code to use long double
    * for the buffers.
    */
   double       a, b;
   char         xdrfile[160];
   FILE	        *infile;
   XDR	        xdrs;
   
   /* xdrfile = (char *)malloc(PATHNAMELEN * sizeof(char)); */
   
  
   sprintf(xdrfile, "%s_t%0"TIMECHARSTR"lu.xdr", params->datafile, timestep);
   if(file_size(xdrfile) == 0) {
      fprintf(stderr, "\n\nError: File %s has zero length\n", xdrfile);
      fflush(stderr);
      free(xdrfile);
      /*  return(1) is decode() OK, 
       *  return(0) is zero filelength, 
       *  return(-1) is error */
      return(0);
   }


   /* Open the data file */

   if (NULL == (infile = fopen(xdrfile,"r"))) {
      fprintf(stderr, "\nError: Could not open data file %s\n",
	      xdrfile);
      free(xdrfile);
      return(-1);
   }
   printf("\nReading file %s\n", xdrfile);
   
   /* Set up an XDR stream */
   
   xdrstdio_create(&xdrs, infile, XDR_DECODE);
   
   
   /* Now read and parse stdin.
      The order of loops matches with how
      ***data is malloc'ed. */
   
   stats->min = stats-> min2 = DBL_MAX;
   stats->max = stats->max2 = -DBL_MAX;

   /* fp = fopen("data.dat", "w"); */

   for(k = 1; k <= params->nz; k++)
      for(j = 1; j <= params->ny;j++)
	 for(i = 1; i <= params->nx; i++)
	 {
	    xdr_double(&xdrs, &a);
	    dat->data[i][j][k] = (double) a;   
	    if(a > stats->max)  stats->max = a;
	    if(a < stats->min)  stats->min = a;
	    

	    if(params->owd) {
	       xdr_double(&xdrs, &b);
	       dat->data2[i][j][k] = (double) b;
	       if(b > stats->max2)  stats->max2 = b;
	       if(b < stats->min2)  stats->min2 = b;
	    }


            /* 
	       IMPORTANT: Would we lose
	       precision if we cast it to float? 
	    */

            /*	    
	    fprintf(fp, "%d %d %d %8.6e\n", i, j, k, data[i][j][k]);
 	    fflush(fp); 
	    */
	 }
  
   /* fclose(fp); */


   /* Close the XDR stream */

   xdr_destroy(&xdrs);


   fclose(infile);
/*   free(xdrfile); */
   return(1);   /* And they all lived happily ever after. */

}






