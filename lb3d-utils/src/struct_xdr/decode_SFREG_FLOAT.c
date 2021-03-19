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




int decode(Params *params, Dat *dat, Stats *stats)
{                      
   int          i, j, k;
   float       a;
   char         xdrfile[160];
   FILE	        *infile;
   XDR	        xdrs;
   
  
   sprintf(xdrfile, "%s", params->inpfile);
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
   stats->min = DBL_MAX;
   stats->max = -DBL_MAX;

   /* fp = fopen("data.dat", "w"); */

   for(k = 1; k <= params->nz; k++)
      for(j = 1; j <= params->ny;j++)
	 for(i = 1; i <= params->nx; i++)
	 {

/* To read DBL precision data just uncomment next */
/* line and comment out the one that follows it:  */

/*	    xdr_double(&xdrs, &a); */
	    xdr_float(&xdrs, &a);
	    dat->data[i][j][k] = (double)a;   
	    if(a > stats->max)  stats->max = (double)a;
	    if(a < stats->min)  stats->min = (double)a;   
	 }


   /* Close the XDR stream */
   xdr_destroy(&xdrs);


   fclose(infile);
   return(1);

}






