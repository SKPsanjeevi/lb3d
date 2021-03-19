 /* Read HDF from the specified file, write ASCII to stdout.
 *
 * Syntax: xdr2all nx ny nz ox oy oz type filename 
 *
 * where type is one of scalar, 2scalar, 3scalar, vector
 *
 */

#include<float.h>
#include<rpc/rpc.h>
#include<sys/types.h>
#include<sys/stat.h>
#include "main.h"

extern off_t   file_size(char *);
extern int hdfread(Params *,Dat *,Stats *,char *);

int decode(Params *params, unsigned long timestep, Dat *dat,
	   Stats *stats)
{                      
   int     i, j, k;
   /* NB: type "double" should be at least 64 bits wide.
    * if not, it's trivial to change this code to use long double
    * for the buffers.
    */
   float        a,b;
   char         hdffile[160];
   FILE	        *infile;

   sprintf(hdffile, "%s_t%0"TIMECHARSTR"lu.h5", params->datafile, timestep);

   if(file_size(hdffile) == 0) {
      fprintf(stderr, "\n\nError: File %s has zero length\n", hdffile);
      fflush(stderr);
      free(hdffile);
      /*  return(1) is decode() OK, 
       *  return(0) is zero filelength, 
       *  return(-1) is error */
      return(0);
   }


   /* Open the data file */

   /*   if (NULL == (infile = fopen(xdrfile,"r"))) {
      fprintf(stderr, "\nError: Could not open data file %s\n",
	      xdrfile);
      free(xdrfile);
      return(-1);
      } 
      printf("\nReading file %s\n", xdrfile); */
   

   /* fp = fopen("data.dat", "w"); */

     /* HDF part */
     hdfread(params,dat,stats,hdffile);

/*   free(xdrfile); */
   return(1);   /* And they all lived happily ever after. */

}






