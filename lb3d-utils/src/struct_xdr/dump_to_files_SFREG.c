#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"


extern int  dump_second_field;

int
dump_to_files(CASF *casf, Params *params, Stats *stats)
{
   unsigned long m;
   double        k;
   char          sf_file[128] = "\0";
   FILE          *sf_fp;


   /* CREATE FILE NAMES */
   /*
     fprintf(stderr, "params->filename in dump_to_files(): \n\t%s\n", 
     params->filename);
     fprintf(stderr, "sf_file in dump_to_files() before sprintf(): \n\t%s\n", 
     sf_file);
   */
   

   /*
     PRINT SOME INFO..
   */
   printf("\n\noutfile is: %s"
	  "\nMax value is %g"
	  "\nMin value is %g"
	  "\nAvg value is %g\n",
	  params->outfile,
	  stats->max,
	  stats->min,
	  stats->avg);

   fflush(stdout);


   strcpy(sf_file, params->outfile);

   /* DELETE FILE IF ALREADY EXISTING...*/
   if(remove(sf_file) == 0) 
      printf("\nAlready existing %s removed.", sf_file);
   printf("\n");


   /* OPEN FILE */
   if(NULL==(sf_fp = fopen(sf_file, "a"))) {
      fprintf(stderr,"\nCould not open sf_file: %s", sf_file); 
      perror("\nCould not open output sf_file");
      return(-1);
   }
   printf("\nWriting onto file %s", sf_file);

   printf("\n");
   

   for(m = 0 ; m < params->kradius ; m++) {
      /* k = DPI / L * m; */
      k = DPI / (double)params->nx * m;   /* IMPORTANT: This is suited
					    for square lattices only.

					    IMPORTANT: We assume the
					    system size is the number
					    of nodes it has.
					    */
      fprintf(sf_fp, 
	      "%15.10e %15.10e %15.10e\n"
	      , k, casf->S[m], casf->errS[m]
	 ); 
      fflush(sf_fp); 
   }
   
   fclose(sf_fp);
   return(1);
   
}







