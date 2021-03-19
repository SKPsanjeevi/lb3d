#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"


extern int  dump_second_field, dump_third_field;

int
dump_to_files(CASF *casf, Params *params, Stats *stats, unsigned long
	      timestep)
{
   unsigned long m;
   double        k;
   char          aux[128] = "\0";
   char          sf_file[128] = "\0";
   char          siz_file[128] = "\0";
   char          lsz_file[128] = "\0";
   FILE          *sf_fp, *siz_fp, *lsz_fp;

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
   printf("\n\nFIELD %s -----------------"
     "\nOutput files' root is: %s"
     "\nMax value is %g"
     "\nMin value is %g"
     "\nAvg value is %g\n",
     params->field,
     params->outfile,
     stats->max,
     stats->min,
     stats->avg);

   fflush(stdout);

   sprintf(sf_file, "%s_t%0"TIMECHARSTR"lu.sf", params->outfile, timestep);
   sprintf(siz_file, "%s.siz", params->outfile);
   sprintf(lsz_file, "%s.lsz", params->outfile);
   

   /* DELETE FILES IF ALREADY EXISTING...*/
   if(remove(sf_file) == 0) 
      printf("\nAlready existing %s removed.", sf_file);
   if(timestep == 0) {
      if(remove(siz_file) == 0) 
	 printf("\nAlready existing %s removed.", siz_file);
      if(remove(lsz_file) == 0) 
	 printf("\nAlready existing %s removed.", lsz_file);
   }
   printf("\n");

   /* OPEN FILES */
   if(NULL==(sf_fp = fopen(sf_file, "a"))) {
      fprintf(stderr,"\nCould not open sf_file: %s", sf_file); 
      perror("\nCould not open output sf_file");
      return(-1);
   }
   printf("\nWriting onto file %s", sf_file);

   if(NULL==(siz_fp = fopen(siz_file, "a"))) {
      fprintf(stderr,"\nCould not open siz_file: %s", siz_file); 
      perror("\nCould not open output siz_file");
      return(-1);
   }
   printf("\nWriting onto file %s", siz_file);

   if(NULL == (lsz_fp = fopen(lsz_file, "a"))) {
      fprintf(stderr, "\nCould not open lsz_file: %s", lsz_file); 
      perror("\nCould not open output lsz_file");
      return(-1);
   }
   printf("\nWriting onto file %s", lsz_file);

   printf("\n");
   

   for(m = 0 ; m < params->kradius ; m++) {
      /* k = DPI / L * m; */
      k = DPI / (double)params->nx * m;   /* IMPORTANT: This is suited
					    for square lattices only.

					    IMPORTANT: We assume the
					    system size is the number
					    of nodes it has.
					    */

/*    DEBUG */
      fprintf(sf_fp, 
	      "%15.10e %15.10e %15.10e\n"
/*	      "%15.10e %15.10e\n"  */
	      ,k, casf->S[m], casf->errS[m]
/*	      ,(casf->sum_sqS[m] / (double)casf->sum_shell[m]) */
/*	      ,(casf->S[m] * casf->S[m]) */

	      ); 
      fflush(sf_fp); 
   }
/*    END DEBUG */


   fprintf(siz_fp, 
	   "%lu "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e %15.10e\n", 
	   timestep, 
	   casf->k1, casf->errk1,
	   casf->k2, casf->errk2,
	   casf->R1, casf->errR1,
	   casf->R2, casf->errR2,
	   stats->min, stats->max, stats->avg); 
   fflush(siz_fp); 

   if(timestep > 0) {
     fprintf(lsz_fp, 
	     "%15.10e %15.10e %15.10e\n", 
	     log((double)timestep), 
	     log(casf->R1),
	     casf->errR1 / casf->R1); 
     fflush(lsz_fp); 
   }

   fclose(sf_fp);
   fclose(siz_fp);
   fclose(lsz_fp);
   return(1);

}







