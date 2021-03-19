#include "main.h"


extern int  dump_second_field;

int
dump_to_files(CASF *casf, XYAVGSF *xyavgsf, Params *params, Stats *stats, 
	      unsigned long timestep)
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
   if(dump_second_field) {
      strcpy(aux, params->outfile2);
      printf("\n\nFIELD %s -----------------"
	     "\noutfile is: %s"
	     "\nMax value is %g"
	     "\nMin value is %g"
	     "\nAvg value is %g\n",
	     "wd",
	     aux,
	     stats->max2,
	     stats->min2,
	     stats->avg2);
   }
   else {
      strcpy(aux, params->outfile);
      printf("\n\nFIELD %s -----------------"
	     "\noutfile is: %s"
	     "\nMax         = %g"
	     "\nMin         = %g"
	     "\nAvg         = %g"
	     "\ncasf->R1    = %g +- %g"
	     "\nxyavgsf->Rx = %g +- %g"
	     "\nxyavgsf->Ry = %g +- %g"
	     "\nxyavgsf->Rz = %g +- %g\n",
	     params->first_field,
	     aux,
	     stats->max, stats->min, stats->avg,
	     casf->R1, casf->errR1,
	     xyavgsf->Rx, xyavgsf->Rxerr,
	     xyavgsf->Ry, xyavgsf->Ryerr,
	     xyavgsf->Rz, xyavgsf->Rzerr);
   }
   fflush(stdout);


   sprintf(sf_file, "%s_t%0"TIMECHARSTR"lu.sf", aux, timestep);
   sprintf(siz_file, "%s.siz", aux);
   sprintf(lsz_file, "%s.lsz", aux);
   

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

// I don't think we need this log data --NGS
/*
   if(NULL == (lsz_fp = fopen(lsz_file, "a"))) {
      fprintf(stderr, "\nCould not open lsz_file: %s", lsz_file); 
      perror("\nCould not open output lsz_file");
      return(-1);
   }
   printf("\nWriting onto file %s", lsz_file);
*/

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
	      ,k, casf->S[m], casf->errS[m]
	      ); 
      fflush(sf_fp); 
   }

   fprintf(siz_fp, 
	   "%lu "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e "
	   "%15.10e %15.10e "
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
	   xyavgsf->Rx, xyavgsf->Rxerr,
	   xyavgsf->Ry, xyavgsf->Ryerr,
	   xyavgsf->Rz, xyavgsf->Rzerr,
	   stats->min, stats->max, stats->avg); 
   fflush(siz_fp); 

// I don't think we need this log data --NGS
/*
   if(timestep > 0) {
     fprintf(lsz_fp, 
	     "%15.10e %15.10e %15.10e\n", 
	     log((double)timestep), 
	     log(casf->R1),
	     casf->errR1 / casf->R1); 
     fflush(lsz_fp); 
   }
*/

   fclose(sf_fp);
   fclose(siz_fp);
//   fclose(lsz_fp);
   return(1);

}







