#include <stdio.h>
#include "main.h"
#include "nrutil.h"

extern int    read_data(Params *);
extern int    decode(Params *, unsigned long, Dat *, Stats *);
extern int    fft(double ***, Params *, Stats *);
extern int    circ_avg(double ***, Params *, CASF *);
extern int    domain_size(CASF *, Params *, Stats *);
extern int    dump_to_files(CASF *, Params *, Stats *, unsigned long);
extern int    *i1t(int n1);
extern double *d1t(int n1);
extern double ***d3tensor(long, long, long, long, long, long);
extern void   free_d3tensor(double ***t, long, long, long, long, long,
			    long);


int           dump_second_field, dump_third_field;

int
caller(Params *params)
{
   unsigned long ts;
   Dat           *dat;
   Stats         *stats;
   CASF          *casf;
   int           dec,iter;
   char		 first_field[20], second_field[20], third_field[20];  
   char          outfile[160], outfile2[160], outfile3[160];

   if(NULL==(dat=(Dat *)malloc(sizeof(Dat)))) {
      fprintf(stderr,"\nError in allocating mem for Dat struct.\n");
      exit(-1);
   }
   if(NULL==(stats=(Stats *)malloc(sizeof(Stats)))) {
      fprintf(stderr,"\nError in allocating mem for Stats struct.\n");
      exit(-1);
   }
   if(NULL==(casf=(CASF *)malloc(sizeof(CASF)))) {
      fprintf(stderr,"\nError in allocating mem for CASF struct.\n");
      exit(-1);
   }

   if(read_data(params) != 1) {
      fprintf(stderr,"Error in read_data() function\n");
      exit(-1);
   }

   dat->data = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   if(params->owd)
      dat->data2 = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   if(params->arr) {
      dat->data2 = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);
      dat->data3 = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   }
   if((params->nx <= params->ny) && (params->nx <= params->nz)) 
      params->kradius = params->nx / 2;
   else if((params->ny <= params->nx) && (params->ny <= params->nz)) 
      params->kradius = params->ny / 2;
   else if((params->nz <= params->nx) && (params->nz <= params->ny)) 
      params->kradius = params->nz / 2;

   casf->S         = d1t(params->kradius);
   casf->errS      = d1t(params->kradius);
   casf->sum_sqS   = d1t(params->kradius);
   casf->sum_shell = (unsigned long *)i1t(params->kradius);


   /*
     CREATE DATA FILE NAME
   */
   sprintf(params->datafile, "%s/%s/%s_%s", params->datapath,
	   params->folder, params->filetype, params->studyname);
   /*  IMPORTANT: sprintf() 
       can't concatenate two strings
       and store the result in one of them. 
   */


   /*
     CREATE OUTPUT FILE NAME
   */
   /*
     I'm writing all taskfarming 
     SF files onto the same 
     directory, as it's simpler..

     if(params->ntask > 1) {
        make_dir(params, folder);
	sprintf(params->outfile, "%s/%s/%s_%s", params->outpath,
	folder, params->filetype, gr_out_file);
     }
     else 
   */
   
   /* ENABLE -res FLAG ONLY FOR NON-TASKFARMING...
    */
   if(params->owd) {
      strcpy(first_field, "od");
      strcpy(second_field, "wd");
      sprintf(outfile, "%s/%s_%s", params->outpath,
	      first_field, params->studyname);
      sprintf(outfile2, "%s/%s_%s", params->outpath,
	      second_field, params->studyname); 
   } 
   else if(params->arr) {
      strcpy(first_field, "arrx");
      strcpy(second_field, "arry");
      strcpy(third_field, "arrz");
      sprintf(outfile, "%s/%s_%s", params->outpath,
              first_field, params->studyname);
      sprintf(outfile2, "%s/%s_%s", params->outpath,
              second_field, params->studyname);
      sprintf(outfile3, "%s/%s_%s", params->outpath,
              third_field, params->studyname);
   } 
   else {
      strcpy(first_field, params->filetype);
      sprintf(outfile, "%s/%s_%s", params->outpath,
	      first_field, params->studyname);
   }


   /*
     TIME LOOP...
   */
   params->sci_start=0;
   for(ts = params->sci_start; 
       ts <= params->n_iteration; 
       ts += params->n_sci) {

      printf("\n\n\nTimestep %lu\t=========================\n", ts);
      do {
	 dec = decode(params, ts, dat, stats);
	 if(dec != 1) {
	    fprintf(stderr, "\n\tSkipping timestep %lu\n", ts);
	    ts += params->n_sci;
	    if(ts > params->n_iteration) {
	       printf("\nNO TIME STEP READ!\n");
	       break;
	    }
	 }
      } while(dec != 1);


      if(ts > params->n_iteration) {
	 break;
      }

//    There are obviously more succint (and decent) ways of calling the 
//    sequence: {fft(); circ_avg(); domain_size(); dump_to_files();}
//    more than once if there is more than one field per register
//    in the .xdr datafile to apply it to!  I can't spare more 
//    time now...	--NGS, Oct.03


      dump_second_field = 0;      
      dump_third_field = 0;
      if(fft(dat->data, params, stats) != 1) {
	 fprintf(stderr,"Error in fft() function\n");
	 exit(-1);
      } 
      if(circ_avg(dat->data, params, casf) != 1) {
	 fprintf(stderr,"Error in circ_avg() function\n");
	 exit(-1);
      } 
      if(domain_size(casf, params, stats) != 1) {
	 fprintf(stderr,"Error in domain_size() function\n");
	 exit(-1);
      }
      strcpy(params->field,first_field);
      strcpy(params->outfile,outfile);
      if(dump_to_files(casf, params, stats, ts) != 1) {
	 fprintf(stderr,"Error in domain_size() function\n");
	 exit(-1);
      }
//
//    DUMP SECOND FIELD IF OWD
//     
      if(params->owd) {
         dump_second_field = 1;
	 dump_third_field  = 0;
	 if(fft(dat->data2, params, stats) != 1) {
	    fprintf(stderr,"Error in fft() function\n");
	    exit(-1);
	 } 
	 if(circ_avg(dat->data2, params, casf) != 1) {
	    fprintf(stderr,"Error in circ_avg() function\n");
	    exit(-1);
	 } 
	 if(domain_size(casf, params, stats) != 1) {
	    fprintf(stderr,"Error in domain_size() function\n");
	    exit(-1);
	 } 
	 strcpy(params->field,second_field);
	 strcpy(params->outfile,outfile2);
	 stats->max = stats->max2;
         stats->min = stats->min2;
         stats->avg = stats->avg2;
	 if(dump_to_files(casf, params, stats, ts) != 1) {
	    fprintf(stderr,"Error in domain_size() function\n");
	    exit(-1);
	 } 
      }
//
//    DUMP SECOND AND THIRD FIELDS IF ARR
//      
      else if(params->arr) {
         dump_second_field = 1;
	 dump_third_field  = 0;
         if(fft(dat->data2, params, stats) != 1) {
            fprintf(stderr,"Error in fft() function\n");
            exit(-1);
         }
         if(circ_avg(dat->data2, params, casf) != 1) {
            fprintf(stderr,"Error in circ_avg() function\n");
            exit(-1);
         }
         if(domain_size(casf, params, stats) != 1) {
            fprintf(stderr,"Error in domain_size() function\n");
            exit(-1);
         }
         strcpy(params->field,second_field);
         strcpy(params->outfile,"\0"); strcpy(params->outfile,outfile2);
         stats->max = stats->max2;
         stats->min = stats->min2;
         stats->avg = stats->avg2;
         if(dump_to_files(casf, params, stats, ts) != 1) {
            fprintf(stderr,"Error in domain_size() function\n");
            exit(-1);
         }
//       
	 dump_second_field = 0;
         dump_third_field  = 1;
         if(fft(dat->data3, params, stats) != 1) {
            fprintf(stderr,"Error in fft() function\n");
            exit(-1);
         }
         if(circ_avg(dat->data3, params, casf) != 1) {
            fprintf(stderr,"Error in circ_avg() function\n");
            exit(-1);
         }
         if(domain_size(casf, params, stats) != 1) {
            fprintf(stderr,"Error in domain_size() function\n");
            exit(-1);
         }
         strcpy(params->field,third_field);
         strcpy(params->outfile,"\0"); strcpy(params->outfile,outfile3);
         stats->max = stats->max3;
         stats->min = stats->min3;
         stats->avg = stats->avg3;
         if(dump_to_files(casf, params, stats, ts) != 1) {
            fprintf(stderr,"Error in domain_size() function\n");
            exit(-1);
         }
      }
   }

   fflush(stderr);
   free(stats); 
   free(casf);
   free_d3tensor(dat->data, 1, params->nx, 1, params->ny, 1,
		 params->nz);  
   if(params->owd)
      free_d3tensor(dat->data2, 1, params->nx, 1, params->ny, 1,
		    params->nz); 
   if(params->arr) {
      free_d3tensor(dat->data2, 1, params->nx, 1, params->ny, 1,
                    params->nz);
      free_d3tensor(dat->data3, 1, params->nx, 1, params->ny, 1,
                    params->nz);
   }
   free(dat);
   free(casf->S);   /* CAREFUL WITH THIS WHEN CALLING caller() MORE
		       THAN ONCE.   WHY???*/
   free(casf->sum_shell);

   fprintf(stdout, "\ncaller() done!\n"); 
   return(1);

}


