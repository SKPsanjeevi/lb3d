#include "main.h"
#include "nrutil.h"

extern int    read_data(Params *);
extern int    decode(Params *, unsigned long, Dat *, Stats *);
extern int    fft(double ***, Params *, Stats *);
extern int    circ_avg(double ***, Params *, CASF *);
extern int    domain_size(CASF *, Params *, Stats *);
extern int    xydomain_size(double ***, Params *, XYAVGSF *);
extern int    dump_to_files(CASF *, XYAVGSF *, Params *, Stats *, unsigned long);
extern unsigned long    *uli1t(unsigned long n1);
extern double *d1t(int n1);
extern double ***d3tensor(long, long, long, long, long, long);
extern void   free_d3tensor(double ***t, long, long, long, long, long,
			    long);


int           dump_second_field;

int
caller(Params *params)
{
   unsigned long ts;
   Dat           *dat;
   Stats         *stats;
   CASF          *casf;
   XYAVGSF       *xyavgsf;
   int           dec;
   




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
   if(NULL==(xyavgsf=(XYAVGSF *)malloc(sizeof(XYAVGSF)))) {
      fprintf(stderr,"\nError in allocating mem for XYAVGSF struct.\n");
      exit(-1);
   }


   /* printf("\nIn caller(params), just before calling read_data(params) \n" */
/* 	  "params->nx = %lu\n" */
/* 	  "params->ny = %lu\n" */
/* 	  "params->nz = %lu\n",  */
/* 	  params->nx, params->ny, params->nz); */
   
   if(read_data(params) != 1) {
      fprintf(stderr,"Error in read_data() function\n");
      exit(-1);
   }

/*    printf("\nIn caller(params), just after returning from read_data(params) \n" */
/* 	  "params->nx = %lu\n" */
/* 	  "params->ny = %lu\n" */
/* 	  "params->nz = %lu\n",  */
/* 	  params->nx, params->ny, params->nz);   */
   

   dat->data = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   if(params->owd)
      dat->data2 = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);

   if((params->nx <= params->ny) && (params->nx <= params->nz)) 
      params->kradius = params->nx / 2;
   else if((params->ny <= params->nx) && (params->ny <= params->nz)) 
      params->kradius = params->ny / 2;
   else if((params->nz <= params->nx) && (params->nz <= params->ny)) 
      params->kradius = params->nz / 2;

   casf->S         = d1t(params->kradius);
   casf->errS      = d1t(params->kradius);
   casf->sum_sqS   = d1t(params->kradius);
   casf->sum_shell = (unsigned long *)uli1t(params->kradius);


// Healthy measure:
   xyavgsf->xxS = 
      xyavgsf->yyS =
      xyavgsf->zzS =
      xyavgsf->xxSerr =
      xyavgsf->yySerr =
      xyavgsf->zzSerr =
      xyavgsf->Rx =
      xyavgsf->Ry =
      xyavgsf->Rz =
      xyavgsf->Rxerr =
      xyavgsf->Ryerr = 
      xyavgsf->Rzerr = 0.0;

   /*
     CREATE DATA FILE NAME
   */
   sprintf(params->datafile, "%s/%s/%s_%s", params->datapath,
	   params->folder, params->filetype, params->studyname);
   /*  IMPORTANT: sprintf() 
       can't concatenate two strings
       and store the result in one of them. 
   */


   if(params->owd) {
      sprintf(params->outfile, "%s/%s_%s", params->outpath,
	      "od", params->studyname);
      sprintf(params->outfile2, "%s/%s_%s", params->outpath,
	      "wd", params->studyname); 
      strcpy(params->first_field, "od");
   }
   else {
      sprintf(params->outfile, "%s/%s_%s", params->outpath,
	      params->filetype, params->studyname);
      strcpy(params->first_field, params->filetype);
   }




   /*
     TIME LOOP...
   */
 //  params->sci_start=0;
   for(ts = params->sci_start; 
       ts <= params->n_iteration; 
       ts += params->n_sci) 
   {
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

      dump_second_field = 0;      
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

      if(xydomain_size(dat->data, params, xyavgsf) != 1) {
	 fprintf(stderr,"Error in domain_size() function\n");
	 exit(-1);
      }

      if(dump_to_files(casf, xyavgsf, params, stats, ts) != 1) {
	 fprintf(stderr,"Error in domain_size() function\n");
	 exit(-1);
      }


      /*
	DUMP SECOND FIELD IF OWD...
      */
      dump_second_field = 1;
      if(params->owd) {
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
	 if(dump_to_files(casf, xyavgsf, params, stats, ts) != 1) {
	    fprintf(stderr,"Error in domain_size() function\n");
	    exit(-1);
	 } 
      }
   }

   fflush(stderr);
   free(stats); 

   free(casf->S);   /* CAREFUL WITH THIS WHEN CALLING caller() MORE
		       THAN ONCE.   WHY???*/
   free(casf->errS);
   free(casf->sum_sqS);
   free(casf->sum_shell);
   free(casf);

   free_d3tensor(dat->data, 1, params->nx, 1, params->ny, 1,
		 params->nz);  
   if(params->owd)
      free_d3tensor(dat->data2, 1, params->nx, 1, params->ny, 1,
		    params->nz); 
   free(dat);


   fprintf(stdout, "\ncaller() done!\n"); 
   return(1);

}


