#include <stdio.h>
#include "main.h"
#include "nrutil.h"

extern int    decode(Params *, Dat *, Stats *);
extern int    fft(double ***, Params *, Stats *);
extern int    circ_avg(double ***, Params *, CASF *);
extern int    domain_size(CASF *, Params *, Stats *);
extern int    dump_to_files(CASF *, Params *, Stats *);
extern int    *i1t(int n1);
extern double *d1t(int n1);
extern double ***d3tensor(long, long, long, long, long, long);
extern void   free_d3tensor(double ***t, long, long, long, long, long,
			    long);

int           dump_second_field, dump_third_field;
int		dec;

int
caller(Params *params)
{
   Dat           *dat;
   Stats         *stats;
   CASF          *casf;
   size_t        l;
   char          aux[160];
   


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


   dat->data = d3tensor(1, params->nx, 1, params->ny, 1, params->nz);

/* Allocate mem for SF vectors: */

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

/* Construct outfile: */

/*   l = strcspn(params->inpfile,"."); */
   l = strlen(params->inpfile) - 4;
   strcpy(aux, params->inpfile);   
   aux[l] = '\0';       /* Important to use single quotes! */
   sprintf(params->outfile, "%s%s", aux, ".sf");
/*   params->outfile[l] = ".sf"; */
   if(VERBOSE)
      printf("\noutfile = %s", params->outfile);


   dec = decode(params, dat, stats);

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
   if(dump_to_files(casf, params, stats) != 1) {
      fprintf(stderr,"Error in domain_size() function\n");
      exit(-1);
   }
   
   fflush(stderr);
   free(stats); 
   free(casf);
   free_d3tensor(dat->data, 1, params->nx, 1, params->ny, 1,
		 params->nz);  
   free(dat);
   free(casf->S);   /* CAREFUL WITH THIS WHEN CALLING caller() MORE
		       THAN ONCE.   WHY???*/
   free(casf->sum_shell);

   fprintf(stdout, "\ncaller() done!\n"); 
   return(1);

}


