#include <stdio.h>
#include<stdlib.h>
#include "main.h"
#include "nrutil.h"

extern int    read_data(Params *);
extern int    decode(char *, Params *, unsigned long, 
		     float ***, Stats *);
extern int    decode_hdf(char *, Params *, unsigned long, 
		     float ***, Stats *);
extern int    extrema(Stats *, float *);
extern int    plinteg(unsigned long, Params *, float *, float *);
extern float  ***f3tensor(long, long, long, long, long, long);
extern float  *vector(long, long);
extern void   free_f3tensor(float ***, long, long, long, long, long, long);
extern void   free_vector(float *, long, long);

int
caller(Params *params)
{
   unsigned long   ts;
   static Stats    *stats;
   float           ***pxx, ***pyy, ***pzz;
   float           *profile;
   float           sigma = 0.0;
   int             ix, iy, iz, dec;
   FILE            *fp;
   char tcheckfile[180]="";



   if(NULL==(stats=(Stats *)malloc(sizeof(Stats)))) {
      fprintf(stderr,"\n\tError in allocating mem for Stats struct.");
      exit(-1);
   }

   if(read_data(params) != 1) {
      fprintf(stderr,"\n\tError in read_data() function.");
      exit(-1);
   }   

   /*  This file will contain "pressperp" data */
   
   if(VERBOSE)
	   printf("\ncaller.c:params->nz=%d", params->nz);

   pxx = f3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   pyy = f3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   pzz = f3tensor(1, params->nx, 1, params->ny, 1, params->nz);
   if(params->direction_z) {
      profile = vector(1, params->nz);
   } else if(params->direction_x) {
      profile = vector(1, params->nx);
   } else if(params->direction_y) {
      profile = vector(1, params->ny);
   }

   /*  profile must be understood as 
    *  Pzz-Pxx for (params->direction_z == TRUE), 
    *  Pxx-Pyy for (params->direction_x == TRUE), and
    *  Pyy-Pzz for (params->direction_y == TRUE).
    *
    */
 
   for(ts = 0; ts <= params->n_iteration; ts += params->n_sci) {
	if(ts< (params->n_iteration)){
		sprintf(tcheckfile, "%s.t%lu", params->checkfile, ts);
		}
   	else{
		sprintf(tcheckfile, "%s", params->checkfile);
		}
	
	if(NULL==(fp = fopen(tcheckfile,"w"))) {
		perror("\nError when opening checkfile!");
		exit(-1);
		}
	   


   /* printf("\n\tts = %lu", ts); */
      do {
	 dec = 0;
	 dec += decode_hdf("pxx", params, ts, pxx, stats);
	 dec += decode_hdf("pyy", params, ts, pyy, stats);
	 dec += decode_hdf("pzz", params, ts, pzz, stats);
	 if(dec != 3) {
	    fprintf(stderr, "\n\tSkipping timestep %lu", ts);
	    ts += params->n_sci;
	 }
      } while((dec != 3) && (ts <= params->n_iteration));
      if(ts > params->n_iteration)
	 continue;

      if(VERBOSE) {
	 printf("\npzz[2][2][32]-pxx[2][2][32]=%g"
		"\npzz[2][2][128]-pxx[2][2][128]=%g",
		pzz[2][2][32]-pxx[2][2][32],
		pzz[2][2][128]-pxx[2][2][128]);
      }

      if(params->direction_z) {
	 ix = params->nx / 2;
	 iy = params->ny / 2;
	 
	 if(VERBOSE)
		 printf("\ncaller.c: (Just before writing on pressperp*)"
				 " params->nz=%d", params->nz);

	 for(iz = 1; iz <= params->nz; iz++) {
	    profile[iz] = pzz[ix][iy][iz] - pxx[ix][iy][iz];
	    //if(ts == params->n_iteration) {
	       fprintf(fp, "\n%d\t%16.8f\t%16.8f\t%16.8f\t%16.8f", 
		       iz, pxx[ix][iy][iz], pyy[ix][iy][iz], pzz[ix][iy][iz],
		       profile[iz]);
	    //}
	 } 
      } else if(params->direction_x) {
	 iy = params->ny / 2;
	 iz = params->nz / 2;
	 
	 for(ix = 1; ix <= params->nx; ix++) {
	    profile[ix] = pxx[ix][iy][iz] - pyy[ix][iy][iz];
	    //if(ts == params->n_iteration) {
	       fprintf(fp, "\n%d\t%16.8f\t%16.8f\t%16.8f\t%16.8f", 
		       ix, pxx[ix][iy][iz], pyy[ix][iy][iz], pzz[ix][iy][iz],
		       profile[ix]);
	    //}
	 }
      }  else if(params->direction_y) {
	 iz = params->nz / 2;
	 ix = params->nx / 2;
	 
	 for(iy = 1; iy <= params->ny; iy++) {
	    profile[iy] = pyy[ix][iy][iz] - pzz[ix][iy][iz];
	    //if(ts == params->n_iteration) {
	       fprintf(fp, "\n%d\t%16.8f\t%16.8f\t%16.8f\t%16.8f", 
		       iy, pxx[ix][iy][iz], pyy[ix][iy][iz], pzz[ix][iy][iz],
		       profile[iy]);
	    //}
	 }
      }
      
      if(plinteg(ts, params, profile, &sigma) != 1) {     
	 fprintf(stderr,"\n\tError in plinteg() function");
	 exit(-1);
      }

   fclose(fp);	      
   } /* Next ts */
   
   

   fflush(stderr);
   free(stats); 
   free_f3tensor(pzz, 1, params->nx, 1, params->ny, 1, params->nz); 
   free_f3tensor(pxx, 1, params->nx, 1, params->ny, 1, params->nz); 
   free_vector(profile, 1, params->nz);

   fprintf(stdout, "\ncaller() done!\n"); 
   return(1);

}


