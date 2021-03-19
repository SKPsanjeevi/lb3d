#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include "main.h"

extern int caller(Params *);
extern int print_help(char []);

/****************************************************************
 ** CODE NAME: POSTPLANAR
 ** (c) Nelido Gonzalez-Segredo
 ****************************************************************
 **
 ** (1) AIM: 
 **     Construct f(z) = (P_perp - P_parall) along axis z;
 **	Compute 2\sigma = \int_0^L f(z)dz, L: length along z;
 **	 
 ** (2) INPUTS:
 **     From input-file:
 **          L           : length along z;
 **          folder      : 
 **    	     gr_out_file :
 **          nt          : No. of timesteps simulated.
 **          meas        : Measurestep.
 **     As command line argument:
 **          [-task N]   : (Only if taskfarming) No. of input files 
 **                        to read in.
 **          -infile INP : Path of directory where code/ and output/ are.
 **          -result RES: Path of directory where results will be.
 **
 ** (3) ERROR HANDLING:
 **     Should a zero-size XDR file be found for certain time step,
 **     it'll be skipped and the file corresponding to the following
 **     time step will be read. An error message will be sent to
 **     stderr.
 **
 *******************************************************************
 */

int 
main(int argc, char *argv[])
{
   int             i;
   Params          *params;
   char            aux[160];


   if(NULL==(params=(Params *)malloc(sizeof(Params)))) {
      fprintf(stderr,"\nError in allocating mem for Params struct.\n");
      exit(-1);
   }

   params->ntask = 1;
   params->nx = params->ny = params->nz = 0;
   params->n_iteration = params->n_sci = 0;
   params->direction_x = params->direction_y = params->direction_z = 0;

   for(i=0; i<argc ; i++) {
      if(strcmp(argv[i],"-h")==0) {
	 printf("\a");
	 print_help(argv[0]);
	 free(params);
	 exit(0);
      }
      if(strcmp(argv[i],"-task")==0) {
	 params->ntask = atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-input")==0) {
	 strcpy(params->inputfile, argv[i+1]); 
      } else if(strcmp(argv[i],"-data")==0) {
	 strcpy(params->datapath, argv[i+1]);
      } else if(strcmp(argv[i],"-result")==0) {
	 strcpy(params->outpath, argv[i+1]);
      } else if(strcmp(argv[i],"-z")==0) {
	 params->direction_z = 1;
      } else if(strcmp(argv[i],"-x")==0) {
	 params->direction_x = 1;
      } else if(strcmp(argv[i],"-y")==0) {
	 params->direction_y = 1;
      }
   }
   
   /*  if all are zero */
   if((params->direction_z==0)&&(params->direction_x==0)&&(params->direction_y==0)) {
      printf("\nPLEASE SPECIFY DIRECTION"
	     "\nPERPENDICULAR TO THE INTERFACE!\n");
      exit(0);
   }

   strcpy(aux, params->inputfile);
   
   for(i=1; i<=params->ntask; i++) {
      if(params->ntask > 1)
	 sprintf(params->inputfile, "%s%d", aux, i);
      
      if(caller(params) != 1) {
	 fprintf(stderr, "\nError in caller(). Exiting to system...\n");
	 free(params);
	 exit(-1);
      }
   }  
   
   
   free(params);
   fprintf(stdout, "\nmain() done!\n"); 
   return(1);
}
   
   


