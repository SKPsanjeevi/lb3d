#include <stdio.h>
#include "main.h"


extern int caller(Params *);
extern int print_help(char []);




int 
main(int argc, char *argv[])
{
   int             i;
   static Params   *params;
   char            aux[160];

   if(NULL==(params=(Params *)malloc(sizeof(Params)))) {
      fprintf(stderr,"\nError in allocating mem for Params struct.\n");
      exit(-1);
   }

   params->ntask = 1;
   params->nx = params->ny = params->nz = 0;
   params->owd = 0;

   /*
     IF FORGOT ONE ARGUMENT
     WARN AND EXIT TO AYATEM...
   */
   if((argc<7) && !(strcmp(argv[1],"-h")==0))
   {
      fprintf(stderr,"\n\nYou need to specify all the parameters\n"
	             "explained in the help screen!\n");
      fflush(stderr);
      exit(0);
   }

   for(i=0; i<argc ; i++) {
      if(strcmp(argv[i],"-h")==0) {
	 printf("\a");
	 print_help(argv[0]);
	 free(params);
	 exit(0);
      }
      if(strcmp(argv[i],"-task")==0) {
	 params->ntask = atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-file")==0) {
	 strcpy(params->filetype, argv[i+1]);
	 if(strcmp(argv[i+1],"owd")==0)
	    params->owd = 1;
	 if(strcmp(argv[i+1],"arr")==0)
            params->arr = 1;
      } else if(strcmp(argv[i],"-nt")==0) {
	 params->n_iteration = (unsigned long)atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-meas")==0) {
	 params->n_sci = (unsigned long)atoi(argv[i+1]);
      } else if(strcmp(argv[i],"-inpfile")==0) {
	 strcpy(params->inpfile, argv[i+1]);
      } else if(strcmp(argv[i],"-data")==0) { 
	 strcpy(params->datapath, argv[i+1]);
	 strcpy(params->outpath, params->datapath);
      } else if(strcmp(argv[i],"-res")==0) { 
	 strcpy(params->outpath, argv[i+1]);
      } 
   }
   
   
   strcpy(aux, params->inpfile);
   for(i=1; i<=params->ntask; i++) {
      if(params->ntask > 1)
	 sprintf(params->inpfile, "%s%d", aux, i);
  
      if(caller(params) != 1) {
	 fprintf(stderr, "\nError in caller(). Exiting to system...\n");
	 free(params);
	 exit(-1);
      }

/*       printf("\nIn main(), just after returning from caller(params) \n" */
/* 	     "params->nx = %lu\n" */
/* 	     "params->ny = %lu\n" */
/* 	     "params->nz = %lu\n",  */
/* 	     params->nx, params->ny, params->nz); */

   }  
   
   
   free(params);
   fprintf(stdout, "\nmain() done!\n"); 
   return(1);
}
   
   




