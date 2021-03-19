#include <stdio.h>
#include "main.h"


extern int caller(Params *);
extern int print_help(char []);




int 
main(int argc, char *argv[])
{
   static Params   *params;

   if(NULL==(params=(Params *)malloc(sizeof(Params)))) {
      fprintf(stderr,"\nError in allocating mem for Params struct.\n");
      exit(-1);
   }

   params->nx = params->ny = params->nz = 0;

   /*
     IF FORGOT ONE ARGUMENT
     WARN AND EXIT TO AYATEM...
   */
   if((argc<4) && !(strcmp(argv[1],"-h")==0))
   {
      fprintf(stderr,"\n\nYou need to specify all the parameters\n"
	             "explained in the help screen!\n");
      fflush(stderr);
      exit(0);
   }

   
   if(strcmp(argv[1],"-h")==0) {
      printf("\a");
      print_help(argv[0]);
      free(params);
      exit(0);
   }
   else {
      params->nx = atoi(argv[1]);
      params->ny = atoi(argv[2]);
      params->nz = atoi(argv[3]);
      strcpy(params->inpfile, argv[4]);	 
   }   
   
   
   if(caller(params) != 1) {
      fprintf(stderr, "\nError in caller(). Exiting to system...\n");
      free(params);
      exit(-1);
   }
   
   
   free(params);
   fprintf(stdout, "\nmain() done!\n"); 
   return(1);
}
   
