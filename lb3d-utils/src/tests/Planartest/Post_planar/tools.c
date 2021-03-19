#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<sys/stat.h>
#include "main.h"


int
make_dir(Params *params, char *dir)
{
   /*  FROM '% man 2 mkdir':
    *
    *  int mkdir(const char *path, mode_t mode);
    */
   char        *path;
   mode_t      mode;

   path = (char *)malloc(PATHNAMELEN * sizeof(char));
   sprintf(path, "%s/%s", params->outpath, dir);
   mode = 00700;

   if(0 != mkdir(path, mode)) {
      fprintf(stderr, "\a\n\n\tError in creating directory %s\n"
	      "\tFunction tools.c:make_dir() called by read_data()\n"
	      "\tExiting to system...\n",
	      path);
      exit(-1);
   }
   
   free(path);
   return(1);
}



off_t
file_size(char *file)
{
   /*  FROM '% man 2 stat':
    *
    *  int stat(const char *path, struct stat *buf);
    */
   struct stat   *info;
   off_t size;
 
   info = (struct stat *)malloc(sizeof(struct stat));
   
   if(0 != stat(file, info)) {
      free(info);
      return(0);
   }
   else{
     size = info->st_size;	
     free(info);
     return(size);
     }
}







