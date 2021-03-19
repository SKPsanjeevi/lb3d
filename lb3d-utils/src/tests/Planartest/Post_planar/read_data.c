 #include <stdio.h>
#include<stdlib.h>
#include <string.h>
#include "main.h"


int read_data(Params *params)
{
  int    nxflag, nyflag, nzflag, folderflag, gr_out_fileflag,
         iteraflag, sciflag;
  size_t l;
  char   *line, *folder, *gr_out_file;
  FILE   *infile;


  char *trimstring(char s[], int lim);

  /* Open input-file */
  if (NULL == (infile = fopen(params->inputfile,"r"))) {
     fprintf(stderr, "Could not open %s\n"
	     "Error in read_data()\n", params->inputfile);
     exit(-1);
  }
  printf("Reading input file %s\n", params->inputfile);

  line = (char *)malloc(MAXLINE * sizeof(char));
  folder = (char *)malloc(NAMELEN * sizeof(char));
  gr_out_file = (char *)malloc(NAMELEN * sizeof(char));
 

  /* If line, folder and gr_out_file ('cs') 
     are defined as char cs[MAXLINE], fgets() 
     and str functions cannot
     assign them the (char *) it returns.
     */

  /* Now use of sscanf as it's the only one that
     deals with blank spaces in format string.
     The upper limit of i in next loop cannot exceed
     no of lines in input-file.
     Use of fgets() avoids using getline().
     */
  /* while(getline(line, sizeof(line)) > 0) */

  /* Difference with read_data_v1.c is in following lines */


  while(!feof(infile)) 
    {  
      line = fgets(line, MAXLINE, infile);
      /* It didn't work with sizeof(line)
	 instead of MAXLINE.
	 Returns line; on end of file or error, returns NULL.
      */
      
      if(sscanf(line, "nx = %d", &(params->nx)) == 1)
	nxflag = 1;
      /*
	{
	perror("Error in reading lattice size nx!\n");
	fprintf(stderr, "In file %s\n", params->inputfile);
	fclose(infile);
	return(-1);
	}
      */
      if(sscanf(line, "ny = %d", &(params->ny)) == 1)
	nyflag = 1;
      /*
	{
	  perror("Error in reading lattice size ny!\n");
	  fprintf(stderr, "In file %s\n", params->inputfile);
	  fclose(infile);
	  return(-1);
	}
      */
      if(sscanf(line, "nz = %d", &(params->nz)) == 1) {
	nzflag = 1;
	if(VERBOSE)
	   printf("\nJust read: nz=%d", params->nz);	
        }
	/*
	{
	  perror("Error in reading lattice size nz!\n");
	  fprintf(stderr, "In file %s\n", params->inputfile);
	  fclose(infile);
	  return(-1);
	}	
	*/   
      if(sscanf(line, "folder = '%s", folder) == 1)
	folderflag = 1;
      /*
	{
	perror("Error in reading folder name 'folder'!\n");
	fprintf(stderr, "In file %s\n", params->inputfile);
	fclose(infile);
	return(-1);
	}
      */
      /* printf("folder = %s", folder); fflush(stdout); */
      
      if(sscanf(line, "gr_out_file = '%s", gr_out_file) == 1)
	{
	  gr_out_fileflag = 1;
	  break;
	}
      if(sscanf(line, "n_iteration = %lu", &(params->n_iteration)) == 1)
	 iteraflag = 1;
      if(sscanf(line, "n_sci = %lu", &(params->n_sci)) == 1)
	 sciflag = 1;


	  /*
	{
	perror("Error in reading file prefix name 'gr_out_file'!\n");
	fprintf(stderr, "In file %s\n", params->inputfile);
	fclose(infile);
	return(-1);
	}
      */
    }

  if(!nxflag) 
    fprintf(stderr, "Coudn't read nx from %s\n", params->inputfile);
  if(!nyflag) 
    fprintf(stderr, "Coudn't read ny from %s\n", params->inputfile);
  if(!nzflag) 
    fprintf(stderr, "Coudn't read nz from %s\n", params->inputfile);
  if(!folderflag) 
    fprintf(stderr, "Coudn't read folder from %s\n", params->inputfile);
  if(!gr_out_fileflag) 
    fprintf(stderr, "Coudn't read gr_out_file from %s\n", params->inputfile);
  if(!iteraflag) 
     fprintf(stderr, "Coudn't read n_iteration from %s\n",
	     params->inputfile);
  if(!sciflag) 
     fprintf(stderr, "Coudn't read n_sci from %s\n",
	     params->inputfile);

  /* REMOVE TRAILING "'" */
  
  l = strcspn(folder, "'");  folder[l] = '\0';
  /* printf("\nl=%d folder=%s", l, folder); */
  l = strcspn(gr_out_file, "'");  gr_out_file[l] = '\0';
  /* printf("\nl=%d gr_out_file=%s", l, gr_out_file); */


  /*  IMPORTANT: sprintf() can't concatenate two strings
     and store the result in one of them. */

  /*  NOTE!!
   *  In this version I'm not using params->datapath
   *  as I'm assuming XDR files are in lbe/output/folder
   *
   */
  /*
    sprintf(params->datafile, "%s/%s_%s", params->datapath,
	  params->filetype, gr_out_file);
	  */

  params->folder = folder;
  params->gr_out_file = gr_out_file;
  if(params->ntask > 1) {
     make_dir(params, folder);
     sprintf(params->outfile, "%s/%s/sigma_%s.dat", params->outpath,
 	     folder,gr_out_file);
  }
  else {
     sprintf(params->outfile, "%s/sigma_%s.dat", params->outpath,
	     gr_out_file);
     sprintf(params->checkfile, "%s/pressprep_%s.dat", params->outpath,
	     gr_out_file);
  }

  printf("\nAs read from input file:\n"
	 "\tnx =             %d\n"
	 "\tny =             %d\n"
	 "\tnz =             %d\n"
	 "\tXDR data folder: %s\n"
	 "\tStudy name:      %s\n"
	 "\tn_iteration =    %lu\n"
	 "\tn_sci =          %lu\n"
	 "\tcheckfile:       %s\n"
	 "\tdirection_z      %d\n"
	 "\tdirection_x      %d\n"
	 "\tdirection_y      %d\n",
	 params->nx, params->ny, params->nz,
	 folder, gr_out_file, params->n_iteration, params->n_sci,
	 params->checkfile, 
	 params->direction_z, params->direction_x, params->direction_y);


  free(line);
  return(1);
}


