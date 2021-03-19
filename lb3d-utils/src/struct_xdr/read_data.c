#include <stdio.h>
#include <string.h>

#include "main.h"


int read_data(Params *params)
{
  int    nxflag=0, nyflag=0, nzflag=0, folderflag=0, 
         gr_out_fileflag=0, sci_startflag=0;
  size_t l;
  char   *line, *folder, *gr_out_file;
  FILE   *infile;

  char *trimstring(char s[], int lim);


/*   printf("\nIn read_data(), just after entering \n" */
/* 	 "params->nx = %lu\n" */
/* 	 "params->ny = %lu\n" */
/* 	 "params->nz = %lu\n",  */
/* 	 params->nx, params->ny, params->nz); */
  /* Open input-file */
  if (NULL == (infile = fopen(params->inpfile,"r"))) {
     fprintf(stderr, "Could not open %s\n"
	     "Error in read_data()\n", params->inpfile);
     exit(-1);
  }
  printf("Reading input file %s\n", params->inpfile);

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
	fprintf(stderr, "In file %s\n", params->inpfile);
	fclose(infile);
	return(-1);
	}
      */
      if(sscanf(line, "ny = %d", &(params->ny)) == 1)
	nyflag = 1;
      /*
	{
	  perror("Error in reading lattice size ny!\n");
	  fprintf(stderr, "In file %s\n", params->inpfile);
	  fclose(infile);
	  return(-1);
	}
      */
      if(sscanf(line, "nz = %d", &(params->nz)) == 1)
	nzflag = 1;
	/*
	{
	  perror("Error in reading lattice size nz!\n");
	  fprintf(stderr, "In file %s\n", params->inpfile);
	  fclose(infile);
	  return(-1);
	}	
	*/   
      if(sscanf(line, "n_sci_start = %lu", &(params->sci_start)) == 1)
	 sci_startflag = 1;
      if(sscanf(line, "folder = '%s", folder) == 1)
	folderflag = 1;
      /*
	{
	perror("Error in reading folder name 'folder'!\n");
	fprintf(stderr, "In file %s\n", params->inpfile);
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
	  /*
	{
	perror("Error in reading file prefix name 'gr_out_file'!\n");
	fprintf(stderr, "In file %s\n", params->inpfile);
	fclose(infile);
	return(-1);
	}
      */

    }
  
  if(!nxflag) 
    fprintf(stderr, "Couldn't read nx from %s\n", params->inpfile);
  if(!nyflag) 
    fprintf(stderr, "Couldn't read ny from %s\n", params->inpfile);
  if(!nzflag) 
    fprintf(stderr, "Couldn't read nz from %s\n", params->inpfile);
  if(!folderflag) 
    fprintf(stderr, "Couldn't read folder from %s\n", params->inpfile);
  if(!gr_out_fileflag) 
    fprintf(stderr, "Couldn't read gr_out_file from %s\n", params->inpfile);
  if(!sci_startflag) 
    fprintf(stderr, "Couldn't read sci_start from %s\n", params->inpfile);

  /* REMOVE TRAILING "'" */
  
  l = strcspn(folder, "'");  folder[l] = '\0';
  /* printf("\nl=%d folder=%s", l, folder); */
  l = strcspn(gr_out_file, "'");  gr_out_file[l] = '\0';
  /* printf("\nl=%d gr_out_file=%s", l, gr_out_file); */


  printf("\nAs read from input file:\n"
	 "\tnx=%d\n"
	 "\tny=%d\n"
	 "\tnz=%d\n"
	 "\tXDR data folder: %s\n"
	 "\tStudy name: %s\n"
	 "\tsci_start=%lu\n",
	 params->nx, params->ny, params->nz,
	 folder, gr_out_file, params->sci_start);

  strcpy(params->folder, folder);
  strcpy(params->studyname, gr_out_file);

  



  /*
    FOLLOWING LINE IS NEEDED TO ENABLE WRITING 
    RESULT FILES ONTO SAME DIRECTORY OF XDR FILES...
  */
  /*
    strcpy(params->outfile, params->datafile);
  */





  
  /* For now, e.g.:      'test/colour_test'
     We want, e.g.:      'test/colour_test_t0001.xdr'
     */


  free(line);
  free(folder);
  free(gr_out_file);
  return(1);
}


