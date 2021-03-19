#include<stdio.h>

int
print_help(char executable[40])
{
   printf("\n\n\t<<<<<<<<<<<<<<< HELP ON sigma"
	  ">>>>>>>>>>>>>>>>>>>>\n\n"
	  "\tSyntax:\n\n"
	  "\t%s [-task N] <-input INP> <-data DAT> <-result RES> <-DIR>\n\n"
	  "\twhere:\n\n"
	  "\t* N is the number of input files to read in\n"
	  "\t  (taskfarming case): input-file1 ... input-fileN.\n"
	  "\t* INP is the input file name (with path).\n"
	  "\t* DAT is the directory where the XDR data files \n"
	  "\t  are (with path).\n", executable);
   printf("\t* RES is the directory where output file\n"
	  "\t  (e.g. sigma_test.dat) is to be written (with path).\n"
	  "\t* DIR is a lowercase character denoting the direction\n"
	  "\t  perpendicular to the planar interface: x, y or z.\n"
	  "\n\n");
    printf("\tPLEASE NOTE:\n\n"
	  "\t* If TASKFARMING, subdirectories of RES are\n"
	  "\t  automatically created and resulting files \n"
	  "\t  (.sf .siz .lsz) stored in them.\n"
	  "\t* If an XDR datafile has zero length or is inaccessible\n"
	  "\t  it will be skipped and the file corresponding to the\n"
	  "\t  next timestep will be read. An error message will be\n"
	  "\t  sent to stderr.\n" );
   printf("\t* Path names: Don't use '~' or the trailing '/';\n"
	  "\t  you can use relative paths.\n"
	  "\t* This program only deals with system sizes being\n"
	  "\t  powers of 2.\n"
	  "\t* It is recommended to redirect stdout and stderr\n"
	  "\t  to a file for later inspection: info on errors and\n"
	  "\t  unread XDR files may be in them.\n\n\n");
   

   return(1);
}
