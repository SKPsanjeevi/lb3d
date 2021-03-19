#include<stdio.h>

int
print_help(char executable[40])
{
   printf("\n\n\t<<<<<<<<<<<<<<< HELP ON sfreg_xdr"
	  ">>>>>>>>>>>>>>>>>>>>\n\n"
	  "\tSyntax:\n\n"
	  "\t%s <nx> <ny> <nz> <inpfile>\n\n"
	  "\twhere:\n\n"
	  "\t* nx,ny,nz are the dimensions of the file.\n"
	  "\t* inpfile is the input file name (with path).\n"
	  "\tPLEASE NOTE:\n\n"
	  "\t* Files produced are .sf.\n"
	  "\t* Path names: Don't use '~' or the trailing '/';\n"
	  "\t  you can use relative paths.\n"
	  "\t* This program only deals with system sizes being\n"
	  "\t  powers of 2.\n"
	  "\t* There is one .sf file per measurestep, \n"
	  "\t  containing (k,SF) data pairs.\n"
	  "\t  \n\n\n",
	  executable);
   

   return(1);
}
