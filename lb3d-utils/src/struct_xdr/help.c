#include<stdio.h>

int
print_help(char executable[40])
{
   printf("\n\n\t<<<<<<<<<<<<<<< HELP ON struct_xdr"
	  ">>>>>>>>>>>>>>>>>>>>\n\n"
	  "\tSyntax:\n\n"
	  "\t%s [-task N] <-file FILETYPE> <-nt NT> <-meas MEAS>\n"
	  "\t\t<-inpfile INP> <-data DATA> <-res RES>\n\n"
	  "\twhere:\n\n"
	  "\t* N is the number of input files to read in\n"
	  "\t  (taskfarming case): input-file1 ... input-fileN.\n"
	  "\t* FILETYPE can be 'colour', 'sur', 'od', 'wd', 'owd' or 'arr'.\n"
	  "\t* NT is no. of time steps simulated (n_iteration).\n"
	  "\t* MEAS is the measurestep.\n"
	  "\t* INP is the input file name (with path).\n"
	  "\t* DATA is the subdirectory where 'folder' is\n"
	  "\t  (with path). ('folder' contains the XDR files.)\n"
	  "\t* RES is the subdirectory where results are to\n"
	  "\t  be written (with path).\n\n\n"
	  "\tPLEASE NOTE:\n\n"
	  "\t* Files produced are .sf and .siz.\n"
	  "\t* If an XDR datafile has zero length or is inaccessible\n"
	  "\t  it will be skipped and the file corresponding to the\n"
	  "\t  next timestep will be read. An error message will be\n"
	  "\t  sent to stderr.\n" 
	  "\t* Path names: Don't use '~' or the trailing '/';\n"
	  "\t  you can use relative paths.\n"
	  "\t* This program only deals with lattice sizes being\n"
	  "\t  powers of 2.\n"
	  "\t* There is one .sf file per measurestep, \n"
	  "\t  containing (k,SF) data pairs.\n"
	  "\t* File .siz contains:\n"
	  "\t\tts k1 errk1 k2 errk2 R1 errR1 R2 errR2 min max avg\n"
	  "\t  where ts is the time step, k? are the ?th moments of\n"
	  "\t  S(k,t) with respect to k, R? the domain size derived\n"
	  "\t  from k?, resp., and err?? are the relevant\n"
	  "\t  uncertainties derived from spherical averaging.\n" 
	  "\t* File .lsz contains (ln k,ln R1) data pairs.\n"
	  "\t* It is recommended to redirect stdout and stderr\n"
	  "\t  to a file for later inspection: info on errors and\n"
	  "\t  unread XDR files may be in them.\n\n\n",
	  executable);
   

   return(1);
}
