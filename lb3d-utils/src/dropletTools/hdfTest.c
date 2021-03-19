/// hdfTest.c
///
/// Author: Sebasstian Schmieschek
///
/// Function to test hdf read/write routine
/// Args: inputFilename, outputFilename
/// Reads inputFile to array. writes array to outputFile


#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <regex.h>

#include "hdfread.h"


int
main(int argc,
     char *argv[])
{
	// This function shall serve testing puposes for the hdf functions only

	char		*FILE, *OUTFILE;
	FLOATNUM	***data;
	int		nx, ny, nz;
	
	FILE = argv[1];
	OUTFILE = argv[2];

	hdfread2(FILE,&data);
	getDims(FILE,&nx,&ny,&nz);
	hdfwrite(OUTFILE,nx,ny,nz,data);

	return(0);
}

