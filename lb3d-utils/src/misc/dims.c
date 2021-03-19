#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Take a number of CPUs on the command line.
 * Return the default dimensions of an MPI Cartesian
 * topology for that number of CPUs.
 *
 * On linux:
 *
 * gcc -Wall -pedantic -o dims dims.c -lmpich
 */

int main(int argc, char *argv[])
{
	int rank;
	int dims[3]={0,0,0};
	int n=0;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (argc>1) {
		n=atoi(argv[1]);
	} else { n=16; }

	MPI_Dims_create(n,3,dims);

	printf("%d = %d x %d x %d\n",n,dims[0],dims[1],dims[2]);
	
	return 0;
}
