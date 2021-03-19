#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* This program can be used to figure out the
 * appropriate job size for HPCx.
 *
 * It takes an optional lattice size on the command line,
 * which defaults to 1024. It then calculates all possible
 * (number of LPARs / no of CPUs per LPAR) combinations
 * which will run that lattice size on LB3D, with an equal chunk
 * of the lattice given to each CPU.
 *
 * On linux:
 * gcc -Wall -pedantic -o hpcx-calc hpcx-calc.c -lmpich
 */

int main(int argc, char *argv[])
{
	int rank;
	int dims[3]={0,0,0};
	int lsize=1080; /* Size of lattice */
	int nlpars,frac,ncpus;
	long totalmem;
	long mem_per_cpu;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (argc>1) {
		lsize=atoi(argv[1]);
	}
	printf("Lattice size %d^3\n",lsize);
	totalmem = lsize*lsize*lsize/1024; /* MB */

	for (nlpars=1;nlpars<=160;nlpars++) {
		for (frac=1;frac<=8;frac++) {
			ncpus = frac*nlpars;
			mem_per_cpu = totalmem/ncpus;
			dims[0]=dims[1]=dims[2]=0;
			MPI_Dims_create(ncpus,3,dims);
			if (
				((lsize%(dims[0]))==0) &&
				((lsize%(dims[1]))==0) &&
				((lsize%(dims[2]))==0)
			   ) {
				printf(
				"%d LPAR, (%d/8)= %d CPUs = %d x %d x %d",
					nlpars,frac,
					ncpus,dims[0],dims[1],dims[2]);
				printf(" Split %d x %d x %d mem %ld\n",
					lsize/dims[0],
					lsize/dims[1],
					lsize/dims[2],
					mem_per_cpu
					);
			}

		}
	}

	
	return 0;
}
