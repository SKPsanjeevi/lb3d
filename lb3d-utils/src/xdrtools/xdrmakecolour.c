/* XDR format both input/output
 * Read od and wd files and make a colour file: od-wd 
 *
 * Syntax: xdrmakecolour nx ny nz od_filename wd_filename col_filename
 *
 */

#include <stdio.h>
#include <rpc/rpc.h>

#define ARGV_NX	1
#define ARGV_NY	2
#define ARGV_NZ	3
#define ARGV_ODFILE 4
#define ARGV_WDFILE 5
#define ARGV_COLFILE 6


int main (int argc, char *argv[])
{
  int nx,ny,nz;
  int x,y,z;
  
  /* NB: type "double" should be at least 64 bits wide.
   * if not, it's trivial to change this code to use long double
   * for the buffers.
   */

#ifdef SGL
  float	od,wd,col;
#else
  double od,wd,col;
#endif

  FILE *odfile,*wdfile,*colfile;
  XDR	xdrs1,xdrs2,xdrs3;


  /* Ensure correct no of arguments */

  if (7!=argc) {
    fprintf(stderr,
	    "Syntax: %s <nx> <ny> <nz> <od_file> <wd_file> <col_file>\n",
	    argv[0]);
    return -1;
  }

  /* Parse the dimension arguments */
  
  nx = atoi(argv[ARGV_NX]);
  ny = atoi(argv[ARGV_NY]);
  nz = atoi(argv[ARGV_NZ]);

  /* Sanity check them */
  
  if ( (nx<=0) || (ny<=0) || (nz<=0) ) {
    fprintf(stderr,"Invalid size %d %d %d\n",nx,ny,nz);
    return -1;
  }

  /* Open the input files */

  if (NULL==(odfile=fopen(argv[ARGV_ODFILE],"r"))) {
    perror("Could not open od_file \n ");
    return -1;

  }

  if (NULL==(wdfile=fopen(argv[ARGV_WDFILE],"r"))) {
    perror("Could not open wd_file \n");
    return -1;
  }

  /* Open the output file */

  if (NULL==(colfile=fopen(argv[ARGV_COLFILE],"w"))) {
    perror("Could not open col_file \n");
    return -1;
  }

  /* Set up an XDR stream */

  xdrstdio_create(&xdrs1,odfile,XDR_DECODE);
  xdrstdio_create(&xdrs2,wdfile,XDR_DECODE);
  xdrstdio_create(&xdrs3,colfile,XDR_ENCODE);

  /* Now read and parse stdin. */


  for (z=1;z<=nz;z++) {
    for (y=1;y<=ny;y++) {
      for (x=1;x<=nx;x++) {
#ifdef SGL
	xdr_float(&xdrs1,&od);
	xdr_float(&xdrs2,&wd);
	col = od-wd;
	xdr_float(&xdrs3,&col);
#else
	xdr_double(&xdrs1,&od);
	xdr_double(&xdrs2,&wd);
	col = od-wd;
	xdr_double(&xdrs3,&col);
#endif
      }
    }
  }


  /* Close the XDR stream */

  xdr_destroy(&xdrs1);
  xdr_destroy(&xdrs2);
  xdr_destroy(&xdrs3);

  return 0; /* And they all lived happily ever after. */

}

