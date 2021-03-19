/*  3-D VERSION OF THE LIKE
 *  PROGRAM BY J. CHIN
 */

/* reads in .xdr rock file and outputs a raw binary (e.g. for */
/* viewing with ImageJ  (import 8bit)*/
/* FR 07 / 2007 */

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <math.h>

#define BUFSIZE 1024 /* So shoot me. */


int main (int argc, char *argv[])
{
  FILE	*f=NULL;
  int	i=0;
  XDR	xdrs;
  char *fname=NULL;
  int nx=0,ny=0,nz=0,cx=0,cy=0,cz=0;
  
  int cnt=0;
  int out=0;
  char outfile[80];
  FILE *rawout=NULL;


#ifdef SGL
  float	buf=0.0;
  float phimax=0.0,phimin=0.0,phisum=0.0,phisum2=0.0;
  float sd=0.0;
#else
  double buf=0.0;
  double phimax=0.0,phimin=0.0,phisum=0.0,phisum2=0.0;
  double sd=0.0;
#endif

  if (argc!=5) {
    fprintf(stderr,"Syntax: %s <filename> <nx> <ny> <nz>\n",argv[0]);
    return -1;
  }
  
  fname=argv[1];
  sprintf(outfile, "%s.bin", fname);
  
  nx=atoi(argv[2]);
  ny=atoi(argv[3]);
  nz=atoi(argv[4]);
  
  /*printf("DIM check: %d %d %d \n",nx,ny,nz);*/
  
  if (NULL==(f=fopen(fname,"r"))) {
    fprintf(stderr, "Unable to open input file %s\n", fname);
    return -1;
  }
  

  if (NULL==(rawout=fopen(outfile,"w"))) {
    fprintf(stderr, "Unable to open output  file %s \n", outfile);
    return -1;
  }

  xdrstdio_create(&xdrs,f,XDR_DECODE);
	
  
  for(cz=0;cz<nz;cz++)
    { 
      for(cy=0;cy<ny;cy++)
	{ 
	  for(cx=0;cx<nx;cx++)
	    { 
	
	      i++;
	      //	    for (i=0;i<(nx*ny*anz);i++) {
#ifdef SGL
	      xdr_float(&xdrs,&buf);
#else
	      xdr_double(&xdrs,&buf);
#endif
	      phisum  += buf;
	      phisum2 += buf*buf;
	      if (buf>phimax) { phimax=buf;}
	      if (buf<phimin) { phimin=buf;}
	      
	      out = ((int) (buf==5.)?254:0);
	      fputc(out,rawout);
	      /*  printf("[nz:%d,ny:%d,nx:%d]:  value: %f --",cz,cy,cx,buf);
		  printf("i: %d value: %d \n",i,out);*/
	      
	    }
	}
    }
  xdr_destroy(&xdrs);
  fclose(rawout);
  
  
#ifdef SGL
  phisum/=(float)(nx*ny*nz);
  phisum2/=(float)(nx*ny*nz);
#else
  phisum/=(double)(nx*ny*nz);
  phisum2/=(double)(nx*ny*nz);
#endif
  sd=sqrt(phisum2-phisum*phisum);
  
  
  printf("Sides copied: %d \n",i);
	printf("Max %f Min %f Mean %f SD %f\n",
	       phimax,phimin,phisum,sd);
	fclose(f);
	return 0;
}
