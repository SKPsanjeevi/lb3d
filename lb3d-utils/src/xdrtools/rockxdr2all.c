#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <math.h>

/* Converts a rock file from xdr to ascii format */
/* Only the rock sites (stone=1) with their lattice coordinates are printed on stndout (sparse rock) */

int main (int argc, char *argv[])
{
  FILE	*f;
  XDR	xdrs;
  unsigned int nx,ny,nz;
  unsigned int x,y,z;
  int c;
  
  if (5!=argc) {
    fprintf(stderr,"Syntax: %s <nx> <ny> <nz> <xdrfile>\n",argv[0]);
    return -1;
  }
  
  nx=(unsigned int)atoi(argv[1]);
  ny=(unsigned int)atoi(argv[2]);
  nz=(unsigned int)atoi(argv[3]);

  if (NULL==(f=fopen(argv[4],"r"))) {
    perror("fopen()");
    return -1;
  }

  xdrstdio_create(&xdrs,f,XDR_DECODE);

  printf("hello\n");

  for (z = 1; z <= nz; z++)
    {
      for (y = 1; y <= ny; y++)
        {
          for (x = 1; x <= nx; x++)
            {
                if (1!=xdr_int(&xdrs,&c)) {
                printf("%d\n",c);
                perror("xdr_int()");
                return -1;   
              }         
              if(c==1) printf("%d %d %d %d\n",x,y,z,c);
	    }
	}
    }

  printf("done\n");

  xdr_destroy(&xdrs);
  fclose(f);

  return 0;
}
