#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "vol2df.h"

void usage(char *s)
{
    fprintf(stderr,"Usage: %s <nx> <ny> <filename>\n",s);
}

int main(int argc, char *argv[])
{
    int nx,ny;
    FILE *f;
    struct binnery *b=NULL;
    float k1,k1_err,l0;
    struct vol2d *v=NULL;
    if (4!=argc) { usage(argv[0]); exit(-1); }

    nx=atoi(argv[1]);
    ny=atoi(argv[2]);

    if (NULL==(f=fopen(argv[3],"rb")))
        { perror("fopen"); exit(-1); }

    if (NULL==(v=vol2df_new_from_xdrdfh(nx,ny,f)))
        { fprintf(stderr,"vol2df_new_fromxdrfh failed.\n"); exit(-1); }

    if (NULL==(b=vol2df_sfactor(v))) { return -1; }

    binnery_scaled_dump(b,stdout,2.0*M_PI/v->nx);
    binnery_scaled_1stmoment(b,2.0*M_PI/v->nx,&k1,&k1_err);
    l0=2.0*M_PI/k1;
    printf("# k1= %.10f err %.10f L= %.10f \n",k1,k1_err,l0);
    binnery_destroy(b);

    return 0;
}
