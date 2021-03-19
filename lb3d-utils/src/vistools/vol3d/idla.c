#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "vol3df.h"


struct vol3d *dla_new(int nx, int ny, int nz)
{
    struct vol3d *self=NULL;

    if (NULL==(self=vol3df_new_zero(nx,ny,nz)))
        { fprintf(stderr,"dla_new: vol3df_new_zero failed.\n"); return NULL; }

    /* Set seed point */

    *(vol3df_wrap_ptr(self,nx/2,ny/2,nz/2)) = 1.0;

    return self;

}

/* Return a random number between 0 and (max-1) inclusive */
int random_whole(int max)
{
    int eta;

    eta = (int)(((float)max)*rand()/(RAND_MAX+1.0));
    return eta;
}

/* Randomly return 1 or 0 */
int coinflip(void)
{
    return (0==(rand()&0x10000000));
}

/* Randomly return +1 or -1 */
int plusminus(void)
{
    if (0==(rand()&0x10000000))
    { return 1; } else { return -1; }
}

/* Return 1 if the given point is next or on top of another point, else
 * return zero.
 */
int dla_testPoint(struct vol3d *self, int x, int y, int z)
{
    int dx,dy,dz;

    for (dz=-1;dz<=1;dz++) {
    for (dy=-1;dy<=1;dy++) {
    for (dx=-1;dx<=1;dx++) {
        if (0!=vol3df_wrap(self,x+dx,y+dy,z+dz)) { return 1; }
    }}}

    return 0;

}

float dla_dist(struct vol3d *self, int x, int y, int z)
{
    float distmax;
    float dist;
    int x0,y0,z0;

    x0=self->nx/2;
    y0=self->ny/2;
    z0=self->nz/2;

    distmax = sqrt(x0*x0+y0*y0+z0*z0);

    dist = sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) );

    return dist/distmax;

}


/* Choose a random point in the volume. If the voxel at that point is
 * nonzero, then increment it and return. Otherwise, take a random step
 * in one of the Cartesian diractions, and repeat.
 */
int dla_iterate(struct vol3d *self)
{
    int x,y,z;
    long int count=0;

    x = self->nx/2;
    y = self->ny/2;
    z = self->nz/2;

    while (0.0!=vol3df_wrap(self,x,y,z)) {
        count++;
        x += plusminus();
        y += plusminus();
        z += plusminus();
    }

    *(vol3df_wrap_ptr(self,x,y,z)) += count;

    return 0;

}


int main(int argc, char *argv[])
{

    struct vol3d *self;
    int t;
    int iostep=1000;
    char buf[1000];

    if (NULL==(self=dla_new(128,128,128)))
        { fprintf(stderr,"dla_new failed.\n"); return -1; }

    srand(time(NULL));

    for (t=0;t<10000;t++) {
        if (0==(t%iostep)) {
            snprintf(buf,100,"dla-%05d.h5",t);
            if (0!=vol3df_write_hdf5(self,buf,"/OutArray")) {
                fprintf(stderr,"Failed to write \"%s\".\n",buf);
            } else {
                fprintf(stderr,"Wrote \"%s\".\n",buf);
            }
        }
        dla_iterate(self);
    }

    return 0;
}
