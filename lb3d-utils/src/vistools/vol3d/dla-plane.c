#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vol3df.h"


struct vol3d *dla_new(int nx, int ny, int nz)
{
    struct vol3d *self=NULL;

    if (NULL==(self=vol3df_new_zero(nx,ny,nz)))
        { fprintf(stderr,"dla_new: vol3df_new_zero failed.\n"); return NULL; }



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


/* Choose a random point in the volume. If the voxel at that point is
 * nonzero, then increment it and return. Otherwise, take a random step
 * in one of the Cartesian diractions, and repeat.
 */
int dla_iterate(struct vol3d *self)
{
    int x,y,z;
    long int count=0;

    x = random_whole(self->nx);
    y = random_whole(self->ny);
    z = random_whole(self->nz);

    while ( (!dla_testPoint(self,x,y,z)) && (0!=x) ) {
        count++;
        x += plusminus();
        y += plusminus();
        z += plusminus();

        if ((self->nx-2) <= x) {
            x = random_whole(self->nx);
            y = random_whole(self->ny);
            z = random_whole(self->nz);
        }
        while (x<0) { x += self->nx; }
        while (y<0) { y += self->ny; }
        while (z<0) { z += self->nz; }

        while (x >= self->nx) { x -= self->nx; }
        while (y >= self->ny) { y -= self->ny; }
        while (z >= self->nz) { z -= self->nz; }


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

    for (t=0;t<100000;t++) {
        if (0==(t%iostep)) {
            snprintf(buf,1000,"dla-%05d.h5",t);
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
