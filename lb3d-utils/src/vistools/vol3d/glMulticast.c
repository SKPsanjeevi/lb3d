#include <flx.h>

/* Write a simple image to the multicast stream */

int glMulticastFrame(int nx, int ny,flx_t *flx_tx )
{
	unsigned char *buffer;
	unsigned char *b;
	double *buf, *p;
        int i;
	double phi, phimin, phimax;

        if (NULL==(buffer=malloc(3*nx*ny))) {
                perror("sim_multicast_frame failed to allocate buffer");
                return -1;
        }

        /* Read the colour field into a buffer of size nx x ny */

        glReadPixels(0,0,nx,ny, GL_RGB,
                GL_UNSIGNED_BYTE,(GLvoid *)buffer);

	flx_image( flx_tx, buffer );
	free( buffer );
	return 0;
} 
