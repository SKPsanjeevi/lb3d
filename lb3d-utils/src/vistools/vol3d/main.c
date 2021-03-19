#include <stdio.h>
#include <stdlib.h>

#include <netinet/in.h>

#include "vol3d.h"

/* Transfer function lookup tables */
uint8_t r[256],g[256],b[256],a[256];

/* Convert unsigned byte input data to
 * network order unsigned byte (rgba)-tuple
 */
void transfer (void *src, void *dst, void *data)
{
	uint8_t r,g,b,a;
	uint16_t data_in;
	uint32_t hostlong;
	uint32_t data_out;

	data_in = *(uint16_t *) src;
	r=g=b=a=(uint8_t) (data_in >> 8);

	if ( a>64) {
		r=g=b=a=200;
	} else {
		r=g=b=a=0;
	}

	hostlong = r & 0xff;
	hostlong <<= 8;
	hostlong |= (g & 0xff);
	hostlong <<= 8;
	hostlong |= (b & 0xff);
	hostlong <<= 8;
	hostlong |= (a & 0xff);

	data_out = htonl(hostlong);
	*(uint32_t *)dst = (uint32_t) data_out;

}


void init_tables(void)
{
	int i;
	for (i=0;i<256;i++) {
		if (i<70) {
			r[i] = g[i]=b[i]=a[i]=0;
		} else if (i<100) {
			r[i]=a[i]=255;
			g[i]=b[i]=0;
		} else if (i<150) {
			g[i]=a[i]=255;
			r[i]=b[i]=0;
		} else  {
			b[i]=a[i]=255;
			g[i]=r[i]=0;
		}

	}
}

void transfer_16to8(void *src,void *dst, void *data)
{
	uint16_t in;

	in = *(uint16_t *)src;
	*(uint8_t *)dst = (uint8_t)(in>>8);
}

int main(int argc, char *argv[])
{
	char infilename[6+1+3+1]; /* 'CThead' + '.' + 'nnn' + \0 */
	char *outfilename="cthead.rgba";

	struct vol3d *v=NULL,*v8bit=NULL,*vrgba=NULL;
	unsigned int nx=256,ny=256,nz=128;
	unsigned int nzslices=113;
	uint16_t zero=0;
	unsigned int z;

	init_tables();

	if (NULL==(v=vol3d_new(nx,ny,nz,VOL3D_UINT16))) {
		fprintf(stderr,"vol3d_new failed\n");
		return -1;
	}

	for (z=0;z<nzslices;z++) {
		snprintf(infilename,sizeof(infilename),
				"CThead.%d",z+1);
		if (0!=vol3d_readzslice_raw(v,z,infilename)) {
			fprintf(stderr,"vol3d_readzslice_raw failed on z=%d\n",
					z);
			vol3d_destroy(v);
			return -1;
		}
		printf("Read <%s>\n",infilename);fflush(stderr);
	}
	/* Fill remaining volume with zeroes */

	if (0!=vol3d_set_subvol(v,
				0,0,nzslices,
				nx,ny,nz-nzslices,
				&zero)) {
		fprintf(stderr,"vol3d_set_subvol failed\n");
		vol3d_destroy(v);
	}

	/* Convert to 8-bits */

	if (NULL==(v8bit=
		vol3d_new_transfer_elsize(v,sizeof(uint8_t),transfer_16to8,NULL))) {
		fprintf(stderr,"vol3d_new_transfer failed\n");
		vol3d_destroy(v);
		return -1;
	}
	vol3d_destroy(v);

	/* v now contains the entire scalar dataset. */

	if (NULL==(vrgba=vol3d_new_uc2rgba_lookup(v8bit,r,g,b,a))) {
		fprintf(stderr,"vol3d_new_uc2rgba_lookup failed\n");
		return -1;
	}

	if (0!=vol3d_write_raw(vrgba,outfilename)) {
		fprintf(stderr,"vol3d_write_raw failed\n");
		vol3d_destroy(v8bit);
		vol3d_destroy(vrgba);
		return -1;
	}

	vol3d_destroy(v8bit);
	vol3d_destroy(vrgba);
	return 0;
}
