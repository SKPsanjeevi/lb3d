#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rpc/rpc.h>

#include "config.h"
#include "vol2df.h"
#include "vol2du8.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

struct vol2d *vol2df_new(int nx, int ny)
{
        return vol2d_new(nx,ny,VOL_TYPE_FLOAT);
}

void vol2df_unsupported_mess(char *featurename)
{
        fprintf(stderr,"%s support was not compiled in.\n",featurename);
}

struct vol2df_stats *vol2df_getstats(struct vol2d *v)
{
        struct vol2df_stats *stat=NULL;
        unsigned long i;
        float phi;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol2df_stats: datatype is not float\n");
                return NULL;
        }

        if (NULL==(stat=malloc(sizeof(struct vol2df_stats)))) {
                perror("vol3df_stats: malloc");
                return NULL;
        }

        stat->max = stat->min = stat->mean = phi = *(float *)v->data;
        stat->variance = phi*phi;

        for (i=1;i<v->nx*v->ny;i++) {
                phi = *( (float *)v->data + i);
                stat->mean += phi;
                stat->variance   += phi*phi;
                if (stat->min>phi) { stat->min=phi; }
                if (stat->max<phi) { stat->max=phi; }
        }

        stat->mean /= (float) (v->nx*v->ny);
        stat->variance /= (float) (v->nx*v->ny);
        stat->variance -= stat->mean*stat->mean;

        return stat;
}

int vol2df_dumpstats(struct vol2d *v, FILE *f)
{
        struct vol2df_stats *stats=NULL;

        if (NULL==(stats=vol2df_getstats(v)))
         {fprintf(stderr,"vol2df_dumpstats: vol2df_stats failed\n"); return -1;}

        fprintf(f,"max %f min %f mean %f variance %f\n",
                        stats->max,
                        stats->min,
                        stats->mean,
                        stats->variance);
        return 0;
}

struct vol2d *vol2df_new_from_2du8(struct vol2d *v)
{
        struct vol2d *vnew=NULL;
        float *ffield=NULL,*out=NULL;
        unsigned char *data=NULL;
        int i;

        if (VOL_TYPE_UINT8 != v->datatype) {
                fprintf(stderr,"vol2df_from_2du8: not a u8 dataset\n");
                return NULL;
        }

        if (NULL==(vnew=vol2df_new(v->nx,v->ny))) {
                perror("vol2df_from_2du8: malloc"); return NULL;
        }

        out=ffield= (float *)vnew->data;
        data = (unsigned char *)v->data;

        for (i=0;i<v->nx*v->ny;i++) {
                *out++ = (float) *data++;
        }

        return vnew;
}

int vol2df_asciidump(struct vol2d *v, FILE *fh)
{
        int x,y;
        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol2df_asciidump: datatype is not float\n");
                return -1;
        }

        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                fprintf(fh,"%+04.04f ",*(float *)vol2d_element(v,x,y));
        }
        fputc('\n',fh);
        }
        return 0;

}

/* NB: corrupts input volume */
struct vol2d *vol2df_new_rfftw(struct vol2d *v)
#ifdef HAVE_FFTW
{
        struct vol2d *vnew=NULL;
        int tnx=0,tny=0;
        fftwf_plan plan;

        tnx=2*(1+(int)(v->nx/2) ); tny=v->ny;

        if (NULL==(vnew=vol2df_new(tnx,tny))) {
                fprintf(stderr,"vol2df_new_rfftw: vol2df_new failed\n");
                return NULL;
        }


        plan = fftwf_plan_dft_r2c_2d (v->ny,v->nx,
                                (float *)v->data, 
                                (fftwf_complex *) vnew->data,
                                FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        return vnew;

}
#else
{
        vol2df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */

struct vol2d *vol2df_new_fullspectrum(struct vol2d *v)
#ifdef HAVE_FFTW
{
        struct vol2d *transform=NULL,*spectrum=NULL;
        int tnx,tny,x,y;
        float *p;

        if (NULL==(transform=vol2df_new_rfftw(v))) {
                fprintf(stderr,
                        "vol2df_new_fullspectrum: vol2df_new_rfftw failed\n");
                return NULL;
        }

        if (NULL==(spectrum=vol2d_new_copytype(v))) {
                fprintf(stderr,
                        "vol2df_new_fullspectrum: vol2d_new_copytype failed\n");
                return NULL;
        }

        tnx=transform->nx/2; tny=transform->ny;
        for (y=0;y<tny;y++) {
        for (x=0;x<tnx;x++) {
                float re,im,mag;
                p=(float *)transform->data + 2*( x + tnx*y);
                re=*p; im=*(p+1); mag=sqrt(re*re+im*im)/(v->nx*v->ny);
                *(float *)vol2d_element_bc(spectrum,x,y)=mag;
                if ((0!=x)&&(0!=y)) {
                *(float *)vol2d_element_bc(spectrum,v->nx-x,v->ny-y)=mag;
                }
        }}
        return spectrum;
}
#else
{
        vol2df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */


/* Check that discrete Parseval theorem obeyed */
int vol2df_rfftw_parsevaltest(struct vol2d *v)
#ifdef HAVE_FFTW
{
        struct vol2d *vnew=NULL;
        int tnx=0,tny=0;
        fftwf_plan plan;
        float suma=0.0,sumb=0.0;
        float *p;
        int i;

        tnx=2*(1+(int)(v->nx/2) ); tny=v->ny;

        if (NULL==(vnew=vol2df_new(tnx,tny))) {
                fprintf(stderr,
                        "vol2df_rfftw_parsevaltest: vol2df_new failed\n");
                return -1;
        }

        p=v->data;
        for (i=0;i<v->nx*v->ny;i++) {
                float phi;
                phi = *p++;
                suma += phi*phi;
        }

        plan = fftwf_plan_dft_r2c_2d (v->ny,v->nx,
                                (float *)v->data, 
                                (fftwf_complex *) vnew->data,
                                FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        {
                int foox,fooy,x,y;
                float re,im;
                float weight;
                foox = tnx/2; fooy=tny;

                if (0==(v->nx%2)) { /* nx is even */
                        for (y=0;y<fooy;y++) {
                                for (x=0;x<foox;x++) {
                                        p=(float *)vnew->data + 2*( x + foox*y);
                                        re = *p;
                                        im = *(p+1);
                                        if ((0==x)||((foox-1)==x)) {
                                                weight=1.0;
                                        } else { weight=2.0; }

                                        sumb += weight *( re*re+im*im); 
                                }
                        }
                } else { /* nx is odd */
                        for (y=0;y<fooy;y++) {
                                for (x=0;x<foox;x++) {
                                        p=(float *)vnew->data + 2*( x + foox*y);
                                        re = *p;
                                        im = *(p+1);
                                        if (0==x) {
                                                weight=1.0;
                                        } else { weight=2.0; }

                                        sumb += weight *( re*re+im*im); 
                                }
                        }
                }
        printf("nx=%d ny=%d\n",v->nx,v->ny);
        printf("foox=%d fooy=%d\n",foox,fooy);

        }

        printf("suma=%f sumb=%f mean=%f\n",suma,sumb,sumb/(v->nx*v->ny));

        return 0;


}
#else
{
        vol2df_unsupported_mess("fftw");
        return -1;
}
#endif /* HAVE_FFTW */


static void vol2df_subtract_streamfunc(void *el, void *data)
{ *(float *)el -= *(float *)data; }

/* Subtract mean value from each voxel, so that overall mean is zero */
int vol2df_zeromean(struct vol2d *v)
{
        struct vol2df_stats *stats=NULL;

        if (NULL==(stats=vol2df_getstats(v))) {
                fprintf(stderr,"vol2df_zeromean: vol2df_getstats failed\n");
                return -1; 
        }

        vol2d_stream(v,vol2df_subtract_streamfunc,&stats->mean);

        free(stats);
        return 0;
}

#ifdef HAVE_FFTW
int is_mirrored_2d(struct vol2d *v, int x)
{
        if (0==(v->nx%2)) {
                if ( (0==x)||(x==(v->nx/2))) { return 0; }
                else { return 1; }
        } else {
                if (0==x) { return 0; } else { return 1; }
        }
}
#endif /* HAVE_FFTW */


struct binnery *vol2df_sfactor(struct vol2d *v)
#ifdef HAVE_FFTW
{
        struct vol2d *transform=NULL;
        int tnx=0,tny=0;
        int i,j,ii,jj;
        int n,nbins;
        int r;
        struct binnery *binnery=NULL;

        if (VOL_TYPE_FLOAT != v->datatype)
        { fprintf(stderr,"vol2df_sfactor: not a float dataset\n");return NULL;}

        if (NULL==(transform=vol2df_new_rfftw(v)))
        { fprintf(stderr,"vol2df_sfactor: new_rfftw failed\n"); return NULL; }

        n=v->nx*v->ny;

        tnx=transform->nx/2; tny=transform->ny;

        /* Set nbins equal to half of the smaller of the system dimensions.  */
        if (v->nx<v->ny) { nbins=v->nx/2; } else { nbins=v->ny/2; }

        if (NULL==(binnery=binnery_new(nbins)))
        { fprintf(stderr,"vol2df_sfactor: binnery_new failed\n");return NULL;}

        for (j=0;j<tny;j++) {
        for (i=0;i<tnx;i++) {
                float re,im,s;

                re = *( (float *)transform->data + 2*(i + j*tnx)   );
                im = *( (float *)transform->data + 2*(i + j*tnx) +1);
                s = (re*re+im*im);

                ii=i;
                if (j>(v->ny/2)) { jj=j-v->ny; } else { jj=j; }
                /* ii>0; jj may be pos or neg; (ii,jj) corresponds to the
                 * unwrapped frequency vector. See notes of 19/5/2004 in
                 * lab book.
                 */

                r = (int)floor(sqrt(ii*ii+jj*jj));
                binnery_add_point(binnery,s,r);

                /* If this point is the complex conjugate of another,
                 * redundant, point, then sum the conjugate point
                 * as well.
                 */
                if (is_mirrored_2d(v,i)) {
                        ii = v->nx - i;
                        jj = v->ny - j;
                        if (jj > (v->ny/2)) { jj -= v->ny; }
                        r = (int)floor(sqrt(ii*ii+jj*jj));
                        binnery_add_point(binnery,s,r);
                }
        }}

        vol2d_destroy(transform);

        return binnery;


        
}
#else
{
        vol2df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */


unsigned int *vol2df_sfactor_old(struct vol2d *v)
#ifdef HAVE_FFTW
{
        struct vol2d *transform=NULL;
        int tnx=0,tny=0;
        int i,j,ii,jj;
        int n,nbins;
        unsigned int *count=NULL;
        float *bin=NULL;
        int r;
        int nvisited=0;

        if (VOL_TYPE_FLOAT != v->datatype)
        { fprintf(stderr,"vol2df_sfactor: not a float dataset\n");return NULL;}

        if (NULL==(transform=vol2df_new_rfftw(v)))
        { fprintf(stderr,"vol2df_sfactor: new_rfftw failed\n"); return NULL; }

        n=v->nx*v->ny;

        tnx=transform->nx/2; tny=transform->ny;

        if (v->nx<v->ny) { nbins=v->nx; } else { nbins=v->ny; }

        if (NULL==(bin=malloc(sizeof(float)*nbins)))
                { perror("vol2df_sfactor: malloc"); return NULL; }
        if (NULL==(count=malloc(sizeof(unsigned int)*nbins)))
                { perror("vol2df_sfactor: malloc"); return NULL; }
        for (i=0;i<nbins;i++) { count[i]=0; bin[i]=0.0; }

        for (j=0;j<tny;j++) {
        for (i=0;i<tnx;i++) {
                float re,im,s;

                re = *( (float *)transform->data + 2*(i + j*tnx)   );
                im = *( (float *)transform->data + 2*(i + j*tnx) +1);
                s = re*re+im*im/n;

                ii=i;
                if (j>(v->ny/2)) { jj=j-v->ny; } else { jj=j; }
                /* ii>0; jj may be pos or neg; (ii,jj) corresponds to the
                 * unwrapped frequency vector. See notes of 19/5/2004 in
                 * lab book.
                 */

                r = (int)floor(sqrt(ii*ii+jj*jj));
                nvisited++;
                if (r<nbins) { bin[r] += s; count[r]++; }

                /* If this point is the complex conjugate of another,
                 * redundant, point, then sum the conjugate point
                 * as well.
                 */
                if (is_mirrored_2d(v,i)) {
                        nvisited++;
                        ii = v->nx - i;
                        jj = v->ny - j;
                        if (jj > (v->ny/2)) { jj -= v->ny; }
                        r = (int)floor(sqrt(ii*ii+jj*jj));
                        if (r<nbins) { bin[r] += s; count[r]++; }
                }
        }}

        printf("nvisited=%d\n",nvisited);

        for (i=0;i<nbins;i++) { printf("%d %f\n",i,bin[i]/count[i]); }

        return NULL;


        
}
#else
{
        vol2df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */



static void accum_streamfunc(void *el,void *data)
{

        struct vol2df_dist *dist;
        int i;
        float phi;

        phi = *(float *)el;
        dist = (struct vol2df_dist *)data;
        i=(int)floor( ( phi - dist->stats->min)/dist->dphi );

        if (i<0) {
                fprintf(stderr,"Warning: phi=%f < phimin=%f\n",
                        phi,dist->stats->min);
                i=0;
        }
        if (i>=dist->nbins) {
                fprintf(stderr,"Warning: phi=%f > phimax=%f\n",
                        phi,dist->stats->max);
                i=dist->nbins-1;
        }

        dist->count[i]++;
}

struct vol2df_dist *vol2df_getdist(struct vol2d *v, int nbins)
{

        struct vol2df_dist *dist=NULL;
        int i;


        if (NULL==(dist=malloc(sizeof(struct vol2df_dist)))) {
                perror("vol2df_getdist: malloc");
                return NULL;
        }

        dist->nbins = nbins;

        if (NULL==(dist->stats=vol2df_getstats(v))) {
                fprintf(stderr,"vol2df_dist: vol2df_getstats failed\n");
                return NULL;
        }

        /* Allocate bins */

        if (NULL==(dist->count=malloc(sizeof(unsigned int)*nbins))) {
                perror("malloc()");
                return NULL;
        }

        /* Zero all counts */

        for (i=0;i<dist->nbins;i++) {
                dist->count[i]=0;
        }

        dist->dphi = (dist->stats->max - dist->stats->min)/(float)(nbins-1);

        vol2d_stream(v,accum_streamfunc,dist);

        return dist;

}

int vol2df_dumpdist(struct vol2d *v, int nbins, FILE *f)
{
        struct vol2df_dist *dist=NULL;
        int i;

        if (VOL_TYPE_FLOAT != v->datatype) {
                fprintf(stderr,"vol2df_dumpdist: volume is not float\n");
                return -1;
        }

        if (NULL==(dist=vol2df_getdist(v,nbins))) {
                fprintf(stderr,"vol2df_dumpdist: vol2df_getdist failed\n");
                return -1;
        }

        for (i=0;i<nbins;i++) {
                fprintf(f,"%f %d\n",
                        dist->stats->min+dist->dphi*(float)i,
                        dist->count[i]
                );
        }

        free(dist);

        return 0;
}

#ifdef HAVE_PNG
int vol2df_write_png(struct vol2d *v, char *filename, float pmax, float pmin)
{
        struct vol2d *v_u8=NULL;

        if (NULL==(v_u8=vol2du8_new_from_2df(v,pmax,pmin))) {
                fprintf(stderr,
                        "vol2df_write_png: vol2du8_new_from_2df failed.\n");
                return -1;
        }

        if (0!=vol2du8_write_png(v_u8,filename)) {
                fprintf(stderr,"vol2df_write_png: vol2du8_write_png failed.\n");
                vol2d_destroy(v_u8);
                return -1;
        }
        vol2d_destroy(v_u8);
        return 0;
}

int vol2df_write_png_normalize(struct vol2d *v, char *filename)
{
        struct vol2d *v_u8=NULL;

        if (NULL==(v_u8=vol2du8_new_from_2df_normalize(v))) {
                fprintf(stderr,
                "vol2df_write_png: vol2du8_new_from_2df_normalize failed.\n");
                return -1;
        }

        if (0!=vol2du8_write_png(v_u8,filename)) {
                fprintf(stderr,"vol2df_write_png: vol2du8_write_png failed.\n");
                vol2d_destroy(v_u8);
                return -1;
        }
        vol2d_destroy(v_u8);
        return 0;
}

#endif /* HAVE_PNG */

void vol2df_xdr_streamfunc(void *el,void *data)
{

        xdr_float((XDR *)data,(float *)el);

}

void vol2df_xdrd_streamfunc(void *el,void *data)
{
        double phi;

        xdr_double((XDR *)data,&phi);
        *(float *)el = (float) phi;

}

struct vol2d *vol2df_new_from_xdrfh (
                unsigned int nx, unsigned int ny,FILE *f)
{
        struct vol2d *v=NULL;
        XDR xdrs;

        if (NULL==(v=vol2df_new(nx,ny))) {
                fprintf(stderr,"vol2df_new_from_xdrfh: vol2df_new failed\n");
                return NULL;
        }

        xdrstdio_create(&xdrs,f,XDR_DECODE);

        vol2d_stream(v,vol2df_xdr_streamfunc,&xdrs);

        xdr_destroy(&xdrs);

        return v;

}
struct vol2d *vol2df_new_from_xdrdfh (
                unsigned int nx, unsigned int ny,FILE *f)
{
        struct vol2d *v=NULL;
        XDR xdrs;

        if (NULL==(v=vol2df_new(nx,ny))) {
                fprintf(stderr,"vol2df_new_from_xdrfh: vol2df_new failed\n");
                return NULL;
        }

        xdrstdio_create(&xdrs,f,XDR_DECODE);

        vol2d_stream(v,vol2df_xdrd_streamfunc,&xdrs);

        xdr_destroy(&xdrs);

        return v;

}
