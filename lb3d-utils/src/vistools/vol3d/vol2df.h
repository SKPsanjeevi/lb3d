#ifndef INCLUDED_VOL2DF_H
#define INCLUDED_VOL2DF_H

#include "vol2d.h"
#include "binnery.h"

struct vol2df_stats {
        float max;
        float min;
        float mean;
        float variance;
};

struct vol2df_dist {
        unsigned int nbins;
        unsigned int *count;
        float dphi;
        struct vol2df_stats *stats;
};


struct vol2d *vol2df_new(int nx, int ny);
struct vol2df_stats *vol2df_getstats(struct vol2d *v);
struct vol2d *vol2df_new_from_2du8(struct vol2d *v);
int vol2df_asciidump(struct vol2d *v, FILE *fh);
struct vol2d *vol2df_new_rfftw(struct vol2d *v);
int vol2df_dumpstats(struct vol2d *v, FILE *f);
struct vol2df_dist *vol2df_getdist(struct vol2d *v, int nbins);
int vol2df_dumpdist(struct vol2d *v, int nbins, FILE *f);
struct vol2d *vol2df_new_fullspectrum(struct vol2d *v);
int vol2df_zeromean(struct vol2d *v);
unsigned int *vol2df_sfactor_old(struct vol2d *v);
struct binnery *vol2df_sfactor(struct vol2d *v);

#ifdef HAVE_PNG
int vol2df_write_png_normalize(struct vol2d *v, char *filename);
int vol2df_write_png(struct vol2d *v, char *filename, float pmax, float pmin);
#endif /* HAVE_PNG */
void vol2df_xdr_streamfunc(void *el,void *data);
void vol2df_xdrd_streamfunc(void *el,void *data);
struct vol2d *vol2df_new_from_xdrfh (
                unsigned int nx, unsigned int ny,FILE *f);
struct vol2d *vol2df_new_from_xdrdfh (
                unsigned int nx, unsigned int ny,FILE *f);
#endif /* INCLUDED_VOL2DF_H */
