#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "binnery.h"

struct binnery *binnery_new(int nbins)
{
        struct binnery *b=NULL;
        int i;
        
        if (NULL==(b=malloc(sizeof(struct binnery))))
        { perror("binnery_new: malloc"); return NULL; }

        if (NULL==(b->count=malloc(sizeof(unsigned int)*nbins)))
        { perror("binnery_new: malloc"); binnery_destroy(b); return NULL; }

        if (NULL==(b->mean=malloc(sizeof(unsigned int)*nbins)))
        { perror("binnery_new: malloc"); binnery_destroy(b); return NULL; }

        if (NULL==(b->variance=malloc(sizeof(unsigned int)*nbins)))
        { perror("binnery_new: malloc"); binnery_destroy(b); return NULL; }

        b->nbins=nbins;

        for (i=0;i<nbins;i++) { b->mean[i]= b->variance[i]=0.0; b->count[i]=0; }

        return b;
}


int binnery_destroy(struct binnery *b)
{
        if (NULL==b) {
                fprintf(stderr,"binnery_destroy: NULL binnery\n"); return -1;
        }

        if (NULL!=b->count) { free(b->count); }
        if (NULL!=b->mean) { free(b->mean); }
        if (NULL!=b->variance) { free(b->variance); }

        free(b);
        return 0;
}

int binnery_add_point(struct binnery *b, float x, unsigned int i)
{
        float mprev,sprev;
        float k;

        if (i>=b->nbins) { return -1; }

        /* Use recurrence formula to ensure that variance is always
         * positive. See Knuth II p232.
         */

        b->count[i]++;
        k=b->count[i];

        mprev = b->mean[i]; sprev = b->variance[i];

        b->mean[i] = mprev + (x-mprev)/k;
        b->variance[i] += (x-mprev)*(x-b->mean[i]);
        return 0;
}

float binnery_mean(struct binnery *b, unsigned int i)
{
        return b->mean[i];
}

float binnery_variance(struct binnery *b, unsigned int i)
{
        return b->variance[i];
}

float binnery_sd(struct binnery *b, unsigned int i)
{
        if (b->count[i]>0) {
                return sqrt(b->variance[i]/b->count[i]);
        } else { return 0; }
}

float binnery_1stmoment(struct binnery *b)
{
        float sum=0.0;
        int i;


        for (i=0;i<b->nbins;i++) {
                sum += b->mean[i] * ((float)i+0.5);
        }
        return sum;
}
float binnery_sum(struct binnery *b)
{
        float sum=0.0;
        int i;

        for (i=0;i<b->nbins;i++) {
                sum += b->mean[i];
        }
        return sum;
}

int binnery_dump(struct binnery *b, FILE *f)
{
        int i;

        for (i=0;i<b->nbins;i++) {
                fprintf(f,"%d %.10f %.10f\n",i,
                                binnery_mean(b,i),
                                binnery_sd(b,i)
                       );
        }
        return 0;
}

int binnery_scaled_dump(struct binnery *b, FILE *f,float scale)
{
        int i;

        for (i=0;i<b->nbins;i++) {
                fprintf(f,"%.10f %.10f %.10f\n",scale*(i+0.5),
                                binnery_mean(b,i),
                                binnery_sd(b,i)
                       );
        }
        return 0;
}

void binnery_scaled_1stmoment(struct binnery *b, float scale,float *mom,float *err)
{
        float sums=0.0,sumsk=0.0;
        float variance_sums=0.0,variance_sumsk=0.0;
        float varmom;
        int i;

        for (i=0;i<b->nbins;i++) {
                float x;
                x=scale*(0.5+i);

                sums += binnery_mean(b,i);
                variance_sums += binnery_variance(b,i);
                sumsk += binnery_mean(b,i)*x;
                variance_sumsk += x*x*binnery_variance(b,i);

        }

        *mom = sumsk/sums;

        varmom = (*mom)*(*mom) * ( variance_sumsk/(sumsk*sumsk)
                                 + variance_sums/(sums*sums)
                        );

        *err = sqrt(varmom);
                           

        return;

}
