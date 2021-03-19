#ifndef INCLUDED_BINNERY_H
#define INCLUDED_BINNERY_H


struct binnery {
        unsigned int nbins;
        unsigned int *count;
        float *mean;
        float *variance;
        float delta;
};

struct binnery *binnery_new(int nbins);
int binnery_destroy(struct binnery *b);
float binnery_1stmoment(struct binnery *b);
float binnery_sum(struct binnery *b);
int binnery_add_point(struct binnery *b, float phi, unsigned int n);
float binnery_mean(struct binnery *b, unsigned int i);
float binnery_variance(struct binnery *b, unsigned int i);
float binnery_sd(struct binnery *b, unsigned int i);
int binnery_dump(struct binnery *b, FILE *f);
int binnery_scaled_dump(struct binnery *b, FILE *f,float scale);
void binnery_scaled_1stmoment(struct binnery *b, float scale,float *mom,float *err);
#endif /* INCLUDED_BINNERY_H */
