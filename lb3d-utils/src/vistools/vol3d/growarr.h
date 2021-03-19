#ifndef INCLUDED_GROWARR_H
#define INCLUDED_GROWARR_H

struct growarr {
        size_t itemsize;        /* Size of one item */
        unsigned int nallocated;         /* Maximum size array can hold */
        unsigned int nitems;             /* Number currently in array */
        void *arr;              /* The array */
};
const char *growarr_cvsid(void) ;
struct growarr *growarr_new(size_t itemsize, unsigned int nitems);
int growarr_destroy(struct growarr *ga);
int growarr_expand(struct growarr *ga, unsigned int nitems);
int growarr_add(struct growarr *ga, void *item);
int growarr_shrink(struct growarr *ga, unsigned int nitems);
void *growarr_get_nth(struct growarr *ga, unsigned int n);
unsigned int growarr_nitems(struct growarr *ga);


#endif /* INCLUDED_GROWARR_H */
