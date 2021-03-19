#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "growarr.h"

/* Code for growable array of arbitrary data.
 * realloc() may be called on the data, so pointers into the array
 * should not be kept around: store indices and use growarr_get_nth instead.
 */

static const char *cvsid =
"$Id: growarr.c,v 1.2 2004/09/03 16:57:48 jon Exp $";

const char *growarr_cvsid(void) {
        return cvsid;
}


/* Create a new growable array to hold items of given size,
 * and at least nitems without reallocating. Return pointer on success,
 * NULL indicates failure.
 */

struct growarr *growarr_new (size_t itemsize, unsigned int nitems)
{
        struct growarr *ga;

        if (NULL==(ga=malloc(sizeof(struct growarr))))
                { perror("growarr_new: malloc"); return NULL; }

        if (NULL==(ga->arr=malloc(itemsize*nitems)))
                { perror("growarr_new: malloc"); free(ga); return NULL; }

        ga->itemsize = itemsize;
        ga->nallocated = nitems;
        ga->nitems = 0;
        
        return ga;
}

/* Deallocate a growarr and all its data. Return 0 on success, else
 * failure 
 */
int growarr_destroy(struct growarr *ga) {
        if (NULL==ga) { fprintf(stderr,"growarr_destroy(NULL)\n"); return -1; }
        if (NULL==ga->arr)
                { fprintf(stderr,"growarr_destroy(arr=NULL)\n"); return -1; }
        free(ga->arr);
        free(ga);
        return 0;
}

/* Expand the array to given size. Return 0 on success, nonzero on failure
 * If nitems==0, then will just expand, at present doubles size.
 */

int growarr_expand(struct growarr *ga, unsigned int nitems)
{
        void *newptr=NULL;

        if (0==nitems) { nitems = ga->nallocated*2; }
        if (nitems==ga->nallocated) { return 0; }
        if (nitems<ga->nallocated) {
                fprintf(stderr,"growarr_expand: attempt to shrink!\n");
                return -1;
        }

        if (NULL==(newptr=realloc(ga->arr,nitems*ga->itemsize)))
                { perror("growarr_expand: realloc"); return -1; }

        ga->arr = newptr;
        ga->nallocated = nitems;

        return 0;
}

/* Copy given item into the growable array, reallocating
 * if necessary. Return the index of the item in the array on success,
 * negative on failure.
 */

int growarr_add(struct growarr *ga, void *item)
{
        void *destptr=NULL;

        if (NULL==ga) { fprintf(stderr,"growarr_add(NULL)\n");return -1;}

        if (ga->nitems>=ga->nallocated) {
                if (ga->nitems>ga->nallocated) {
                        fprintf(stderr,
                               "WARNING: growarr_add: items=%d allocated=%d\n",
                                        ga->nitems,ga->nallocated);
                }

                /* List too small, so grow it */

                if (0!=growarr_expand(ga,2*ga->nallocated)) {
                        fprintf(stderr,
                                "growarr_add: growarr_expand failed\n");
                        return -1;
                }
        }

        /* Enough entries free */

        if (ga->nitems>=ga->nallocated) {
                /* Can't happen */
                fprintf(stderr,"growarr_add: too small after expansion!\n");
                return -1;
        }


        destptr = (unsigned char *) ga->arr
                + ga->itemsize * ga->nitems;

        destptr = memcpy(destptr,item,ga->itemsize);

        ga->nitems++;
        return ga->nitems-1;
}

/* Shrink the array to nitems, or to its smallest possible size if
 * nitems is zero.
 * Return zero on success, nonzero on failure.
 */

int growarr_shrink(struct growarr *ga, unsigned int nitems)
{
        void *newptr;

        if (0==nitems) { nitems = ga->nallocated; }

        if (NULL==(newptr=realloc(ga->arr,nitems*ga->itemsize)))
                { perror("growarr_shrink: realloc"); return -1; }

        ga->nallocated = nitems;
        
        return 0;
}



/* Return pointer to nth element. n is zero-based. */
void *growarr_get_nth(struct growarr *ga,unsigned int n)
{
        if (n>=ga->nitems) {
                fprintf(stderr,"growarr_get_nth(%u) in list of %u!\n",
                                n,ga->nitems);
                fflush(stderr);
                return NULL;
        } /* Out of bounds */
        return (unsigned char *)ga->arr + n*ga->itemsize;
}


/* Return number of items in array.
 * This is the same as one plus the index of the last item.
 */
unsigned int growarr_nitems(struct growarr *ga) {
        if (NULL==ga) { fprintf(stderr,"growarr_nitems(NULL)\n"); return -1; }
        return ga->nitems;
}

