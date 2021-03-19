#include <stdio.h>
#include <stdlib.h>
#include "growarr.h"


/*** small test suite */

struct foo {
        char str[16];
        int bar;
};

int test_new_destroy(void)
{
        struct growarr *ga=NULL;

        if (NULL==(ga=growarr_new(sizeof(struct foo),12))) {
                fprintf(stderr,"growarr_new failed\n");
                return -1;
        }

        if (0!=growarr_destroy(ga)) {
                fprintf(stderr,"growarr_destroy failed\n");
                return -1;
        }
        return 0;
}

int test_allocation(void)
{
        struct growarr *ga=NULL;
        struct foo foo;
        int copy;
        const unsigned int nit = 13;
        int i,j;

        if (NULL==(ga=growarr_new(sizeof(struct foo),nit))) {
                fprintf(stderr,"growarr_new failed\n");
                return -1;
        }

        for (i=0;i<nit;i++) {
                for (j=0;j<16;j++) { foo.str[j] = 'A'+j ; }
                foo.bar = i;

                if (0>(copy=growarr_add(ga,&foo))) {
                        fprintf(stderr,"growable_add failed on %d\n",i);
                        return -1;
                }
        }

        if (0!=growarr_destroy(ga)) {
                fprintf(stderr,"growarr_destroy failed\n");
                return -1;
        }

        return 0;
}

int test_expansion(void)
{
        struct growarr *ga=NULL;
        struct foo foo;
        struct foo *copy;
        const unsigned int nit = 128;
        int i,j;

        if (NULL==(ga=growarr_new(sizeof(struct foo),1))) {
                fprintf(stderr,"growarr_new failed\n");
                return -1;
        }

        for (i=0;i<nit;i++) {
                for (j=0;j<16;j++) { foo.str[j] = 'A'+j ; }
                foo.bar = i;

                if (0>growarr_add(ga,&foo)) {
                        fprintf(stderr,"growable_add failed on %d\n",i);
                        return -1;
                }
        }

        for (i=0;i<nit;i++) {
                if (NULL==(copy=growarr_get_nth(ga,i))) {
                        fprintf(stderr,"growarr_get_nth failed on %d\n",i);
                        return -1;
                }

                if (i != copy->bar) {
                        fprintf(stderr,"Bad element %d=%d\n",i,copy->bar);
                        return -1;
                }
        }

        if (NULL!=growarr_get_nth(ga,nit)) {
                fprintf(stderr,"Returned non-NULL for (nitems)th element!\n");
                return -1;
        }

        if (0!=growarr_destroy(ga)) {
                fprintf(stderr,"growarr_destroy failed\n");
                return -1;
        }

        return 0;

}

int main(int argc, char *argv[])
{
        printf("Testing constructor and destructor.."); fflush(stdout);
        if (0!=test_new_destroy()) {printf("FAILED\n");} else { printf("ok\n");}
        printf("Testing allocation.."); fflush(stdout);
        if (0!=test_allocation()) {printf("FAILED\n");} else { printf("ok\n");}
        printf("Testing expansion.."); fflush(stdout);
        if (0!=test_expansion()) {printf("FAILED\n");} else { printf("ok\n");}
	return 0;
}
