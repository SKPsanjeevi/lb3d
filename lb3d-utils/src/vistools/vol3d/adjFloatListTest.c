#include <stdio.h>
#include <stdlib.h>

#include "adjFloatList.h"

#define VERBOSE(x) \
    x \
    fprintf(stdout,#x "\n");

int main(int argc, char *argv[])
{
    struct adjFloatList *fl=NULL;
    GLfloat a=1,b=2,c=3;

    if (NULL==(fl=adjFloatList_new())) {
        fprintf(stderr,"adjFloatList_new failed.\n");
        return -1;
    }

    adjFloatList_dump(fl,stdout);
    adjFloatList_destroy(fl); /* Destroy empty list */

    if (NULL==(fl=adjFloatList_new())) {
        fprintf(stderr,"adjFloatList_new failed.\n");
        return -1;
    }

    VERBOSE(adjFloatList_wrapOn(fl);)
    VERBOSE(adjFloatList_addFloat(fl,&a,10,0,1,"a");)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_addFloat(fl,&b,10,0,1,"b");)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_addFloat(fl,&c,10,0,1,"c");)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_next(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_incFloat(fl);)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_incFloat(fl);)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_decFloat(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_next(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_next(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_next(fl);)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_prev(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_prev(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_wrapOff(fl);)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_next(fl);)
    adjFloatList_dump(fl,stdout);
    VERBOSE(adjFloatList_prev(fl);)
    VERBOSE(adjFloatList_prev(fl);)
    VERBOSE(adjFloatList_prev(fl);)
    adjFloatList_dump(fl,stdout);

    VERBOSE(adjFloatList_toggleWrap(fl);)
    adjFloatList_dump(fl,stdout);

    adjFloatList_destroy(fl);

    fprintf(stdout,"All done!\n");

    return 0;
}
