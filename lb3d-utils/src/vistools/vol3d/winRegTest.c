#include <stdio.h>
#include <stdlib.h>
#include "winReg.h"

struct winRegEnt *winReg_getHead(void);

int winReg_dumpStrings(void)
{
    struct winRegEnt *ent=NULL;

    ent = winReg_getHead();

    while (NULL!=ent) {
        printf("win=%03d ptr=%p str=\"%s\"\n",ent->win,ent->ptr,
                (char *)ent->ptr);
        ent=ent->next;
    }

    return 0;
}

int main(int argc, char *argv[])
{
    winReg_setPtr(1,"foo");
    winReg_setPtr(23,"bar");
    winReg_setPtr(456,"knobs");
    winReg_setPtr(8,"bananana");
    winReg_dumpStrings();

    printf("remove(123)\n");fflush(stderr);fflush(stdout);
    winReg_remove(123);
    winReg_dumpStrings();

    printf("remove(23)\n");fflush(stderr);fflush(stdout);
    winReg_remove(23);
    winReg_dumpStrings();

    printf("set(1,newfoo)\n");fflush(stderr);fflush(stdout);
    winReg_setPtr(1,"newfoo");
    winReg_dumpStrings();

    printf("set(23,newBar)\n");fflush(stderr);fflush(stdout);
    winReg_setPtr(23,"newBar");
    winReg_dumpStrings();

    {
        int tries[]={1,456,23,8,999,-1};
        int i=0;

        while (tries[i]>0) {
                printf("get(%d)=%s\n",tries[i],(char *)winReg_getPtr(tries[i]));
                i++;
        }
    }


    return 0;
}
