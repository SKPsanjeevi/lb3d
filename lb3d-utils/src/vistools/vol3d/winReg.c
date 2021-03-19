#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "winReg.h"

/*
 * When programming callbacks in GLUT, it's useful for the callback to
 * be able to access state data for the corresponding window.
 * Unfortunately GLUT will not let you pass information through to a
 * callback in the way that, eg, GTK does. However, it is possible to
 * find the window ID using glutGetWindow(), so you can associate
 * information with that.
 *
 * This code provides a window registry to make callback programming
 * easier. When setting up a window's callbacks, winReg_add(void *ptr)
 * will associate the arbitrary pointer ptr with that window. In the
 * callback, winReg_get() will retrieve the pointer for the
 * corresponding window.
 */


static struct winRegEnt *winRegHead = NULL;

struct winRegEnt *winReg_getHead(void)
{
    return winRegHead;
}

/* Ensure that "win" is a valid window */
static int winReg_getWindow(int win)
{
    if ((0==win)&&(0==(win=glutGetWindow())))
        { fprintf(stderr,"winReg_getWindow: glutGetWindow()==0"); return -1; }
    return win;
}

/* Find the entry with the corresponding id. If found, return
 * a pointer to the "next" pointer to this entry. The "next" pointer
 * is either the "next" element of the previous element, or the head
 * pointer.
 *
 * Return NULL if not found.
 */

static struct winRegEnt **winReg_find(int win)
{
    struct winRegEnt **node;

    win = winReg_getWindow(win);
    node=&winRegHead; /* never NULL */

    while (NULL != *node) {
        if (win == (*node)->win) {
            return node;
        }
        node = &( (*node)->next );
    }
    return NULL;
}


/* If win==0, use current window ID. Set the pointer for the
 * corresponding ID, overwriting if necessary.
 */

int winReg_setPtr(int win, void *ptr)
{
    struct winRegEnt **node;

    win = winReg_getWindow(win);

    if (NULL==(node=winReg_find(win))) {

        /* This window has no entry; create one and insert it
         * at the list head.
         */

        struct winRegEnt *newent=NULL;

        if (NULL==(newent=malloc(sizeof(struct winRegEnt))))
            { perror("winReg_setPtr: malloc"); return -1; }

        newent->ptr = ptr;
        newent->win = win;
        newent->next = winRegHead;
        winRegHead = newent;

        return 0;
    }

    /* node points to a pointer to the entry */

    if (win != (*node)->win) {
        /* This should never happen. */
        fprintf(stderr,"winReg_setPtr: INCONSISTENT RESULT; exiting.\n");
        exit(-1);
    }

    (*node)->ptr = ptr;

    return 0;
}

/* Destroy the entry for the given window ID, or current if win==0 */

int winReg_remove(int win)
{
    struct winRegEnt **node,*ent;
    win = winReg_getWindow(win);

    if (NULL==(node=winReg_find(win))) {
        fprintf(stderr,"winReg_remove: no such window %d\n",win);
        return -1;
    }

    /* node points to a pointer to the entry to unlink. */

    ent = *node;
    *node = ent->next; /* entry is now unlinked. */

    free(ent);

    return 0;
}

/* Retrieve the entry for the given ID, or current if win==0.
 * Return NULL if none exists.
 */

void *winReg_getPtr(int win)
{
    struct winRegEnt **node;

    win = winReg_getWindow(win);

    if (NULL==(node=winReg_find(win))) {
        /* not found */
        return NULL;
    }

    /* This should never happen. */
    if ( win != (*node)->win ) {
        fprintf(stderr,"winReg_getPtr: INCONSISTENT RESULT; exiting.\n");
        exit(-1);
    }

    return (*node)->ptr;
}


