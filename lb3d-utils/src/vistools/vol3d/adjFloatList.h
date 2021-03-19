#ifndef INCLUDED_ADJFLOATLIST_H
#define INCLUDED_ADJFLOATLIST_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <stdio.h>
#include <stdlib.h>


struct adjFloat {
    struct adjFloat *next;
    struct adjFloat *prev;

    GLfloat *val;
    GLfloat max,min;
    GLfloat delta;
    char *name;

    void (*callback)(float, void *);
    void *callbackData;
};

struct adjFloatList {
    struct adjFloat *head;
    struct adjFloat *current;
    int wrapOn;
};


struct adjFloatList *adjFloatList_new(void);
int adjFloatList_destroy(struct adjFloatList *self);
struct adjFloat *adjFloat_new(GLfloat *val, GLfloat max, GLfloat min, GLfloat delta, char *name);
int adjFloatList_addFloat(struct adjFloatList *self, GLfloat *val, GLfloat min, GLfloat max, GLfloat delta, char *name);
int adjFloatList_dump(struct adjFloatList *self, FILE *fh);
struct adjFloat *adjFloatList_find(struct adjFloatList *self, GLfloat *val);
int adjFloatList_setCallback(struct adjFloatList *self, GLfloat *val, void (*callback)(float, void *), void *callbackData);
int adjFloatList_removeFloat(struct adjFloatList *self, GLfloat *val);
GLfloat adjFloat_clip(struct adjFloat *self);
GLfloat adjFloatList_incFloat(struct adjFloatList *self);
GLfloat adjFloatList_decFloat(struct adjFloatList *self);
int adjFloatList_next(struct adjFloatList *self);
int adjFloatList_prev(struct adjFloatList *self);
int adjFloatList_wrapOn(struct adjFloatList *self);
int adjFloatList_wrapOff(struct adjFloatList *self);
int adjFloatList_toggleWrap(struct adjFloatList *self);

char *adjFloatList_getCurrentName(struct adjFloatList *self);
GLfloat adjFloatList_getCurrentValue(struct adjFloatList *self);
#endif /* INCLUDED_ADJFLOATLIST_H */
