#include "adjFloatList.h"

/*
 * Points to note:
 * (1) The array of adjFloat elements may be reallocated or moved
 *     at *any* time.
 * (2) currently, the "name" element of each adjFloat is not freed
 *     at destroy time. I could add a "freeName" flag if needs be.
 */

/*
 * Create a new adjFloatList.
 * Return a pointer to it; return NULL on error.
 */

struct adjFloatList *adjFloatList_new(void)
{
    struct adjFloatList *self=NULL;
    if (NULL==(self=malloc(sizeof(struct adjFloatList))))
    { perror("adjFloatList_new: malloc"); return NULL; }

    self->head = self->current = NULL;

    return self;

}

/* Destroy an adjFloatList, deallocating all of its entries.
 * Return 0 on success, nonzero on error.
 */

int adjFloatList_destroy(struct adjFloatList *self)
{
    struct adjFloat *el,*prev;

    if (NULL==self)
    { fprintf(stderr,"adjFloatList_destroy(NULL)!\n"); return -1; }




    if (NULL!=self->head) {

        /* Walk the list and free every element. */

        el = self->head;

        do {
            prev=el;
            el = el->next;
            free(prev);
        } while (self->head != el);
    }

    free(self);
    return 0;
}

struct adjFloat *adjFloat_new(GLfloat *val, GLfloat max, GLfloat min, GLfloat delta, char *name)
{
    struct adjFloat *self=NULL;

    if (NULL==(self=malloc(sizeof(struct adjFloat))))
    { perror("adjFloat_new: malloc"); return NULL; }

    self->val = val;
    self->max = max;
    self->min = min;
    self->delta = delta;
    self->name = name;
    self->next = self->prev = NULL;
    self->callback = NULL;
    self->callbackData = NULL;
    return self;
}

/* Take a pointer to a float value and an adjFloatList. Make the float
 * value adjustable.
 * Return 0 on success, nonzero on error.
 */

int adjFloatList_addFloat(struct adjFloatList *self,
        GLfloat *val, GLfloat min, GLfloat max, GLfloat delta, char *name)
{
    struct adjFloat *head,*tail,*newel=NULL;

    if (NULL==self) {
        fprintf(stderr,"adjFloatList_addFloat(self==NULL)!\n");
        fflush(stderr);
        return -1;
    }

    if (NULL==(newel=adjFloat_new(val,max,min,delta,name))) {
        fprintf(stderr,"adjFloatList_addFloat: adjFloat_new failed.\n");
        return -2;
    }
    
    if (NULL==(head=self->head)) {
        /* First element in the list */

        self->head = self->current = newel->next = newel->prev = newel;

        return 0;
    }

    /* self->head is non-NULL, so nonzero list size.
     * Link new element onto end of list.
     */

    tail = head->prev;

    tail->next = newel; newel->prev = tail;
    head->prev = newel; newel->next = head;

    return 0;

}


/* Dump an ASCII description of an adjFloatList to a filehandle */

int adjFloatList_dump(struct adjFloatList *self, FILE *fh)
{
    struct adjFloat *el=NULL;
    fprintf(fh,"\nadjFloatList: wrapOn=%d\n",self->wrapOn);

    if (NULL==(el=self->head)) {
        fprintf(fh,"No elements.\n"); return 0;
    }

    do {
        fprintf(fh,"    (%f) max %f min %f delta %f name \"%s\" %s\n",
                *(el->val),
                el->max,
                el->min,
                el->delta,
                el->name,
                (self->current == el) ? "<-- Current" : ""
                );
        el = el->next;

    } while (self->head != el);

    fprintf(fh,"Dump done.\n\n");
    fflush(fh);

    return 0;
}

/* Take a pointer to a float. Return the corresponding adjFloat structure,
 * or NULL if not found.
 */

struct adjFloat *adjFloatList_find(struct adjFloatList *self, GLfloat *val)
{
    struct adjFloat *el=NULL;

    if (NULL==(el=self->head)) { return NULL; }

    do {
        if (val == el->val) {
            return el;
        }
        el = el->next;
    } while (self->head != el);

    return NULL;

}

/* Set a callback on the given float. This callback is called, with the
 * given data, when the float is changed.
 * Return 0 if callback was successfully changed.
 * Return nonzero otherwise, including in the case where the float
 * was not registered as adjustable.
 */

int adjFloatList_setCallback(struct adjFloatList *self, GLfloat *val,
        void (*callback)(float, void *), void *callbackData)
{
    struct adjFloat *el=NULL;

    if (NULL==(el=adjFloatList_find(self,val))) {
        fprintf(stderr,"adjFloatList_setCallback: no element found.\n");
        return -1;
    }

    el->callback = callback;
    el->callbackData = callbackData;

    return 0;

}

/* Remove a float from the adjustable list. Return 0 if it was
 * successfully removed, else nonzero.
 */

int adjFloatList_removeFloat(struct adjFloatList *self, GLfloat *val)
{
    struct adjFloat *el=NULL;

    if (NULL==(el=adjFloatList_find(self,val))) {
        fprintf(stderr,"adjFloatList_removeFloat: no element found.\n");
        return -1;
    }

    if (el->next == el) {
        /* last element of list */
        free(el);
        self->head = self->current = NULL;
        return 0;
    }

    el->next->prev = el->prev;
    el->prev->next = el->next;
    free(el);

    return 0;

}

GLfloat adjFloat_clip(struct adjFloat *self)
{

    if (*self->val >self->max) {
        *self->val  = self->max;
    }
    if (*self->val <self->min) {
        *self->val  = self->min; 
    }
    *self->val = *self->val ;

    return *self->val;

}

/* Increment the current adjustable float by delta */

GLfloat adjFloatList_incFloat(struct adjFloatList *self)
{
    if (NULL==self->current) {
        fprintf(stderr,"adjFloatList_incFloat: no elements!\n");
        return -1;
    }
    *(self->current->val) += self->current->delta;
    return adjFloat_clip(self->current);
}


/* Decrement the current adjustable float by delta */
GLfloat adjFloatList_decFloat(struct adjFloatList *self)
{
    if (NULL==self->current) {
        fprintf(stderr,"adjFloatList_decFloat: no elements!\n");
        return -1;
    }
    *(self->current->val) -= self->current->delta;
    return adjFloat_clip(self->current);
}

/* Select next adjustable float */
int adjFloatList_next(struct adjFloatList *self)
{
    if (NULL==self->current) {
        fprintf(stderr,"adjFloatList_next: no elements!\n");
        return -1;
    }
    if (self->head == self->current->next) {
        if (!self->wrapOn) { return 0; }
    }


    self->current = self->current->next;

    return 0;
}


/* Select previous adjustable float */
int adjFloatList_prev(struct adjFloatList *self)
{
    if (NULL==self->current) {
        fprintf(stderr,"adjFloatList_prev: no elements!\n");
        return -1;
    }
    if (self->head==self->current) {
        if (!self->wrapOn) {
            return 0;
        }
    }

    self->current = self->current->prev;
    return 0;
}

/* Adjust wrapping mode of adjfloat. If wrapping is on, then "next" after
 * the last element selects the first. If off (default), then it remains
 * unchanged.
 */
int adjFloatList_wrapOn(struct adjFloatList *self)
{ self->wrapOn = 1; return 0; }

int adjFloatList_wrapOff(struct adjFloatList *self)
{ self->wrapOn = 0; return 0; }

int adjFloatList_toggleWrap(struct adjFloatList *self)
{ self->wrapOn = self->wrapOn ? 0 : 1; return 0;}

char *adjFloatList_getCurrentName(struct adjFloatList *self)
{
    return self->current->name;
}
GLfloat adjFloatList_getCurrentValue(struct adjFloatList *self)
{
    return *(self->current->val);
}
