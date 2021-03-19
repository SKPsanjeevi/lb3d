#ifndef INCLUDED_VOLRENDER_H
#define INCLUDED_VOLRENDER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "vol3d.h"
#include "vol3df.h"
#include "lut.h"
#include "quat.h"
#include "lutWidget.h"
#include "glutViewer.h"
#include "adjFloatList.h"

#include "dirscan.h"

#include "screenshot.h"

#include "volRender.h"

enum volRenderFileMode {
    VOLRENDER_READ = 0,
    VOLRENDER_WRITE = 1
};

struct volRender {
    GLuint texid;
    struct glutViewer *gv;
    struct adjFloatList *afl;

    unsigned int nVolumes;
    int currVol;
    int drawBox;
    struct vol3d **quantized ,*vrgba;
    char **name;

    struct u8_rgba_lut *lut;
    struct u8_freqtable **ft;
    struct lutWidget *lw;

    float rangeMax,rangeMin;
    int useRange;

    /* Current number of plane slices through volume, and
     * low-res and high-res values of same.
     */
    int slices,slicesLow, slicesHigh;
    /* Times stored as unsigned 32-bit milliseconds since the epoch.
     * lastClickTime is the time at which the last event was performed
     * at low LOD. lodTime is the time to wait before rendering at high LOD.
     */
    long lastClickTime,lodTime;


    GLfloat fogStart,fogEnd,fogDensity;

    GLfloat globalAlpha;

    int fogOn;
    int wrapFrames;

    time_t newestTime,minAge; /* Parameters used for dirscanning mode */
    char *prefix;             /* Filename prefix in dirscanning mode */

};

/* volRender.c */
matrix matrix_fromGLarray(GLfloat *m);
void getNearestCubeVecs(vector *forwardp, vector *upp, GLfloat *m);
void glVertex_vec(vector v);
void snapToCube(struct glutViewer *gv);
void initCubeVectors(void);
int volRender_xdr(struct volRender *self, XDR *xdrp);
int volRender_serializeFH(struct volRender *self, FILE *f, enum volRenderFileMode mode);
int volRender_serializeFile(struct volRender *self, char *fname, enum volRenderFileMode mode);
long gettime_ms(void);
int checkGLerror(char *mess);
int vol3d_f_to_u8_ftable(struct vol3d *v, struct vol3d **quantp, struct u8_freqtable **tablep, float phimax, float phimin);
int rgba_write_png_slice(struct vol3d *v, char *filename, unsigned int z);
void vertex(GLfloat tx, GLfloat ty, GLfloat tz, GLfloat x, GLfloat y, GLfloat z);
void myDisplay(struct glutViewer *gv, void *data);
void vrgba_cheat(struct vol3d *quantized, struct vol3d *vrgba);
int volRender_setTexture(struct volRender *self);
int volRender_reload(struct volRender *self);
void myKeyCallback(struct lutWidget *self, unsigned char c, void *ptr);
int volRender_destroy(struct volRender *self);
struct volRender *volRender_new(int nvols);
void volRender_keyboardFunc(unsigned char key, int x, int y);
void volRender_advanceFrame(struct volRender *self, int delta);
void wheelUpFunc(struct glutViewer *gv, void *data);
void wheelDownFunc(struct glutViewer *gv, void *data);
int volRender_quantizeVolume(struct volRender *self, struct vol3d *vol, int i);
int volRender_loadVolume(struct volRender *self, char *filename, int i);
void volRender_scannerTimerFunc(int win);
void volRender_lodTimerFunc(int win);
void volRender_motionFunc(int x, int y);
#endif /* INCLUDED_VOLRENDER_H */
