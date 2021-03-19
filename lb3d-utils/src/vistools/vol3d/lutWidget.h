#ifndef INCLUDED_LUTWIDGET_H
#define INCLUDED_LUTWIDGET_H

#include "vol3d.h"
#include "vol3df.h"

#include "vol2df.h"

#include "winReg.h"
#include "vector.h"

#include "lut.h"
#include "vol3df.h"

#define WHEEL_UP 3
#define WHEEL_DOWN 4

#define LUTWIDGET_LEFT_BUTTON 0
#define LUTWIDGET_MIDDLE_BUTTON 1
#define LUTWIDGET_RIGHT_BUTTON 2

struct lutWidget {
    int win; /* GLUT window ID */
    GLuint cursorList,freqGraphList;
    int cursorOn;
    GLuint tfList[4]; /* display list for R,G,B,A TF graphs */
    struct u8_freqtable *freqtable;
    struct u8_rgba_lut *lut;
    GLuint width,height;
    GLint tfFlags[4]; /* flags for R, G, B, A modification */
    GLfloat clickx,clicky; /* coordinates of last mouse click */
    GLfloat cursx,cursy; /* mouse cursor coordinates */
    int buttonState[3];
    void (*keyCallback)(struct lutWidget *,unsigned char, void *);
    void *keyCallbackData;
};

struct u8_rgba_lut *u8_rgba_lut_default(void);
struct vol3d *vol3d_u8_lut_to_rgba(struct vol3d *u8, struct u8_rgba_lut *lut);
void lutWidget_postRedisplay(struct lutWidget *self);
int lutWidget_setFT(struct lutWidget *self, struct u8_freqtable *freqtable);
int lutWidget_makeFreqGraph(struct lutWidget *self);
int lutWidget_makeTFGraph(struct lutWidget *self);
int lutWidget_drawTFGraph(struct lutWidget *self);
int lutWidget_drawFreqGraph(struct lutWidget *self);
void lutWidget_win2ortho(struct lutWidget *self, int xw, int yw, GLfloat *xgp, GLfloat *ygp);
void lutWidget_ortho2win(struct lutWidget *self, GLfloat xg, GLfloat yg, int *xwp, int *ywp);
void lutWidget_setTFPoint(struct lutWidget *self, unsigned int voxelval, unsigned int tfval);
void lutWidget_leftClick(struct lutWidget *self, GLfloat tx, GLfloat ty);
void lutWidget_leftDrag(struct lutWidget *self, GLfloat tx, GLfloat ty);
void lutWidget_motionFunc(int x, int y);
void lutWidget_passiveMotionFunc(int x, int y);
void lutWidget_reshapeFunc(int w, int h);
void lutWidget_wheelUp(struct lutWidget *self);
void lutWidget_wheelDown(struct lutWidget *self);
void lutWidget_mouseFunc(int button, int state, int x, int y);
void lutWidget_drawCursor(struct lutWidget *self);
void lutWidget_newCursor(struct lutWidget *self);
void lutWidget_changeCursor(struct lutWidget *self);
void lutWidget_toggleTFComponent(struct lutWidget *self, int component);
void lutWidget_keyboardFunc(unsigned char key, int x, int y);
void lutWidget_displayFunc(void);
void lutWidget_entryFunc(int state);
int lutWidget_initWindow(struct lutWidget *self);
void lutWidget_init(struct lutWidget *widget);
struct lutWidget *lutWidget_new(struct u8_freqtable *freqtable, struct u8_rgba_lut *lut);
void lutWidget_setKeyCallback(struct lutWidget *self, void (*callback)(struct lutWidget *, unsigned char, void *), void *data);

#endif /* INCLUDED_LUTWIDGET_H */
