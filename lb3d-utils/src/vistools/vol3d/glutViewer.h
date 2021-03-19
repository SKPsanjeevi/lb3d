#ifndef INCLUDED_GLUTVIEWER_H
#define INCLUDED_GLUTVIEWER_H
/* glutViewer.c */

struct glutViewer {
    
    int win; /* GLUT window ID */

    GLuint width,height; /* Window geometry */

    int clickx,clicky; /* coordinates of last mouse click */
    int cursx,cursy; /* mouse cursor coordinates */
    int orthographic; /* Set for orthographic projection */
    GLfloat orthoSize; /* Orthographic projection size */
    int buttonState[3]; /* left, middle, right; 1==down */
    int stereoFlag; /* nonzero = GLX_STEREO mode */
    GLfloat eyesize; /* distance between eyes in GLX_STEREO mode */

    GLfloat bgColor[4];

    GLfloat fovy; /* Field of view (degrees) */
    GLfloat nearClip; /* Distance to near clipping plane */
    GLfloat geomx,geomy,geomz; /* Geometry position */

    quat_t orquat; /* Orientation quaternion */
    GLfloat sphererad; /* radius of trackball sphere */

    int alphaBufferFlag; /* Set to true if alpha buffer enabled */

    void *data; /* Arbitrary pointer */

    void (*displayCallback)(struct glutViewer *, void *);
    void *displayCallbackData;
    void (*wheelUpFunc)(struct glutViewer *, void *);
    void *wheelUpData;
    void (*wheelDownFunc)(struct glutViewer *, void *);
    void *wheelDownData;


};

void glutViewer_mouseRotate(struct glutViewer *self, int x0, int y0, int x1, int y1);
void glutViewer_mouseZoom(struct glutViewer *self, int dz);
void glutViewer_mouseTranslate(struct glutViewer *self, int dx, int dy);
void glutViewer_reshapeFunc(int width, int height);
void glutViewer_displayFunc(void);
void glutViewer_mouseFunc(int button, int state, int x, int y);
void glutViewer_motionFunc(int x, int y);
int glutViewer_init(struct glutViewer *self);
struct glutViewer *glutViewer_new(void);
struct glutViewer *glutViewer_new_stereo(void);
int glutViewer_setData(struct glutViewer *self, void *data);
void *glutViewer_getData(struct glutViewer *self);
quat_t glutViewer_getQuat(struct glutViewer *self);
void glutViewer_wheelUpFunc(struct glutViewer *self, void (*func)(struct glutViewer *, void *));
void glutViewer_wheelUpData(struct glutViewer *self, void *data);
void glutViewer_wheelDownFunc(struct glutViewer *self, void (*func)(struct glutViewer *, void *));
void glutViewer_wheelDownData(struct glutViewer *self, void *data);
void (*glutViewer_setDisplayCallback(struct glutViewer *self, void (*newCallback)(struct glutViewer *, void *)))(struct glutViewer *, void *);
void (*glutViewer_getDisplayCallback(struct glutViewer *self))(struct glutViewer *, void *);
void *glutViewer_setDisplayCallbackData(struct glutViewer *self, void *data);
void *glutViewer_getDisplayCallbackData(struct glutViewer *self, void *data);
void glutViewer_postRedisplay(struct glutViewer *self);
void glutViewer_toggleOrtho(struct glutViewer *self);

#endif /* INCLUDED_GLUTVIEWER_H */
