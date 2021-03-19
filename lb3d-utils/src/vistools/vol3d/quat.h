#ifndef INCLUDED_QUAT_H
#define INCLUDED_QUAT_H

struct quat {
        GLfloat q0,q1,q2,q3;
};
typedef struct quat quat_t;

void quatmatrix(quat_t q, GLfloat *m);
void quatmatrixandinv(quat_t q, GLfloat *m, GLfloat *mi);
quat_t quatrotation(GLfloat theta, GLfloat x, GLfloat y, GLfloat z);
quat_t quatmultiply(quat_t a, quat_t b);
quat_t quatfromrotmatrix(GLfloat *mat);

#endif /* INCLUDED_QUAT_H */
