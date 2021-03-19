

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <math.h>
#include "quat.h"


/*********** Quaternion stuff ****************/

/* Turn a quaternion into a matrix */

void quatmatrix(quat_t q,GLfloat *m )
{
        GLfloat w,x,y,z;

        w = q.q0; x = q.q1; y=q.q2; z=q.q3;

        m[ 0] = 1 - 2*y*y - 2*z*z;
        m[ 1] = 2*x*y + 2*w*z;
        m[ 2] = 2*x*z - 2*w*y;
        m[ 3] = 0;

        m[ 4] = 2*x*y - 2*w*z;
        m[ 5] = 1 - 2*x*x - 2*z*z;
        m[ 6] = 2*y*z + 2*w*x;
        m[ 7] = 0;

        m[ 8] = 2*x*z + 2*w*y;
        m[ 9] = 2*y*z - 2*w*x;
        m[10] = 1 - 2*x*x - 2*y*y;
        m[11] = 0;

        m[12] = 0;
        m[13] = 0;
        m[14] = 0;
        m[15] = 1;

        return;

}

/* Turn a quaternion into a matrix and its inverse (ie transpose) */

void quatmatrixandinv(quat_t q,GLfloat *m,GLfloat *mi )
{
        GLfloat w,x,y,z;

        w = q.q0; x = q.q1; y=q.q2; z=q.q3;

        mi[ 0] = m[ 0] = 1 - 2*y*y - 2*z*z;
        mi[ 4] = m[ 1] = 2*x*y + 2*w*z;
        mi[ 8] = m[ 2] = 2*x*z - 2*w*y;
        mi[12] = m[ 3] = 0;

        mi[ 1] = m[ 4] = 2*x*y - 2*w*z;
        mi[ 5] = m[ 5] = 1 - 2*x*x - 2*z*z;
        mi[ 9] = m[ 6] = 2*y*z + 2*w*x;
        mi[13] = m[ 7] = 0;

        mi[ 2] = m[ 8] = 2*x*z + 2*w*y;
        mi[ 6] = m[ 9] = 2*y*z - 2*w*x;
        mi[10] = m[10] = 1 - 2*x*x - 2*y*y;
        mi[14] = m[11] = 0;

        mi[ 3] = m[12] = 0;
        mi[ 7] = m[13] = 0;
        mi[11] = m[14] = 0;
        mi[15] = m[15] = 1;

        return;

}

/* Construct the quaternion corresponding to a right-handed rotation of
 * theta radians about the given vector.
 */

quat_t quatrotation(GLfloat theta, GLfloat x, GLfloat y, GLfloat z)
{
        GLfloat st,ct;
        quat_t q;

        st = sin(0.5*theta);
        ct = cos(0.5*theta);

        q.q0 = ct;
        q.q1 = x*st;  
        q.q2 = y*st;  
        q.q3 = z*st;  
        return q;
}
/*
 * Return the quaternion product of two quaternions.
 *
 * ( a0, a ) ( b0, b) =  ( a0*b0 - a.b , a0*a + b0*b + a x b )
 *
 */

quat_t quatmultiply(quat_t a, quat_t b)
{
        quat_t q;

        q.q0 = a.q0 * b.q0 
                - a.q1*b.q1
                - a.q2*b.q2
                - a.q3*b.q3;
        
        q.q1 =    a.q0 * b.q1
                + b.q0 * a.q1
                + a.q2 * b.q3 - a.q3 * b.q2;

        q.q2 =    a.q0 * b.q2
                + b.q0 * a.q2
                + a.q3 * b.q1 - a.q1 * b.q3;

        q.q3 =    a.q0 * b.q3
                + b.q0 * a.q3
                + a.q1 * b.q2 - a.q2 * b.q1;
        return q;
}

/* Take a 3x3 rotation matrix. Return the corresponding
 * quaternion.
 */
quat_t quatfromrotmatrix(GLfloat *mat)
{
    GLfloat x,y,z,w,x2,y2,z2,w2;
    GLfloat max;
    quat_t q;
    int largest=0; /* Index of largest quat component */

    /* Squared values can be taken from the diagonal. */

    x2 = ( 1.0 + mat[0] - mat[4] - mat[8] ) /4.0;
    y2 = ( 1.0 - mat[0] + mat[4] - mat[8] ) /4.0;
    z2 = ( 1.0 - mat[0] - mat[4] + mat[8] ) /4.0;
    w2 = ( 1.0 + mat[0] + mat[4] + mat[8] ) /4.0;

    max=w2;
    if (x2>max) { max=x2; largest=1; }
    if (y2>max) { max=y2; largest=2; }
    if (z2>max) { max=z2; largest=3; }

    switch(largest) {
        case 0: /* w */
            w = sqrt(w2);
            x = (mat[5] - mat[7])/(4.0*w);
            y = (mat[6] - mat[2])/(4.0*w);
            z = (mat[1] - mat[3])/(4.0*w);
            break;

        case 1: /* x */
            x = sqrt(x2);
            w = (mat[5] - mat[7])/(4.0*x);
            y = (mat[1] + mat[3])/(4.0*x);
            z = (mat[2] + mat[6])/(4.0*x);
            break;

        case 2: /* y */
            y = sqrt(y2);
            w = (mat[6] - mat[2])/(4.0*y);
            x = (mat[1] + mat[3])/(4.0*y);
            z = (mat[5] + mat[7])/(4.0*y);
            break;

        case 3: /* z */
        default:
            z = sqrt(z2);
            w = (mat[1] - mat[3])/(4.0*z);
            x = (mat[2] + mat[6])/(4.0*z);
            y = (mat[5] + mat[7])/(4.0*z);
            break;

    }

    q.q0=w; q.q1=x; q.q2=y; q.q3=z;

    return q;

}
