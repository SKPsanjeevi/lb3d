
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"

/* 3-vector manipulation */

vector vector_new(float x, float y, float z)
{
    vector v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

vector vector_add (vector a, vector b)
{
        vector c;
        c.x = a.x + b.x; c.y = a.y + b.y; c.z = a.z + b.z;
        return c;
}

vector vector_sub (vector a, vector b)
{
        vector c;
        c.x = a.x - b.x; c.y = a.y - b.y; c.z = a.z - b.z;
        return c;
}

vector vector_cross(vector a, vector b)
{
        vector c;
        c.x = a.y * b.z - a.z * b.y;
        c.y = a.z * b.x - a.x * b.z;
        c.z = a.x * b.y - a.y * b.x;
        return c;
}

float vector_modulus(vector v)
{
        return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}


/* Return the area of the triangle defined by three vertices */

double triangle_area(vector a, vector b, vector c)
{
        vector v1,v2;

        v1 = vector_sub(c,a);
        v2 = vector_sub(b,a);

        return 0.5*vector_modulus(vector_cross(v1,v2));
}

vector vector_scale(vector a, float s)
{
        vector b;
        b.x = s * a.x;
        b.y = s * a.y;
        b.z = s * a.z;
        return b;
}

/* Take two vertices, the associated scalars, and an isovalue.
 * Calculate, by linear interpolation, the point on the line
 * at which the scalar is equal to the isovalue. Return this point.
 */
vector vector_interpolate(vector ra, float phia, 
                vector rb, float phib,
                float isoval)
{
        float lambda;
        vector delta;

        lambda = (isoval - phia)/(phib - phia);

        delta = vector_sub(rb,ra);
        delta = vector_scale(delta,lambda);
        return vector_add(ra,delta);
}

vector vector_normalized(vector v)
{
        float mag;
        mag = vector_modulus(v);
        return vector_scale(v,1.0/mag);

}

vector vector_negnormalized(vector v)
{
        float mag;
        mag = vector_modulus(v);
        return vector_scale(v,-1.0/mag);

}

matrix matrix_zero(void)
{
    matrix m;
    int i;
    for (i=0;i<9;i++) { m.el[i]=0.0; }
    return m;
}

matrix matrix_identity(void)
{
    matrix m;
    m.el[0] = m.el[4] = m.el[8] = 1.0;
    m.el[1] = m.el[2] = m.el[3] =
    m.el[5] = m.el[6] = m.el[7] = 0.0;
    return m;
}

matrix matrix_outer_vectors(vector a, vector b)
{
    matrix m;
    m.el[0] = a.x * b.x;
    m.el[1] = a.x * b.y;
    m.el[2] = a.x * b.z;
    m.el[3] = a.y * b.x;
    m.el[4] = a.y * b.y;
    m.el[5] = a.y * b.z;
    m.el[6] = a.z * b.x;
    m.el[7] = a.z * b.y;
    m.el[8] = a.z * b.z;

    return m;
}

float vector_dot(vector a, vector b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

vector vector_matrix_dot(matrix m, vector a)
{
    vector p;

    p.x = m.el[0]*a.x + m.el[1]*a.y + m.el[2]*a.z;
    p.y = m.el[3]*a.x + m.el[4]*a.y + m.el[5]*a.z;
    p.z = m.el[6]*a.x + m.el[7]*a.y + m.el[8]*a.z;

    return p;
}


matrix matrix_product(matrix a, matrix b)
{
    int i,j,k;
    matrix m;

    m = matrix_zero();

    for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
        for (k=0;k<3;k++) {
            m.el[3*i+j] += a.el[3*i+k] * b.el[3*k+j];
        }
    }}

    return m;

}

float matrix_determinant(matrix m) {

        return   m.el[0]*m.el[4]*m.el[8]  
               - m.el[0]*m.el[5]*m.el[7] 
               - m.el[1]*m.el[3]*m.el[8] 
               + m.el[1]*m.el[5]*m.el[6] 
               + m.el[2]*m.el[3]*m.el[7] 
               - m.el[2]*m.el[4]*m.el[6];
}

float matrix_trace(matrix m) {
    return m.el[0]+m.el[4]+m.el[8];
}

/* Frobenius norm */
float matrix_frob(matrix m) {
    float sum=0.0;
    int i;

    for (i=0;i<9;i++) {
        sum += m.el[i]*m.el[i];
    }
    return sqrt(sum);
}

matrix matrix_add(matrix a, matrix b)
{
    int i; matrix c;
    for (i=0;i<9;i++) { c.el[i] = a.el[i] + b.el[i]; }
    return c;
}

matrix matrix_sub(matrix a, matrix b)
{
    int i; matrix c;
    for (i=0;i<9;i++) { c.el[i] = a.el[i] - b.el[i]; }
    return c;
}

matrix matrix_scale(matrix a, float s)
{
    int i; matrix m;
    for (i=0;i<9;i++) { m.el[i] = a.el[i]*s; }
    return m;
}

