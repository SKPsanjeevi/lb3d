#ifndef INCLUDED_VECTOR_H
#define INCLUDED_VECTOR_H

struct vector {
        float x,y,z;
};

typedef struct vector vector;
vector vector_new(float x, float y, float z);
vector vector_add(vector a, vector b);
vector vector_sub(vector a, vector b);
vector vector_cross(vector a, vector b);
float vector_modulus(vector v);
double triangle_area(vector a, vector b, vector c);
vector vector_scale(vector a, float s);
vector vector_interpolate(vector ra, float phia, vector rb, float phib, float isoval);
vector vector_normalized(vector v);
vector vector_negnormalized(vector v);

typedef struct { float el[9]; } matrix;
matrix matrix_zero(void);
matrix matrix_identity(void);
matrix matrix_outer_vectors(vector a, vector b);
vector vector_matrix_dot(matrix m, vector a);
float vector_dot(vector a, vector b);
matrix matrix_product(matrix a, matrix b);
float matrix_determinant(matrix m) ;
float matrix_trace(matrix m) ;
float matrix_frob(matrix m) ;
matrix matrix_add(matrix a, matrix b);
matrix matrix_sub(matrix a, matrix b);
matrix matrix_scale(matrix a, float s);

#endif /* INCLUDED_VECTOR_H */
