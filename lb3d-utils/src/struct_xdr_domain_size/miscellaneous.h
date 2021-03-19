#include <stdio.h>
void failed(char message[]);
void readeol(FILE *fp);

int myipow(int a, int n);
float myfpow(float a, int n);
double mydpow(double a, int n);

int *i1t(int n1);
int **i2t(int n1, int n2);
int ***i3t(int n1, int n2, int n3);

float *f1t(int n1);
float **f2t(int n1, int n2);
float ***f3t(int n1, int n2, int n3);
float ****f4t(int n1, int n2, int n3, int n4);
float *****f5t(int n1, int n2, int n3, int n4, int n5);

double *d1t(int n1);
double **d2t(int n1, int n2);
double ***d3t(int n1, int n2, int n3);

char *stradd(int n, ...);

void pdarray(int n, int m, double **a);
void pdvector(int n, double *a);
void pfarray(int n, int m, float **a);
void pfvector(int n, float *a);
void piarray(int n, int m, int **a);
void pivector(int n, int *a);

void zdarray(int n, int m, double **a);
void zdvector(int n, double *a);
void zfarray(int n, int m, float **a);
void zfvector(int n, float *a);
void ziarray(int n, int m, int **a);
void zivector(int n, int *a);
