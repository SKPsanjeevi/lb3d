#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stddef.h>

#define NR_END 1
#define FREE_ARG char*

#include "miscellaneous.h"
#include "nrutil.h"



int *i1t(int n1)
{
  int *p,*a;
  int i;
  if((p = (int *)malloc((size_t)n1*sizeof(int))) == NULL)
    failed("i1t: failed");
  for(i = 0, a = p; i < n1; i++)*a++ = 0;
  return p;
}

int **i2t(int n1, int n2)
{
  int **p, *a;
  int i;
  if((p = (int **)malloc((size_t)n1*sizeof(int *))) == NULL)
    failed("i2t: failed n1");
  if((p[0] = (int *)malloc((size_t)n1*n2*sizeof(int))) == NULL)
    failed("i2t: failed n2");
  
  /* having allocated the memory, now do pointer arithmetic to 
     make it a 2-D array. Not the most obvious way to do this! -BG */  
  for(i = 0; i < n1-1; i++)p[i+1] = p[i] + n2;
  for(i = 0, a = p[0]; i < n1*n2; i++)*a++ = 0;
  return p;
}

int ***i3t(int n1, int n2, int n3)
{
  int ***p, *a;
  int i,j;
  if((p = (int ***)malloc((size_t)n1*sizeof(int **))) == NULL)
    failed("i3t: failed n1");
  if((p[0] = (int **)malloc((size_t)n1*n2*sizeof(int *))) == NULL)
    failed("i3t: failed n2");
  if((p[0][0] = (int *)malloc((size_t)n1*n2*n3*sizeof(int))) == NULL)
    failed("i3t: failed n3");
  for(i = 0; i < n2-1; i++)p[0][i+1] = p[0][i] + n3;
  for(i = 0; i < n1-1; i++)
    {
      p[i+1] = p[i] + n2;
      p[i+1][0] = p[i][0] + n2*n3;
      for(j = 0; j < n2-1; j++)p[i+1][j+1] = p[i+1][j] + n3;
    }
  for(i = 0, a = p[0][0]; i < n1*n2*n3; i++)*a++ = 0;
  return p;
}
	
float *f1t(int n1)
{
  float *p,*a;
  int i;
  if((p = (float *)malloc((size_t)n1*sizeof(float))) == NULL)
    failed("f1t: failed");
  for(i = 0, a = p; i < n1; i++)*a++ = 0;
  return p;
}

float **f2t(int n1, int n2)
{
  float **p, *a;
  int i;
  if((p = (float **)malloc((size_t)n1*sizeof(float *))) == NULL)
    failed("f2t: failed n1");
  if((p[0] = (float *)malloc((size_t)n1*n2*sizeof(float))) == NULL)
    failed("f2t: failed n2");
  for(i = 0; i < n1-1; i++)p[i+1] = p[i] + n2;
  for(i = 0, a = p[0]; i < n1*n2; i++)*a++ = 0;
  return p;
}

float ***f3t(int n1, int n2, int n3)
{
  float ***p, *a;
  int i,j;
  if((p = (float ***)malloc((size_t)n1*sizeof(float **))) == NULL)
    failed("f3t: failed n1");
  if((p[0] = (float **)malloc((size_t)n1*n2*sizeof(float *))) == NULL)
    failed("f3t: failed n2");
  if((p[0][0] = (float *)malloc((size_t)n1*n2*n3*sizeof(float))) == NULL)
    failed("f3t: failed n3");
  for(i = 0; i < n2-1; i++)p[0][i+1] = p[0][i] + n3;
  for(i = 0; i < n1-1; i++)
    {
      p[i+1] = p[i] + n2;
      p[i+1][0] = p[i][0] + n2*n3;
      for(j = 0; j < n2-1; j++)p[i+1][j+1] = p[i+1][j] + n3;
    }
  for(i = 0, a = p[0][0]; i < n1*n2*n3; i++)*a++ = 0;
  return p;
}

double *d1t(int n1)
{
  double *p,*a;
  int i;
  if((p = (double *)malloc((size_t)n1*sizeof(double))) == NULL)
    failed("d1t: failed");
  for(i = 0, a = p; i < n1; i++)*a++ = 0;
  return p;
}

double **d2t(int n1, int n2)
{
  double **p, *a;
  int i;
  if((p = (double **)malloc((size_t)n1*sizeof(double *))) == NULL)
    failed("d2t: failed n1");
  if((p[0] = (double *)malloc((size_t)n1*n2*sizeof(double))) == NULL)
    failed("d2t: failed n2");
  for(i = 0; i < n1-1; i++)p[i+1] = p[i] + n2;
  for(i = 0, a = p[0]; i < n1*n2; i++)*a++ = 0;
  return p;
}

double ***d3t(int n1, int n2, int n3)
{
  double ***p, *a;
  int i,j;
  if((p = (double ***)malloc((size_t)n1*sizeof(double **))) == NULL)
    failed("d3t: failed n1");
  if((p[0] = (double **)malloc((size_t)n1*n2*sizeof(double *))) == NULL)
    failed("d3t: failed n2");
  if((p[0][0] = (double *)malloc((size_t)n1*n2*n3*sizeof(double))) == NULL)
    failed("d3t: failed n3");
  for(i = 0; i < n2-1; i++)p[0][i+1] = p[0][i] + n3;
  for(i = 0; i < n1-1; i++)
    {
      p[i+1] = p[i] + n2;
      p[i+1][0] = p[i][0] + n2*n3;
      for(j = 0; j < n2-1; j++)p[i+1][j+1] = p[i+1][j] + n3;
    }
  for(i = 0, a = p[0][0]; i < n1*n2*n3; i++)*a++ = 0;
  return p;
}

	
void failed(char message[])
{
  printf("\n *** %s ***\n \n", message);
  exit(1);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}





float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        float ***t;

        /* allocate pointers to pointers to rows */
        t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
        if (!t) nrerror("allocation failure 1 in f3tensor()");
        t += NR_END;
        t -= nrl;

        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
        if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
        t[nrl] += NR_END;
        t[nrl] -= ncl;

        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
        if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
        t[nrl][ncl] += NR_END;
        t[nrl][ncl] -= ndl;

        for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++) {
                t[i]=t[i-1]+ncol;
                t[i][ncl]=t[i-1][ncl]+ncol*ndep;
                for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

        /* return pointer to array of pointers to rows */
        return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */

