/* 
 *  PROGRAM BY E. BREITMOSER, EPCC, University of Edinburgh
 *  developped on HPCx, August 2004
 *
 *  writes xdr and ascii file as outoput file
 *
 * 3D, for straight channels (0 or 90 degree angle within gid) and 
 * tilted 3D channels described through a line in 3D
 * If angle ne 45 degree (or 0 or 90), some grid points are missed, 
 * result looks staggered.
 */

/*
 * produces tesfile.xdr and testfile.txt as output
 * Input file 'testfile':
 * nx
 * ny
 * nz
 * n = # of channels
 * (x1,y1,z1) for channel 1
 * (x2,y2,z2) for channel 1
 * width/radius in gridpoints for channel 1
 * ...
 * (x1,y1,z1) for channel n
 * (x2,y2,z2) for channel n
 * width/radius in gridpoints for channel n
 */

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <math.h>
float function1(int c1, int c2, int x1, int x2, int x) ;
void swap(int *p, int *q);

#define MAXSTRING 160

int main (int argc, char *argv[])
{
        int         status_n;
	FILE	    *f0,*f1,*f2;
	int         i,j,k,chan;
	XDR	    xdrs;
	char        *fname;
        char        extxdr[] = ".xdr",ext[] = ".txt";
        char        xdrfile[MAXSTRING],txtfile[MAXSTRING];
        char        command[MAXSTRING];
	int         nx,ny,nz,channels;
        float       ***pos;
	int       *tx1,*tx2,*ty1,*ty2,*tz1,*tz2 ;/*3d channels */
	int       *width; /* 'Width' of the channel*/
	float     x,y,z;
	if (argc!=2) {
		fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
		return -1;
	}

	fname=argv[1];

	/* Input file */
        f0 = fopen(fname,"r");
        /* 3D grid */
        fscanf(f0,"%d",&nx);
        fscanf(f0,"%d",&ny);
        fscanf(f0,"%d",&nz);  
        /* Number of channels in total */
        fscanf(f0,"%d",&channels);  

	/* Allocate 3D points for all channels */
	tx1  = malloc(channels * sizeof(int));
	tx2  = malloc(channels * sizeof(int));
	ty1  = malloc(channels * sizeof(int));
	ty2  = malloc(channels * sizeof(int));
	tz1  = malloc(channels * sizeof(int));
	tz2  = malloc(channels * sizeof(int));
        width = malloc(channels * sizeof(int));

	/* Read channels*/
	for(chan=0;chan<channels;++chan){
        fscanf(f0,"%d %d %d", &tx1[chan], &ty1[chan],&tz1[chan]);
        fscanf(f0,"%d %d %d", &tx2[chan], &ty2[chan],&tz2[chan]);
        fscanf(f0,"%d", &width[chan]);
        printf("Channels: %d %d %d %d %d %d %d\n",tx1[chan],ty1[chan],tz1[chan],tx2[chan],ty2[chan],tz2[chan],width[chan]);
}
        fclose(f0);

        sprintf(txtfile,"%s%s",fname,ext);

        /* Allocate tensor (i.e. 3d array) and set pointers */
          pos=(float ***)malloc(nx*sizeof(float **));
          if(!pos){fprintf(stderr,"Error: (1) can't allocate memory \n");
          exit(1);}
          pos[0]=(float **)malloc(nx*ny*sizeof(float *));
          if(!pos[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(
1);}
          pos[0][0]=(float *)calloc(nx*ny*nz,sizeof(float));
          if(!pos[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}

          
          for(j=1;j<ny;j++){pos[0][j]=pos[0][j-1]+nz;}

            for(i=1;i<nx;i++){
            pos[i]=pos[i-1]+ny;
            pos[i][0]=pos[i-1][0]+ny*nz;
            for(j=1;j<ny;j++){
              pos[i][j]=pos[i][j-1]+nz;
            }
            } 
        sprintf(txtfile,"%s%s",fname,ext);
        sprintf(xdrfile,"%s%s",fname,extxdr);
	/* Open ASCII file */
        f1 = fopen(txtfile,"w");
	/* Create xdrfile, stop if it already exists/overwrite? */
        sprintf(command,"touch %s ",xdrfile);
        system(command);
        if(NULL==(f2=fopen(xdrfile,"w"))){
	  perror("Unable to open output file");
	  return -1;
	}        
        printf("%s \n",xdrfile);

	  /* In the beginning set all points to 1 == closed/wall */ 
   
          for(k=0;k<nz;++k){
           for(j=0;j<ny;++j){
            for(i=0;i<nx;++i){
	      pos[i][j][k] = 1.0; }}}

	  /* Set all points to zero within channels channel, 3d */
	  for(chan=0;chan<channels;++chan){
           for(k=0;k<nz;++k){
            for(j=0;j<ny;++j){
             for(i=0;i<nx;++i){
	       /* swap points if for all coords if x1 > x2 */
	       if(tx1[chan] > tx2[chan]){
               swap(&tx1[chan],&tx2[chan]);
               swap(&ty1[chan],&ty2[chan]);
               swap(&tz1[chan],&tz2[chan]);}
	       /* End of swap */

	       /* Case 1 */
	       if((ty1[chan] == ty2[chan]) && (tx1[chan] != tx2[chan])){
      	        if(i>=tx1[chan] && i<=tx2[chan]){         
		  y = ty1[chan];
		  z = function1(tz1[chan],tz2[chan],tx1[chan],tx2[chan],i);
		  /*	 printf("Function: %d %d %f %f\n",i,j,y,z); */
	        if((y+(float)width[chan])>=(float)j &&   
                  (float)j>=(y-(float)width[chan]) &&
                  (z+(float)width[chan])>=(float)k &&
                  (float)k>=(z-(float)width[chan]))
		 {pos[i][j][k]=0.0; }}
	       }
	       /* Case 2 */
	       else if((tx1[chan] == tx2[chan]) && (ty1[chan] != ty2[chan])){
		 /*Swap points y and z points if y1>y2 */
                if(ty1[chan] > ty2[chan]){
                 swap(&ty1[chan],&ty2[chan]);
                 swap(&tz1[chan],&tz2[chan]);}
	       /* End of swap */
      	        if(j>=ty1[chan] && j<=ty2[chan]){
		  x =tx1[chan];
		  z = function1(tz1[chan],tz2[chan],ty1[chan],ty2[chan],j);
		  /*	 printf("Function: %d %d %f %f\n",i,j,y,z); */
	        if((x+(float)width[chan])>=(float)i &&   
                  (float)i>=(x-(float)width[chan]) &&
                  (z+(float)width[chan])>=(float)k &&
                  (float)k>=(z-(float)width[chan]))
		 {pos[i][j][k]=0.0; }}
}
	       /* Case 3 */
	       else if((tx1[chan] == tx2[chan]) && (tz1[chan] != tz2[chan])){
		 /*Swap points y and z points if z1>z2 */
                if(tz1[chan] > tz2[chan]){
                 swap(&ty1[chan],&ty2[chan]);
                 swap(&tz1[chan],&tz2[chan]);}
	       /* End of swap */
      	        if(k>=tz1[chan] && k<=tz2[chan]){
		  x = tx1[chan];
		  y = function1(ty1[chan],ty2[chan],tz1[chan],tz2[chan],k);
		  /*	 printf("Function: %d %d %f %f\n",i,j,y,z); */
	        if((x+(float)width[chan])>=(float)i &&   
                  (float)i>=(x-(float)width[chan]) &&
                  (y+(float)width[chan])>=(float)j &&
                  (float)j>=(y-(float)width[chan]))
		 {pos[i][j][k]=0.0; }}
}
	       /* Case 4*/
	       else{
   	        /* If x point is between start and end x */
      	        if(i>=tx1[chan] && i<=tx2[chan]){
		  y = function1(ty1[chan],ty2[chan],tx1[chan],tx2[chan],i);
		  z = function1(tz1[chan],tz2[chan],ty1[chan],ty2[chan],j);
		  /*	 printf("Function: %d %d %f %f\n",i,j,y,z); */
	        if((y+(float)width[chan])>=(float)j &&   
                  (float)j>=(y-(float)width[chan]) &&
                  (z+(float)width[chan])>=(float)k &&
                  (float)k>=(z-(float)width[chan]))
		 {pos[i][j][k]=0.0; }}

 		  }}}}
	  }

          xdrstdio_create(&xdrs,f2,XDR_ENCODE); 

	  for(k=0;k<nz;++k){
	   for(j=0;j<ny;++j){
            for(i=0;i<nx;++i){
	      /* Check output of xdr file*/    
             xdr_float(&xdrs,&pos[i][j][k]); 
             fprintf(f1," %d %d %d %f \n", i,j,k,pos[i][j][k]);
	}     
	}
	}
	     xdr_destroy(&xdrs); 
       
       	     fclose(f1);
	     fclose(f2); 
	return 0;
	  
	     }

/* function for 3D line with
 * y=y1+(x-x1)/(x2-x1)*(y2-y1) and z=z1+(y-y1)/(y2-y1)*(z2-z1) etc */

float function1(int c1, int c2, int x1, int x2, int x)
{
  if((c2-c1) == 0)
   {return ((float)c1);}
  else if ((c2-c1) !=0 && (x2-x1) == 0){printf("Attention, wrong function1, check code\n");exit(1);}
  else{
    return ((float)c1+(float)(x-x1)*(float)(c2-c1)/(float)(x2-x1));}
}

/*Swap points x1,x2 or y1,y2 or z1,z2 */
void swap(int *p, int*q)
{
  int tmp;

  tmp = *p;
  *p = *q;
  *q = tmp;
}

