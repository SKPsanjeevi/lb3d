/* 

ascci2xdr tool v 0.0 
date: 20.5.09
Ariel Narváez

The txt2xdr reads a ascci file with the columns
nx ny nz
i j k pore
...
...
...
...

using this information write a xdr file


compile the program:

gcc scci2xdr.c -o ascci2xdr

running:
./block input-data


output: input-data.xtr

to visualize with paraview, create first a vtk-file, using the 
ACSII-file head, included in this directory. The size and the nodes 
number must be in agreement.
linux comand:
cat head input-data.xtr > data.vtk
*/

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>


#define MAXSTRING 160

int main (int argc, char *argv[]){
  int       status_n;
  FILE	    *f0,*f2;
  int       i,j,k,i0,j0,k0,pr;
  XDR	    xdrs;
  char      *fname;
  char      extxdr[] = ".xdr"; 
  char      xdrfile[MAXSTRING];
  char      command[MAXSTRING];
  int       nx,ny,nz;
  float     ***pos;
  float	    color_pr;


  if (argc!=2){
    fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
    return -1;
  }
  
  fname=argv[1];
  
  /* Input file */
  f0 = fopen(fname,"r");

  /* reading 3D grid */
  fscanf(f0,"%d",&nx);
  fscanf(f0,"%d",&ny);
  fscanf(f0,"%d",&nz);  

  printf("nx = %d, ny = %d, nz = %d\n",nx,ny,nz);

  /* Allocate tensor (i.e. 3d array) and set pointers */
  pos=(float ***)malloc(nx*sizeof(float **));
  if(!pos){
    fprintf(stderr,"error: (1) can't allocate memory \n");
    exit(1);
  }
  pos[0]=(float **)malloc(nx*ny*sizeof(float *));
  if(!pos[0]){
    fprintf(stderr,"error: (2) can't allocate memory \n"); 
    exit(1);
  }
  pos[0][0]=(float *)calloc(nx*ny*nz,sizeof(float));
  if(!pos[0][0]){
    fprintf(stderr,"error: (3) can't allocate memory \n");
    exit(1);
  }
  
  for(j=1;j<ny;j++){
    pos[0][j]=pos[0][j-1]+nz;
  }

  for(i=1;i<nx;i++){
    pos[i]=pos[i-1]+ny;
    pos[i][0]=pos[i-1][0]+ny*nz;
    for(j=1;j<ny;j++){
      pos[i][j]=pos[i][j-1]+nz;
    }
  } 
   
  sprintf(xdrfile,"%s%s",fname,extxdr); 

  /* Create xdrfile, stop if it already exists/overwrite? */
  sprintf(command,"touch %s ",xdrfile);
  system(command);
  if(NULL==(f2=fopen(xdrfile,"w"))){
    perror("unable to open output file");
    return -1;
  }        
  printf("%s \n",xdrfile);

  /* Reading pores */
  for(pr=0;pr<nx*ny*nz;++pr){
  fscanf(f0,"%d %d %d %f", &i, &j, &k, &color_pr); 
  
  if (pr == 0) 
    {         
      i0 = i;
      j0 = j;
      k0 = k;
    }

  i=i-i0;
  j=j-j0;
  k=k-k0;
    
  pos[i][j][k] = color_pr;
  }

  fclose(f0);

  xdrstdio_create(&xdrs,f2,XDR_ENCODE); 
  
  for(k=0;k<nz;++k){
    for(j=0;j<ny;++j){
      for(i=0;i<nx;++i){
	/* Check output of xdr file*/
	xdr_float(&xdrs,&pos[i][j][k]); 
      }     
    }
  }
  xdr_destroy(&xdrs); 
       
  fclose(f2); 
  return 0;
}






