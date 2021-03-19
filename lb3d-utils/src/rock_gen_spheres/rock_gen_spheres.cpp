#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <math.h>
#include <time.h>
#include <set>

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#define MAXSTRING 160

#define DATASET "OutArray"
//#define H5FILE "test.h5"

#define SC 0
#define BCC 1
#define FCC 2
#define RND 3

using namespace std;

int       nx,ny,nz,ns,rs,ndist;
float     ***pos;
float     w;

#ifdef ALLOWNOOVERLAP
set<int>  allowed;
int overlap;
#endif

int amod(int x, int m) {
    int r = x%m;
    return r<0 ? r+m : r;
}

void generate_sphere(int z,int y,int x,int r) {
  //fprintf(stdout,"Sphere at (%d,%d,%d)\n",x,y,z);
  //fprintf(stdout,"Checking volume (%d %d %d) - (%d %d %d)\n",x-r,y-r,z-r,x+r,y+r,z+r);
  for (int k = z-r; k <= z+r; k++) {
    //fprintf(stdout,"i%d\n",i);
    for (int j = y-r; j <= y+r; j++) {
      //fprintf(stdout,"j%d\n",j);
      for (int i = x-r; i <= x+r; i++) {
        //fprintf(stdout,"k%d\n",k);
        //fprintf(stdout,"Checking site %d %d %d\n",i,j,k);
        //fprintf(stdout,"Checking %d vs %d\n",(i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z), r*r);
        if ((i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z) <= r*r) {
          //fprintf(stdout,"Creating rock at (x,y,z) (%d,%d,%d) [%d,%d,%d]\n",amod(i,nx),amod(j,ny),amod(k,nz),i,j,k);
          pos[amod(k,nz)][amod(j,ny)][amod(i,nx)] = w;
        }
      }
    }
  }

#ifdef ALLOWNOOVERLAP
  // If overlap is not allowed, remove sites
  if (overlap == 0) {
    for (int k = z-2*r; k <= z+2*r; k++) {
      for (int j = y-2*r; j <= y+2*r; j++) {
        for (int i = x-2*r; i <= x+2*r; i++) {
          if ((i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z) <= 4*r*r) {
            //fprintf(stdout,"Disallowing site (x,y,z) (%d,%d,%d) [%d %d %d] %d\n",amod(i,nx),amod(j,ny),amod(k,nz),i,j,k,amod(k,nz)+nz*(amod(j,ny)+ny*amod(i,nx)));
            allowed.erase(amod(k,nz)+nz*(amod(j,ny)+ny*amod(i,nx)));
          }
        }
      }
    }
  }
#endif
}

int main (int argc, char *argv[]) {
#ifdef WRITETXT
  FILE   *f1;
  char   txtfile[MAXSTRING];
#endif
#ifdef WRITEXDR
  FILE   *f2;
  char   xdrfile[MAXSTRING];
  char   command[MAXSTRING];
  XDR    xdrs;
#endif
  int    x, y, z, tns, cnt;
  char   *fname;
  char   h5file[MAXSTRING];
#ifdef ALLOWNOOVERLAP
  set<int>::iterator it;
  int    p, c;
#endif
  clock_t t1, t2;

  //unsigned int szip_options_mask;
  //unsigned int szip_pixels_per_block;

#ifdef ALLOWNOOVERLAP
  overlap = 1;
#endif

  fprintf(stdout,"=== COmposite SPheres ROck GENerator === \n");

  srand ( time(NULL) );
  if (argc!=2) {
    fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
    exit(1);
  }

  fprintf(stdout,"Enter nx: \n");
  scanf("%d",&nx);
  fprintf(stdout,"Enter ny: \n");
  scanf("%d",&ny);
  fprintf(stdout,"Enter nz: \n");
  scanf("%d",&nz);
  fprintf(stdout,"Enter wettability: \n");
  scanf("%f",&w);
  fprintf(stdout,"Enter requested distribution (0 - SC, 1 - BCC, 2 - FCC, 3 - RND): \n");
  scanf("%d",&ndist);
  
  if ( ndist == SC ) {
    fprintf(stdout,"Generating spheres on a simple cubic lattice \n");
  }
  else if ( ndist == BCC ) {
    fprintf(stdout,"Generating spheres on a body-centered cubic lattice \n");
  }
  else if ( ndist == FCC ) {
    fprintf(stdout,"Generating spheres on a face-centered cubic lattice \n");
  }
  else if ( ndist == RND ) {
    fprintf(stdout,"Generating randomly distributed spheres \n");
  }
  else {
    fprintf(stderr,"Unknown distribution, aborting \n");
    exit(1);
  }

  if (ndist == SC || ndist == BCC || ndist == FCC) {
    fprintf(stdout,"Enter cell edge size (recommended: a divisor of the lattice size): \n");
    scanf("%d",&ns);
    fprintf(stdout,"Enter sphere size: \n");
    scanf("%d",&rs);
  }
  else if (ndist == RND) {
    fprintf(stdout,"Enter number of spheres: \n");
    scanf("%d",&ns);
    fprintf(stdout,"Enter sphere size: \n");
    scanf("%d",&rs);
#ifdef ALLOWNOOVERLAP
    fprintf(stdout,"Allow overlap (0/1)? \n");
    scanf("%d",&overlap);
#endif
  }

  fname = argv[1];

  fprintf(stdout,"Creating geometry for a %d x %d x %d system... \n",nx,ny,nz);
  fprintf(stdout,"Output file(s) will be: \n");
#ifdef WRITETXT
  fprintf(stdout,"\t- %s.txt\n",fname);
  sprintf(txtfile,"%s%s",fname,".txt");
#endif
#ifdef WRITEXDR
  fprintf(stdout,"\t- %s.xdr\n",fname);
  sprintf(xdrfile,"%s%s",fname,".xdr");
#endif
  fprintf(stdout,"\t- %s.h5\n",fname);
  sprintf(h5file,"%s%s",fname,".h5");

#ifdef WRITETXT
  /* Open ASCII file */
  f1 = fopen(txtfile,"w");
#endif
#ifdef WRITEXDR
  /* Create xdrfile, stop if it already exists/overwrite? */
  sprintf(command,"touch %s ",xdrfile);
  system(command);
  if(NULL==(f2=fopen(xdrfile,"w"))){
    fprintf(stderr,"Unable to open output file, aborting... ");
    exit(1);
  }
#endif
  /* Allocate tensor (i.e. 3d array) and set pointers */
  fprintf(stdout,"Setting up 3d array... \n");
  pos = (float ***) malloc(nz*sizeof(float **));
  if (!pos) {
    fprintf(stderr,"Error: (1) can't allocate memory, aborting... \n");
    exit(1);
  }
  pos[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!pos[0]) { 
    fprintf(stderr,"Error: (2) can't allocate memory, aborting... \n");
    exit(1);
  }
  pos[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!pos[0][0]) {
    fprintf(stderr,"Error: (3) can't allocate memory, aborting... \n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    pos[0][j] = pos[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    pos[k]=pos[k-1]+ny;
    pos[k][0]=pos[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      pos[k][j]=pos[k][j-1]+nx;
    }
  }

  fprintf(stdout,"Finished setting up 3d array \n");

  cnt = 0;
  fprintf(stdout,"Filling lattice with %d zeroes... \n",nz*ny*nx);
  t1 = clock();
#ifdef ALLOWNOOVERLAP
  if (overlap == 0) { 
    fprintf(stdout,"(and generating initial list of allowed sites)\n");
  }
#endif
  fprintf(stdout,"Filled %d",cnt);
  for(int k=0;k<nz;k++) {
    fprintf(stdout,"\rFilled %d",cnt);
    for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
        pos[k][j][i] = 0.0;
#ifdef ALLOWNOOVERLAP
        if (overlap == 0) {
          //fprintf(stdout,"Inserting (%d,%d,%d) %d into allowed list\n",i,j,k,k+nz*(j+ny*i));
          allowed.insert(k+nz*(j+ny*i));
        }
#endif
        cnt++;
      }
    }
  }
  fprintf(stdout,"\rFilled %d!\n",cnt);
  t2 = clock();
  fprintf(stdout,"Filling took %f seconds\n",((float)t2 - (float)t1)/1000000.0);
  t1 = clock();
/*
  fprintf(stdout,"Checking decomposition routine... \n");
  it = allowed.begin();
  p = allowed.size();
  for (int i=0;i<p;i++) {
    c = *it;
    x = 0;
    y = 0;
    z = 0;
    while (c >= nz*ny) {
      x++;
      c-= nz*ny;
    }
    while (c >= nz) {
      y++;
      c-= nz;
    }
    z = c;
    fprintf(stdout,"Decomposed coordinate %d into (%d, %d, %d)\n",*it,x,y,z);
    it++;
  }
*/

  cnt = 0;
  switch(ndist) {
    case SC:
      tns = floor((float)nz/(float)ns)*floor((float)ny/(float)ns)*floor((float)nx/(float)ns);
      fprintf(stdout,"Attempting to place %d spheres... \n",tns);
      fprintf(stdout,"Completed %d",cnt);
      for (int k=0; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=0; j < ny; j += ns) {
          for (int i=0; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      fprintf(stdout,"\rCompleted %d!\n",cnt);
      break;
    case BCC:
      tns = 2*floor((float)nz/(float)ns)*floor((float)ny/(float)ns)*floor((float)nx/(float)ns);
      fprintf(stdout,"Attempting to place %d spheres... \n",tns);
      fprintf(stdout,"Completed %d",cnt);
      for (int k=0; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=0; j < ny; j += ns) {
          for (int i=0; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      for (int k=ns/2; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=ns/2; j < ny; j += ns) {
          for (int i=ns/2; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      fprintf(stdout,"\rCompleted %d!\n",cnt);
      break;
    case FCC:
      tns = 4*floor((float)nz/(float)ns)*floor((float)ny/(float)ns)*floor((float)nx/(float)ns);
      fprintf(stdout,"Attempting to place %d spheres... \n",tns);
      fprintf(stdout,"Completed %d",cnt);
      for (int k=0; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=0; j < ny; j += ns) {
          for (int i=0; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      for (int k=ns/2; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=ns/2; j < ny; j += ns) {
          for (int i=0; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      for (int k=ns/2; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=0; j < ny; j += ns) {
          for (int i=ns/2; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      for (int k=0; k < nz; k += ns) {
#ifdef SHOWPROGRESS
        fprintf(stdout,"\rCompleted %d",cnt);
#endif
        for (int j=ns/2; j < ny; j += ns) {
          for (int i=ns/2; i < nx; i += ns) {
            generate_sphere(k,j,i,rs);
            cnt++;
          }
        }
      }
      fprintf(stdout,"\rCompleted %d!\n",cnt);
      break;
    case RND:
      tns = ns;
      fprintf(stdout,"Attempting to place %d spheres... \n",tns);
      fprintf(stdout,"Completed %d",cnt);
#ifdef ALLOWNOOVERLAP
      if (overlap == 0) {
        for (int n = 0; n < ns; n++) {
#ifdef SHOWPROGRESS
          fprintf(stdout,"\rCompleted %d",cnt);
#endif
          if (allowed.size() == 0) {
            fprintf(stderr,"No more valid placements for spheres, aborting... \n");
            exit(2);
          }
          it = allowed.begin();
          p = rand() % allowed.size();
          for (int i = 0; i < p; i++) {
            it++;
          }
          c = *it;
          x = 0;
          y = 0;
          z = 0;
          while (c >= nz*ny) {
            x++;
            c-= nz*ny;
          }
          while (c >= nz) {
            y++;
            c-= nz;
          }
          z = c;
          //fprintf(stdout,"Decomposed coordinate %d into (%d, %d, %d)\n",*it,x,y,z);
          generate_sphere(z,y,x,rs);
          //fprintf(stdout,"Created sphere at (%d, %d, %d)\n",x,y,z);
          cnt++;
        }
      }
      else {
#endif
        for (int n = 0; n < ns; n++) {
#ifdef SHOWPROGRESS
          fprintf(stdout,"\rCompleted %d",cnt);
#endif
          x = rand() % nx;
          y = rand() % ny;
          z = rand() % nz;
          generate_sphere(z,y,x,rs);
          cnt++;
        }
#ifdef ALLOWNOOVERLAP
      }
#endif
      fprintf(stdout,"\rCompleted %d\n",cnt);
      break;
    default:
      fprintf(stderr,"Unknown switch, aborting... \n");
      exit(1);
  }
  t2 = clock();
  fprintf(stdout,"Succesfully placed %d spheres \n",tns);
  fprintf(stdout,"Placing took %f seconds\n",((float)t2 - (float)t1)/1000000.0);

  fprintf(stdout,"Writing data... \n");
#ifdef WRITEXDR
  xdrstdio_create(&xdrs,f2,XDR_ENCODE); 
#endif

#if defined(WRITETXT) || defined(WRITEXDR)
  for(int k=0;k<nz;++k) {
    for(int j=0;j<ny;++j) {
      for(int i=0;i<nx;++i) {
        /* Check output of xdr file*/
#ifdef WRITETXT
        fprintf(f1," %d %d %d %f \n", k,j,i,pos[k][j][i]);
#endif
#ifdef WRITEXDR
        xdr_float(&xdrs,&pos[k][j][i]);
#endif
      }
    }
  }
#endif

#ifdef WRITEXDR
  xdr_destroy(&xdrs); 
  fclose(f2);
#endif

#ifdef WRITETXT
  fclose(f1);
#endif
  //HDF variable declarations
  hid_t           file, space, dset, dcpl;    /* Handles */
  herr_t          status;
  H5D_layout_t    layout;
  hsize_t         dims[3] = {nz, ny, nx},
  // Size limit on chunks of 4GB, so use some magic numbers to try and cap it
                  chunk[3] = {1,1,1},
                  start[3] = {0, 0, 0},
                  stride[3] = {1, 1, 1},
                  count[3] = {1, 1, 1},
                  block[3] = {nz, ny, nx};
  hsize_t         chunk_read[3];

  file = H5Fcreate (h5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (3, dims, NULL);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  //status = H5Pset_chunk (dcpl, 3, chunk);
  //szip_options_mask=H5_SZIP_NN_OPTION_MASK;
  //szip_pixels_per_block=32;
  //status = H5Pset_szip (dcpl, szip_options_mask, szip_pixels_per_block);
  dset = H5Dcreate (file, DATASET, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  //dset = H5Dcreate (file, DATASET, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sselect_hyperslab (space, H5S_SELECT_SET, start, stride, count, block);
  status = H5Dwrite (dset, H5T_NATIVE_FLOAT, H5S_ALL, space, H5P_DEFAULT, pos[0][0]);
  //status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  fprintf(stdout,"Done! \n");


  fprintf(stdout,"Checking hdf output \n");
  file = H5Fopen (h5file, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset = H5Dopen (file, DATASET, H5P_DEFAULT);

  dcpl = H5Dget_create_plist (dset);
  layout = H5Pget_layout (dcpl);
    printf ("Storage layout for %s is: ", DATASET);
  switch (layout) {
      case H5D_COMPACT:
        printf("H5D_COMPACT\n");
        break;
      case H5D_CONTIGUOUS:
        printf("H5D_CONTIGUOUS\n");
        break;
      case H5D_CHUNKED:
        printf("H5D_CHUNKED\n");
        H5Pget_chunk(dcpl,3,chunk_read);
        printf("Chunk dimensions are %d x %d x %d\n",(int) chunk_read[0], (int) chunk_read[1], (int) chunk_read[2]);
        break;
      case H5D_LAYOUT_ERROR:
        printf("H5D_LAYOUT_ERROR\n");
        break;
      case H5D_NLAYOUTS:
        printf("H5D_NLAYOUTS\n");
        break;
  }

  //status = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos[0][0]);

//   printf ("\nData as written to disk by hyberslabs:\n");
//   for (int k=0; k<nz; k++) {
//       printf (" z = %d [",k);
//       for (int j=0; j<ny; j++) {
//         printf ( " y = %d (",j);
//         for (int i=0; i < nx; i++) {
//           printf (" %f", pos[k][j][i]);
//         }
//         printf(") ");
//       }
//       printf ("]\n");
//   }

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  //status = H5Sclose (space);
  status = H5Fclose (file);
  
  return 0;
}

