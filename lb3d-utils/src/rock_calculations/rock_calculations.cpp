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

#ifdef H5_PARALLEL
#include "hdf5.h"
#define H5_NO_NAMESPACE
#else
#include "H5Cpp.h"
#endif

#define DATASET "OutArray"

#define MAXSTRING 160
#define RANK 3

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

#define DARCYM2 9.869233E-13

using namespace std;

int main (int argc, char *argv[]) {
  char od_fname[MAXSTRING], vel_fname[MAXSTRING], rock_fname[MAXSTRING];
  int wall, iolet;
  int nx, ny, nz;
  float ***vel, ***od, ***rock;
  float od_sum, od_sum_z;
  float vel_sum, vel_sum_z;
  float mf_sum, mf_sum_z;
  int sites, rocksites, rocksites_z, fluidsites, fluidsites_z;
  bool verbose = true;
  bool dumppermeability = false;
  bool repeat = true;
  int irepeat;

  float avg_grad_p_S, avg_rho_top_P, avg_rho_bottom_P, avg_rho_top, avg_rho_bottom, avg_rho_PS;
  float avg_v_S;
  float nu, eta;

  float tau, delta_x, delta_t;

  float perm, poro;
//HDF variable declarations
  hid_t           od_file, vel_file, rock_file, od_space, vel_space, rock_space, dset_od, dset_vel, dset_rock; // Handles
  herr_t          status;
  hsize_t         dims_od[RANK],
                  dims_vel[RANK],
                  dims_rock[RANK],
                  dims[RANK],
                  offset[RANK] = {0, 0, 0},
                  count[RANK] = {1, 1, 1};

  if (argc!=4 && argc!=5) {
    fprintf(stderr,"Syntax: %s <od_filename> <vel_filename> <rock_filename> [--noverbose/--dumppermeability]\n", argv[0]);
    exit(1);
  }

  sprintf(od_fname,"%s",argv[1]);
  sprintf(vel_fname,"%s",argv[2]);
  sprintf(rock_fname,"%s",argv[3]);
  
  if (argc == 5) {
    if ( strstr(argv[4],"--noverbose") != NULL ) {
      verbose = false;
    }
    if ( strstr(argv[4],"--dumppermeability") != NULL ) {
      verbose = false;
      dumppermeability = true;
    }

  }

  if (!dumppermeability) fprintf(stdout,"=== ROck CALculator  === \n");

  if (verbose) fprintf(stdout,"Attempting to read OD file %s \n", od_fname);
  od_file = H5Fopen (od_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (od_file < 0) {
    fprintf(stderr,"Invalid OD file, aborting \n");
    exit(1);
  }

  dset_od = H5Dopen(od_file, DATASET, H5P_DEFAULT);
  od_space = H5Dget_space(dset_od);
  H5Sget_simple_extent_dims(od_space,dims_od,NULL);
  if (verbose) fprintf (stdout,"OD dataspace dimensions are %d x %d x %d \n", (int) dims_od[2], (int) dims_od[1], (int) dims_od[0]);

  if (verbose) fprintf(stdout,"Attempting to read VEL file %s \n", vel_fname);
  vel_file = H5Fopen (vel_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (vel_file < 0) {
    fprintf(stderr,"Invalid VEL file, aborting \n");
    exit(1);
  }

  dset_vel = H5Dopen(vel_file, DATASET, H5P_DEFAULT);
  vel_space = H5Dget_space(dset_vel);
  H5Sget_simple_extent_dims(vel_space,dims_vel,NULL);
  if (verbose) fprintf (stdout,"VEL dataspace dimensions are %d x %d x %d \n", (int) dims_vel[2], (int) dims_vel[1], (int) dims_vel[0]);

  if (verbose) fprintf(stdout,"Attempting to read ROCK file %s \n", rock_fname);
  rock_file = H5Fopen (rock_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (rock_file < 0) {
    fprintf(stderr,"Invalid ROCK file, aborting \n");
    exit(1);
  }

  dset_rock = H5Dopen(rock_file, DATASET, H5P_DEFAULT);
  rock_space = H5Dget_space(dset_rock);
  H5Sget_simple_extent_dims(rock_space,dims_rock,NULL);
  if (verbose) fprintf (stdout,"ROCK dataspace dimensions are %d x %d x %d \n", (int) dims_rock[2], (int) dims_rock[1], (int) dims_rock[0]);
  
  
  if ( dims_od[0] != dims_vel[0] ||
       dims_od[1] != dims_vel[1] ||
       dims_od[2] != dims_vel[2] ||
       dims_od[0] != dims_rock[0] ||
       dims_od[1] != dims_rock[1] ||
       dims_od[2] != dims_rock[2]    ) {
    fprintf(stderr,"Dimensions of OD, VEL and ROCK files don't match, aborting\n");
    exit(1);
  }
  else {
    dims[0] = dims_od[0];
    dims[1] = dims_od[1];
    dims[2] = dims_od[2];
  }

  /* Check sizes - ask for offsets */

#ifndef WALL
  if (verbose) fprintf(stdout,"Enter wall thickness: \n");
  scanf("%d",&wall);
#else
  wall = WALL;
#endif

#ifndef IOLET
  if (verbose) fprintf(stdout,"Enter inlet/outlet size: \n");
  scanf("%d",&iolet);
#else
  iolet = IOLET;
#endif

  nx = dims[2];
  ny = dims[1];
  nz = dims[0];

  if (verbose) fprintf(stdout,"Actual system size will be %d x %d x %d\n", nx-2*wall, ny-2*wall, nz-2*iolet);

  /* Allocate tensor (i.e. 3d array) and set pointers */
  if (verbose) fprintf(stdout,"Setting up 3d arrays... \n");

  if (verbose) fprintf(stdout,"\t - OD... \n");
  od = (float ***) malloc(nz*sizeof(float **));
  if (!od) {
    fprintf(stderr,"Error: (1) can't allocate OD memory, aborting... \n");
    exit(1);
  }
  od[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!od[0]) { 
    fprintf(stderr,"Error: (2) can't allocate OD memory, aborting... \n");
    exit(1);
  }
  od[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!od[0][0]) {
    fprintf(stderr,"Error: (3) can't allocate OD memory, aborting... \n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    od[0][j] = od[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    od[k]=od[k-1]+ny;
    od[k][0]=od[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      od[k][j]=od[k][j-1]+nx;
    }
  }

//   od = (float ***) malloc(nz*sizeof(float **));
//   if (!od) {
//     fprintf(stderr,"Error: (1) can't allocate OD memory, aborting... \n");
//     exit(1);
//   }
//   for (int k=0;k<nz;k++) {
//     od[k] = (float **) malloc(ny*sizeof(float *));
//     if (!od[k]) {
//       fprintf(stderr,"Error: (1) can't allocate OD memory, aborting... \n");
//       exit(1);
//     }
//     for (int j=0;j<ny;j++) {
//       od[k][j] = (float *) malloc(nx*sizeof(float));
//       if (!od[k][j]) {
//         fprintf(stderr,"Error: (1) can't allocate OD memory, aborting... \n");
//         exit(1);
//       }
//     }
//   }

  /* Allocate tensor (i.e. 3d array) and set pointers */
  if (verbose) fprintf(stdout,"\t - VEL... \n");
  vel = (float ***) malloc(nz*sizeof(float **));
  if (!vel) {
    fprintf(stderr,"Error: (1) can't allocate VEL memory, aborting... \n");
    exit(1);
  }
  vel[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!vel[0]) { 
    fprintf(stderr,"Error: (2) can't allocate VEL memory, aborting... \n");
    exit(1);
  }
  vel[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!vel[0][0]) {
    fprintf(stderr,"Error: (3) can't allocate VEL memory, aborting... \n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    vel[0][j] = vel[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    vel[k]=vel[k-1]+ny;
    vel[k][0]=vel[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      vel[k][j]=vel[k][j-1]+nx;
    }
  }

  /* Allocate tensor (i.e. 3d array) and set pointers */
  if (verbose) fprintf(stdout,"\t - ROCK... \n");
  rock = (float ***) malloc(nz*sizeof(float **));
  if (!rock) {
    fprintf(stderr,"Error: (1) can't allocate ROCK memory, aborting... \n");
    exit(1);
  }
  rock[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!rock[0]) { 
    fprintf(stderr,"Error: (2) can't allocate ROCK memory, aborting... \n");
    exit(1);
  }
  rock[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!rock[0][0]) {
    fprintf(stderr,"Error: (3) can't allocate ROCK memory, aborting... \n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    rock[0][j] = rock[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    rock[k]=rock[k-1]+ny;
    rock[k][0]=rock[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      rock[k][j]=rock[k][j-1]+nx;
    }
  }

  if (verbose) fprintf(stdout,"Finished setting up 3d arrays \n");

  if (verbose) fprintf(stdout,"Loading OD into array\n");
  status = H5Sselect_hyperslab(od_space, H5S_SELECT_SET, offset, NULL, count, dims);
  status = H5Dread(dset_od, H5T_NATIVE_FLOAT, H5S_ALL, od_space, H5P_DEFAULT, od[0][0]);

  if (verbose) fprintf(stdout,"Loading VEL into array\n");
  status = H5Sselect_hyperslab(vel_space, H5S_SELECT_SET, offset, NULL, count, dims);
  status = H5Dread(dset_vel, H5T_NATIVE_FLOAT, H5S_ALL, vel_space, H5P_DEFAULT, vel[0][0]);

  if (verbose) fprintf(stdout,"Loading ROCK into array\n");
  status = H5Sselect_hyperslab(rock_space, H5S_SELECT_SET, offset, NULL, count, dims);
  status = H5Dread(dset_rock, H5T_NATIVE_FLOAT, H5S_ALL, rock_space, H5P_DEFAULT, rock[0][0]);

  /*
  fprintf(stdout,"Checking iolets...\n");
  if (iolet > 0) {
    float rocksum_z;
    for(int k=0;k<nz;k++) {
      if (k==iolet) k=nz-iolet;
      rocksum_z = 0.0;
      for(int j=wall;j<ny-wall;j++) {
        for(int i=wall;i<nx-wall;i++) {
          rocksum_z += rock[k][j][i];
        }
      }
      if (rocksum_z > 0.0) printf("WARNING: Rocksum at z = %d: %f != 0.0\n",k,rocksum_z);
    }
  }
  fprintf(stdout,"Done checking iolets...\n");
  */

  if (verbose) fprintf(stdout,"Starting measurements, looping over sites...\n");

  vel_sum = 0.0;
  od_sum = 0.0;
  mf_sum = 0.0;
  rocksites = 0;
  fluidsites = 0;
  sites = 0;

  if (verbose) {
    printf("[1]: z plane index\n");
    printf("[2]: Average porespace velocity on plane z\n");
    printf("[3]: Average porespace density on plane z\n");
    printf("[4]: Average porespace massflow =avg(rho(x,y)*v(x,y)) on plane z\n");
    printf("[5]: Plane porosity (porespace/x*y, walls not included) \n");
    printf("[1]    \t[2]    \t[3]    \t[4]    \t[5]\n");
  }

  for(int k=0;k<iolet;k++) {
    vel_sum_z = 0.0;
    od_sum_z = 0.0;
    mf_sum_z = 0.0;
    rocksites_z = 0;
    fluidsites_z = 0;

    for(int j=wall;j<ny-wall;j++) {
      for(int i=wall;i<nx-wall;i++) {
        //if (verbose) printf("Considering site (%d, %d, %d)\n",i,j,k);
        if (rock[k][j][i] > 0.0) {
          rocksites_z++;
        }
        else {
          fluidsites_z++;
          vel_sum_z += vel[k][j][i];
          od_sum_z += od[k][j][i];
          mf_sum_z += vel[k][j][i] * od[k][j][i];
        }
      }
    }
    if (verbose) printf("%d \t %d \t %d \t %e \t %e \t %e \t %e\n",k, rocksites_z, fluidsites_z, vel_sum_z/(float) fluidsites_z, od_sum_z/(float) fluidsites_z, mf_sum_z/(float) fluidsites_z, (float) fluidsites_z / (float) (fluidsites_z + rocksites_z));
  }


  for(int k=iolet;k<nz-iolet;k++) {
    vel_sum_z = 0.0;
    od_sum_z = 0.0;
    mf_sum_z = 0.0;
    rocksites_z = 0;
    fluidsites_z = 0;

    for(int j=wall;j<ny-wall;j++) {
      for(int i=wall;i<nx-wall;i++) {
        sites++;
        if (rock[k][j][i] > 0.0) {
          rocksites++;
          rocksites_z++;
        }
        else {
          fluidsites++;
          fluidsites_z++;
        }
        vel_sum_z += vel[k][j][i];
        vel_sum += vel[k][j][i];
        od_sum_z += od[k][j][i];
        od_sum += od[k][j][i];
        mf_sum_z += vel[k][j][i] * od[k][j][i];
        mf_sum += vel[k][j][i] * od[k][j][i];
      }
    }

    if (verbose) printf("%d \t %d \t %d \t %e \t %e \t %e \t %e\n",k, rocksites_z, fluidsites_z, vel_sum_z/(float) fluidsites_z, od_sum_z/(float) fluidsites_z, mf_sum_z/(float) fluidsites_z, (float) fluidsites_z / (float) (fluidsites_z + rocksites_z));
    
    if (k==iolet) {
      // z = LI+1
      avg_rho_bottom_P = od_sum_z/(float) fluidsites_z;
      if (verbose) printf("Setting avg_rho_bottom_P = %f \n",avg_rho_bottom_P);
      avg_rho_bottom = od_sum_z/(float) (fluidsites_z + rocksites_z);
      if (verbose) printf("Setting avg_rho_bottom = %f \n",avg_rho_bottom);
    }

    if (k==nz-iolet-1) {
      // z = L3-LO
      avg_rho_top_P = od_sum_z/(float) fluidsites_z;
      if (verbose) printf("Setting avg_rho_top_P = %f \n",avg_rho_top_P);
      avg_rho_top = od_sum_z/(float) (fluidsites_z + rocksites_z);
      if (verbose) printf("Setting avg_rho_top = %f \n",avg_rho_top);
    }

  }
  
  for(int k=nz-iolet;k<nz;k++) {
    vel_sum_z = 0.0;
    od_sum_z = 0.0;
    mf_sum_z = 0.0;
    rocksites_z = 0;
    fluidsites_z = 0;

    for(int j=wall;j<ny-wall;j++) {
      for(int i=wall;i<nx-wall;i++) {
        if (rock[k][j][i] > 0.0) {
          rocksites_z++;
        }
        else {
          fluidsites_z++;
          vel_sum_z += vel[k][j][i];
          od_sum_z += od[k][j][i];
          mf_sum_z += vel[k][j][i] * od[k][j][i];
        }
      }
    }
    if (verbose) printf("%d \t %d \t %d \t %e \t %e \t %e \t %e\n",k, rocksites_z, fluidsites_z, vel_sum_z/(float) fluidsites_z, od_sum_z/(float) fluidsites_z, mf_sum_z/(float) fluidsites_z, (float) fluidsites_z / (float) (fluidsites_z + rocksites_z));
  }

  if (avg_rho_top_P > avg_rho_bottom_P) {
    if (verbose) fprintf(stdout,"WARNING: pressure is higher at the top than at the bottom.\n");
  }

  if (sites != (nx-2*wall) * (ny-2*wall) * (nz-2*iolet)) {
    fprintf(stderr,"ERROR: total sites counted unequal to total sites calculated\n");
    exit(1);
  }

  if (verbose) printf("Finished looping over sites\n");

  /* Close all the HDF5 handles */
  status = H5Dclose (dset_od);
  status = H5Dclose (dset_vel);
  status = H5Dclose (dset_rock);
  status = H5Sclose (od_space);
  status = H5Sclose (vel_space);
  status = H5Sclose (rock_space);
  status = H5Fclose (od_file);
  status = H5Fclose (vel_file);
  status = H5Fclose (rock_file);

  while (repeat) {

#ifndef TAU
    if (verbose) fprintf(stdout,"Enter tau (s): \n");
    scanf("%f",&tau);
#else
    tau = TAU;
#endif

#ifndef DELTA_X
    if (verbose) fprintf(stdout,"Enter delta_x (m): \n");
    scanf("%f",&delta_x);
#else
    delta_x = DELTA_X;
#endif

#ifndef DELTA_T
    if (verbose) fprintf(stdout,"Enter delta_t (s): \n");
    scanf("%f",&delta_t);
#else
    delta_t = DELTA_T;
#endif

    if (!dumppermeability) printf("Sample size (%d x %d x %d): %d sites, %d rock sites, %d fluid sites \n", nx-2*wall, ny-2*wall, nz-2*iolet, sites, rocksites, fluidsites);
    if (!dumppermeability) printf("Tau = %e, delta_x = %e, delta_t = %e \n", tau, delta_x, delta_t);

    poro = (float) fluidsites / (float) sites;
    if (!dumppermeability) printf("\tCalculated porosity: %e\n", poro);

    // avg_grad_p_S should be < 0 for g_accn > 0
    avg_grad_p_S = (delta_x * (avg_rho_top_P - avg_rho_bottom_P) ) / (3.0 * delta_t * delta_t * (float) (nz-2*iolet - 1 ) ); // (Ls - 1)
    if (verbose) printf("Setting avg_grad_p_S = %e \n",avg_grad_p_S);
    avg_v_S = vel_sum / (float) sites;
    if (verbose) printf("Setting avg_v_S = %e \n",avg_v_S);
#ifdef USETOTALAVGRHO
    avg_rho_PS = od_sum / (float) fluidsites;
#else
    avg_rho_PS = 0.5*(avg_rho_top_P + avg_rho_bottom_P);
#endif
    if (verbose) printf("Setting avg_rho_PS = %e \n",avg_rho_PS);
    nu = ( (tau / delta_t) - 0.5 ) * (delta_x * delta_x) / (3.0 * delta_t);
    if (verbose) printf("Setting nu = %e\n",nu);
    eta = nu * avg_rho_PS;
    if (verbose) printf("Setting eta = %e\n",eta);
    perm = -eta * avg_v_S / avg_grad_p_S;
    if (!dumppermeability) { 
      printf("\tCalculated permeability ( (delta x)^2 ): %e\n", perm);
      printf("\tCalculated permeability (mD): %e\n", 1000.0*perm/DARCYM2);
    }
    else {
      printf("%e %e\n",perm, 1000.0*perm/DARCYM2);
    }

    if (verbose) fprintf(stdout,"Done! \n");

#ifdef NOREPEAT
    repeat = false;
#endif

    if (verbose && repeat) {
      fprintf(stdout,"Change parameters (0/1)? \n");
      scanf("%d",&irepeat);
      if (irepeat <= 0) repeat = false;
    }
  }

  return 0;
}
