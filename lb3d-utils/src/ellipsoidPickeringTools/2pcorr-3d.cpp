#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

int main (int argc, char *argv[]) {

  int Gx =  0; // System size X
  int Gy =  0; // System size Y
  int Gz =  0; // System size Z
  int NP = -1; // Number of particles
  int t  = -1; // Timestep
 
  // File name strings
  char* fname;

  // Flags
  bool fname_set = false;

  // Parsing CLI
  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-f") == 0 ) {
      if ( i+1 < argc ) {
	fname = argv[++i];
	//fprintf(stdout,"  MD file <%s>\n",fname);
	fname_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -f flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-n") == 0 ) {
      if ( i+1 < argc ) {
	NP = atoi(argv[++i]);
	//fprintf(stdout,"  Number of particles <%d>\n",NP);
      }
      else {
	fprintf(stderr,"  Missing argument to -n flag. \n");
	exit(1);
      }
    }

    else if ( strcmp(argv[i],"-t") == 0 ) {
      if ( i+1 < argc ) {
	t = atoi(argv[++i]);
	//fprintf(stdout,"  Timestep <%d>\n",t);
      }
      else {
	fprintf(stderr,"  Missing argument to -t flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-x") == 0 ) {
      if ( i+1 < argc ) {
	Gx = atoi(argv[++i]);
	//fprintf(stdout,"  nx <%d>\n",Gx);
      }
      else {
	fprintf(stderr,"  Missing argument to -x flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-y") == 0 ) {
      if ( i+1 < argc ) {
	Gy = atoi(argv[++i]);
	//fprintf(stdout,"  ny <%d>\n",Gy);
      }
      else {
	fprintf(stderr,"  Missing argument to -y flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-z") == 0 ) {
      if ( i+1 < argc ) {
	Gz = atoi(argv[++i]);
	//fprintf(stdout,"  nz <%d>\n",Gz);
      }
      else {
	fprintf(stderr,"  Missing argument to -z flag. \n");
	exit(1);
      }
    }
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }
  }

  // Checking for mandatory parameters
  if (! fname_set ) {
    fprintf(stderr,"ERROR: MD file -f is mandatory.\n");
    fprintf(stderr,"  -f               <md-cfg file>\n");
    fprintf(stderr,"  -n               <number of particles>\n");
    fprintf(stderr,"  -t               <timestep>\n");
    fprintf(stderr,"  -x               <nx>\n");
    fprintf(stderr,"  -y               <ny>\n");
    fprintf(stderr,"  -z               <nz>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  // Arrays that will hold positons
  double X[NP];
  double Y[NP];
  double Z[NP];

  // Open file for reading
  FILE * pFile;     
  pFile = fopen(fname,"r");
  if (!pFile) {
    fprintf(stderr, "Failed to open file <%s>... \n",fname);
    exit(1);
  }

  // Read data
  for(int p = 0; p < NP; p++) {
    fscanf( pFile, "%le %le %le %*e %*e %*e %*e %*e %*e %*d\n", &X[p], &Y[p], &Z[p]);

    //fprintf(stdout,"NEW %d %lE %lE %lE\n", p,X[p],Y[p],Z[p]);
  }
                
  // We're done with this file
  fclose(pFile);

  fprintf(stdout,"File <%s> read succesfully.",fname);

  int    Ngr[NP][Gx/2]; // Number of particles in a shell around a particular particle
  double gr[NP][Gx/2];  // Correlations per particle
  double Gr[Gx/2];      // Correlations, binned for distance only
  double dNgr[Gx/2];

  // Zero out some arrays
  for (int p = 0; p < NP; p++) {
    for (int i = 0; i < Gx/2; i++) {
      Ngr[p][i]  = 0;
      gr[p][i]   = 0.0;
      Gr[i]      = 0.0;
      dNgr[i]    = 0.0;
    }
  }

  // Distances
  double Dx, Dy, Dz, Dr;

  // Comment this set of arrays...

  int iDr = 0;
  int Ngrp[NP][Gx/2];
  double grp[NP][Gx/2];
  double Grp[Gx/2];

  // Loop over all particles
  for (int p1 = 0; p1 < NP; p1++) {
    for (int p2 = 0; p2 < NP; p2++) {
      // Looping over all distinct particle pairs (double loop?)
      if (p1 != p2) {

	// Calculated boxed X-distance
	Dx = X[p2] - X[p1];
	if (Dx < 0) Dx = -Dx;
	if (Dx > Gx/2) Dx = Gx - Dx;

	// Calculated boxed Y-distance
	Dy = Y[p2] - Y[p1];
	if (Dy < 0) Dy = -Dy;
	if (Dy > Gy/2) Dy = Gy - Dy;

	// Calculated boxed Z-distance
	Dz = Z[p2] - Z[p1];
	if (Dz < 0) Dz = -Dz;
	if (Dz > Gz/2) Dz = Gz - Dz;

	// Calculated boxed total distance (shells of finite thickness)
	Dr = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
	iDr = int(Dr+0.5);

	// Increase number of particles in the annulus iDr	    
	if(Dr < Gx/2) Ngr[p1][iDr]++;
      }
    }
    // Loop over annuli
    for (int i = 1; i < Gx/2; i++) {
      // Rescale correlation function
      gr[p1][i] = Ngr[p1][i]/(2*M_PI*i);
      // Sum over all particles
      Gr[i] += gr[p1][i];

      dNgr[i] += Ngr[p1][i];
      grp[p1][i] = Ngrp[p1][i]/(2*M_PI*i);
      Grp[i] += grp[p1][i];
    }
  }

  int NdGrb    = 0;
  int NdGrpb   = 0;
  double dGrpb = 0.0;
  double dGrb  = 0.0;

  for(int i = 1; i < Gx/2; i++)	{
    // Rescale for number of particles
    Gr[i] = Gr[i]/NP;
    dNgr[i] = dNgr[i]/NP;
    Grp[i] = Grp[i]/NP;
    if(i > Gx/4) {
      // Count / sum totals for latter half of the range
      dGrb += Gr[i];
      NdGrb++;
      dGrpb += Grp[i];
      NdGrpb++;
    }
  }
  // Take averages of the later half of the range
  dGrb = dGrb/NdGrb;
  dGrpb = dGrpb/NdGrpb;
  
  // Declare arrays for weighted averages
  // Numbers denote different methods (see code below)
  // "n" stands for "normed"
 
  double Gr_avg5[Gx/2];
  double Gr_avg5n[Gx/2];
  
  // Taking weighted average over 7 sites, special cases exist for tail ends of the curve (truncation of the averaging)
  Gr_avg5[0] = (4*Gr[0]+3*Gr[1]+2*Gr[2]+Gr[3])/10;
  Gr_avg5[1] = (3*Gr[0]+4*Gr[1]+3*Gr[2]+2*Gr[3]+Gr[4])/13;
  Gr_avg5[2] = (2*Gr[0]+3*Gr[1]+4*Gr[2]+3*Gr[3]+2*Gr[4]+Gr[5])/15;;
  for(int i = 3; i<Gx/2 - 3; i++)	{
    Gr_avg5[i] = (Gr[i-3]+2*Gr[i-2]+3*Gr[i-1]+4*Gr[i]+3*Gr[i+1]+2*Gr[i+2]+Gr[i+3])/16;
    if(Gr[i] == 0) Gr_avg5[i] = 0.0;
  }
  Gr_avg5[Gx/2 - 3] = (Gr[Gx/2 - 6]+2*Gr[Gx/2 - 5]+3*Gr[Gx/2 - 4]+4*Gr[Gx/2 - 3]+3*Gr[Gx/2 - 2]+2*Gr[Gx/2 - 1])/15;
  Gr_avg5[Gx/2 - 2] = (Gr[Gx/2 - 5]+2*Gr[Gx/2 - 4]+3*Gr[Gx/2 - 3]+4*Gr[Gx/2 - 2]+3*Gr[Gx/2 - 1])/13;
  Gr_avg5[Gx/2 - 1] = (Gr[Gx/2 - 4]+2*Gr[Gx/2 - 3]+3*Gr[Gx/2 - 2]+4*Gr[Gx/2 - 1])/10;

  // Declare arrays for normed values
  double Grn[Gx/2];

  // Open file for writing G(r)
  FILE * pFile_Gr;     
  char pFile_Gr_fname[128];
  sprintf(pFile_Gr_fname,"Gr-t%.8d.asc",t);
  pFile_Gr = fopen(pFile_Gr_fname,"w");
  if (!pFile_Gr) {
    fprintf(stderr, "Failed to open file <%s>... \n",pFile_Gr_fname);
    exit(1);
  }

  // Open file for writing G_avg5(r)
  FILE * pFile_Gr5;     
  char pFile_Gr5_fname[128];
  sprintf(pFile_Gr5_fname,"Gr_avg5-t%.8d.asc",t);
  pFile_Gr5 = fopen(pFile_Gr5_fname,"w");
  if (!pFile_Gr5) {
    fprintf(stderr, "Failed to open file <%s>... \n",pFile_Gr5_fname);
    exit(1);
  }

  fprintf(pFile_Gr ,"#? %5s %8s %12s %12s %12s %12s\n","t","shell","G(r)","G_n(r)","G(r)*NP","dNgr");
  fprintf(pFile_Gr5,"#? %5s %8s %12s %12s\n","t","shell","G(r)","G_n(r)");

  for(int i = 1; i < Gx/2; i++)	{
    // Rescale by averages to get ~1 at the tails
    Grn[i] = Gr[i] / dGrb;
    Gr_avg5n[i] = Gr_avg5[i] / dGrb;
    fprintf(stdout   ,"%.8d %.8d %lE %lE %lE %lE\n", t, i, Gr[i] , Grn[i] , Gr[i]*NP , dNgr[i]);
    fprintf(pFile_Gr ,"%.8d %.8d %lE %lE %lE %lE\n", t, i, Gr[i] , Grn[i] , Gr[i]*NP , dNgr[i]);
    fprintf(pFile_Gr5,"%.8d %.8d %lE %lE\n", t, i, Gr_avg5[i], Gr_avg5n[i]);
  }
 
  fclose(pFile_Gr);
  fclose(pFile_Gr5);

  return 0;

} 
