#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>

#include "droplet-surfacetension.h"
#include "hdf5_helper.h"
#include "lbe_version.h"

using namespace std;

double psi(double rho) {
  return (1.0 - pow(M_E,-rho));
}

double calc_r_mass(double od_in, double od_corner, int np, double Rp, double M_d) {
  double rpow3 = M_d / ( 4.0/3.0 * M_PI * ( od_in - od_corner ) );
  rpow3 += 0.5*np*pow(Rp,3);
  return pow(rpow3,1.0/3.0);
}

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz;
  int cx, cy, cz;

  int avg = 2; // Average over how many lattice sites in each direction?

  double p_in,   p_corner;
  double pc_in,  pc_corner;
  double rho_in, rho_corner;
  double od_in,  od_corner,  od_total;
  double wd_in,  wd_corner,  wd_total;
  double sur_in, sur_corner, sur_total;

  int np_sphere = 0;
  int np_sphere_cap = 0;
  int np = 0; // Number of particles
  double Rp = 0.0;
  double g_br = 0.0;
  double g_bs = 0.0;
  double g_ss = 0.0;

  double r_mass;

  // File names
  char *fname_od, *fname_wd, *fname_sur;

  // Array pointers
  float ***od, ***wd, ***sur;

  hsize_t *dims;

  int t = 0;

  // Flags
  bool fname_od_set = false;
  bool fname_sur_set = false;
  bool fname_wd_set = false;
  bool header = false;
  // bool verbose = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);
    if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--header") == 0    ) {
      header = true;
      //fprintf(stdout,"  Header\n");
    }
    else if ( strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--odfile") == 0     ) {
      if ( i+1 < argc ) {
	fname_od = argv[++i];
	//fprintf(stdout,"  OD file <%s>\n",fname_od);
	fname_od_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -o flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--surfile") == 0    ) {
      if ( i+1 < argc ) {
	fname_sur = argv[++i];
	//fprintf(stdout,"  Surfactant file <%s>\n",fname_sur);
	fname_sur_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -s flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--timestep") == 0   ) {
      if ( i+1 < argc ) {
	t = atoi(argv[++i]);
	//fprintf(stdout,"  Time step <%d>\n",t);
      }
      else {
	fprintf(stderr,"  Missing argument to -t flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-w") == 0 || strcmp(argv[i],"--wdfile") == 0     ) {
      if ( i+1 < argc ) {
	fname_wd = argv[++i];
	//fprintf(stdout,"  WD file <%s>\n",fname_wd);
	fname_wd_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -w flag. \n");
	exit(1);
      }

    }
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }

  }

  if (! ( fname_od_set && fname_wd_set) ) {
    fprintf(stdout,"droplet-surfacetension (%s)\n\n", GIT_DESC);
    fprintf(stderr,"od file -o, wd file -w are mandatory.\n\n");
    fprintf(stderr,"  -h --header      \n");
    fprintf(stderr,"  -o --odfile      <odfile>\n");
    fprintf(stderr,"  -s --surfile     <surfile>\n");
    fprintf(stderr,"  -t --timestep    <timestep>\n");
    fprintf(stderr,"  -w --wdfile      <wdfile>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  if (hdfReadAll(fname_od    ,&od    ,&dims) != 0) {
    exit(1);
  }
  hdfFreeDims(dims);
  if (hdfReadAll(fname_wd    ,&wd    ,&dims) != 0) {
    exit(1);
  }
  if (fname_sur_set) {
    hdfFreeDims(dims);
    if (hdfReadAll(fname_sur   ,&sur   ,&dims) != 0) {
      exit(1);
    }
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  // Get attributes
  LB3dAttributes attr = LB3dAttributes(fname_od);
  attr.getDouble("g_br",&g_br);
  attr.getDouble("g_bs",&g_bs);
  attr.getDouble("g_ss",&g_ss);
  attr.getInteger("np_sphere",&np_sphere);
  attr.getInteger("np_sphere_cap",&np_sphere_cap);
  attr.getDouble("R_para",&Rp);
  np = np_sphere - 2*np_sphere_cap;

  // Calculate center positions
  cx = nx/2 - 1;
  cy = ny/2 - 1;
  cz = nz/2 - 1;

  // Calculate densities

  rho_in = 0.0;
  od_in  = 0.0;
  wd_in  = 0.0;
  sur_in = 0.0;
  for(int k=-avg;k<=avg;k++) {
    for(int j=-avg;j<=avg;j++) {
      for(int i=-avg;i<=avg;i++) {
        rho_in += od[cz+k][cy+j][cx+i];
        rho_in += wd[cz+k][cy+j][cx+i];
        if (fname_sur_set) {
          rho_in += sur[cz+k][cy+j][cx+i];
          sur_in += sur[cz+k][cy+j][cx+i];
        }
        od_in += od[cz+k][cy+j][cx+i];
        wd_in += wd[cz+k][cy+j][cx+i];
      }
    }
  }
  rho_in /= pow((2*avg + 1), 3);
  p_in    = rho_in / 3.0; // speed of sound ^2
  od_in  /= pow((2*avg + 1), 3);
  wd_in  /= pow((2*avg + 1), 3);
  sur_in /= pow((2*avg + 1), 3);

  //fprintf(stdout,"P_in: %f\n",pin);

  rho_corner = 0.0;
  od_corner  = 0.0;
  wd_corner  = 0.0;
  sur_corner = 0.0;
  for(int k=0;k<=2*avg;k++) {
    for(int j=0;j<=2*avg;j++) {
      for(int i=0;i<=2*avg;i++) {
        rho_corner += od[k][j][i];
        rho_corner += wd[k][j][i];
        if (fname_sur_set) {
          rho_corner += sur[k][j][i];
          sur_corner += sur[k][j][i];
        }
        od_corner += od[k][j][i];
        wd_corner += wd[k][j][i];
      }
    }
  }
  rho_corner /= pow((2*avg + 1), 3);
  p_corner    = rho_corner / 3.0; // speed of sound ^2
  od_corner  /= pow((2*avg + 1), 3);
  wd_corner  /= pow((2*avg + 1), 3);
  sur_corner /= pow((2*avg + 1), 3);

  od_total  = 0.0;
  wd_total  = 0.0;
  sur_total = 0.0;
  for(int k=0;k<nz;k++) {
    for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
        od_total  += od[k][j][i];
        wd_total  += wd[k][j][i];
	if (fname_sur_set) {
	  sur_total += sur[k][j][i];
	}
      }
    }
  }

  // Calculate r thru mass bunching

  pc_in     = p_in     + 18.0*2.0*(g_br*psi(od_in    )*psi(wd_in    ))/3.0 + 18.0*2.0*(g_bs*psi(od_in    )*psi(sur_in    ))/3.0 + 18.0*2.0*(g_bs*psi(wd_in    )*psi(sur_in    ))/3.0 + 18.0*(g_ss*psi(sur_in    )*psi(sur_in    ))/3.0;
  pc_corner = p_corner + 18.0*2.0*(g_br*psi(od_corner)*psi(wd_corner))/3.0 + 18.0*2.0*(g_bs*psi(od_corner)*psi(sur_corner))/3.0 + 18.0*2.0*(g_bs*psi(wd_corner)*psi(sur_corner))/3.0 + 18.0*(g_ss*psi(sur_corner)*psi(sur_corner))/3.0;

  double M_d = (od_total - nx*ny*nz*od_corner);

  r_mass = calc_r_mass(od_in, od_corner, np, Rp, M_d);

  double DPc = pc_in - pc_corner;
  double sigma = 0.5*DPc*r_mass;
  double rho_m = wd_corner;

  double chi = ( (double) np ) * ( r_mass - sqrt( ( r_mass * r_mass ) - ( Rp * Rp ) ) ) / (2.0* r_mass) ;

  if (fname_sur_set) {
    rho_m += sur_corner;
  }

  if (header) {
    fprintf(stdout,"#\n");
    fprintf(stdout,"# GENERATED BY: droplet-surfacetension (%s) \n", GIT_DESC);
    fprintf(stdout,"#\n");
    fprintf(stdout,"# g_br = %f\n",g_br);
    fprintf(stdout,"# g_bs = %f\n",g_bs);
    fprintf(stdout,"# g_ss = %f\n",g_ss);
    fprintf(stdout,"# np_sphere = %d , np_sphere_cap = %d, np = %d\n",np_sphere,np_sphere_cap,np);
    fprintf(stdout,"# Rp = %f\n",Rp);
    fprintf(stdout,"#\n");

    fprintf(stdout,"#? %5.5s %12.12s %12.12s %12.12s %12.12s %12.12s %10.10s %10.10s %10.10s %12.12s %12.12s %12.12s\n",
	    "t","R","Pc_in","Pc_corner","DPc","sigma","g_br","g_bs","g_ss","M_d","rho_m","chi");
    fprintf(stdout,"# ===============================================================================================================================================\n");
  }

  fprintf(stdout,"%#8.8d %E %E %E %E %E % 9.3E % 9.3E % 9.3E %E %E %E\n",t,r_mass,pc_in,pc_corner,DPc,sigma,g_br,g_bs,g_ss,M_d,rho_m,chi);

  // fprintf(stdout,"# %#8.8d wd_total = %E , wd_in = %E , wd_corner = %E \n", t, wd_total, wd_in, wd_corner);
  // fprintf(stdout,"# %#8.8d od_total = %E , od_in = %E , od_corner = %E => droplet mass = %E , droplet excess density = %E\n",t, od_total, od_in, od_corner, M_d, (od_in-od_corner));

  //fprintf(stdout,"Freeing od!\n");
  hdfFree3DArray(od);
  //fprintf(stdout,"Freeing wd!\n");
  hdfFree3DArray(wd);
  if (fname_sur_set) {
    //fprintf(stdout,"Freeing sur!\n");
    hdfFree3DArray(sur);
  }

  return 0;
}
