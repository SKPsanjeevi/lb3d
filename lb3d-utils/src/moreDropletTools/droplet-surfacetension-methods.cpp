#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>

#ifndef NOMKL
#include "droplet_fit.h"
#else
#include "hdf5_helper.h"
#endif

#include "droplet-surfacetension-methods.h"

using namespace std;

double psi(double rho) {
  return (1.0 - pow(M_E,-rho));
}

double calc_r_mass(double od_in, double od_top, double od_corner, double od_total, int nx, int ny, int nz) {
  double rpow3_corner = (od_total - nx*ny*nz*od_corner) / ( 4.0/3.0 * M_PI * ( od_in - od_corner ));
  double rpow3_top    = (od_total - nx*ny*nz*od_top   ) / ( 4.0/3.0 * M_PI * ( od_in - od_top    ));

  // Dummy operation on rpow3_top to prevent compiler warnings. Bad me!
  rpow3_top /= 1.0;

  return pow(rpow3_corner,1.0/3.0);
}

#ifndef NOMKL
double calc_r_fit(char *fname_colour) {
  EllipseFit fitX = getEllipseFitFromSlice(fname_colour,X);
  return fitX.getLength(); // Calculate radius for a sphere
}
#endif

double calc_r_lin(float ***colour, int cx, int cy, int nz) {
  bool indroplet = false;
  double d_lin = 0.0;
  float c;
  
  // Linear interpolation between the two sides at the edge of the droplet.
  for (int k=0;k<nz;k++) {
    c = colour[k][cy][cx];
    if (indroplet) {
      if (c < 0.0) { // Going from red (drop) to blue (bulk)
	d_lin += ( colour[k-1][cy][cx] / ( colour[k-1][cy][cx] - c ) );
	indroplet = false;
      }
      else {
      d_lin += 1.0;
      }
    }
    else {
      if (c > 0.0) { // Going from blue (bulk) to red (drop)
	d_lin += ( c / ( c - colour[k-1][cy][cx] ) );
	indroplet = true;
      }
    }
  }
  
  // Radius = diameter / 2.0
 
  return (0.5*d_lin);
}

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz;
  int cx, cy, cz;

  int avg = 2;

  double r_fit = 0.0,  r_lin,   r_mass;
  double p_in,   p_top,   p_corner;
  double pc_in,  pc_top,  pc_corner;
  double scp_in, scp_top, scp_corner;
  double rho_in, rho_top, rho_corner;
  double od_in,   od_top,   od_corner,  od_total;
  double wd_in,   wd_top,   wd_corner,  wd_total;
  double sur_in,  sur_top,  sur_corner, sur_total;
  double g_br = 0.0;
  double g_bs = 0.0;
  double g_ss = 0.0;

  // File names
  char *fname_colour, *fname_od, *fname_wd, *fname_sur, *fname_scp;

  // Array pointers
  float ***colour, ***od, ***wd, ***sur, ***scp;

  hsize_t *dims;

  int t = 0;

  // Flags
  bool fname_colour_set = false;
  bool fname_od_set = false;
  bool fname_scp_set = false;
  bool fname_sur_set = false;
  bool fname_wd_set = false;
  bool profile = false;
  bool header = false;
  bool verbose = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if      ( strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--colourfile") == 0 ) {
      if ( i+1 < argc ) {
	fname_colour = argv[++i];
	//fprintf(stdout,"  Colour file <%s>\n",fname_colour);
	fname_colour_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -c flag. \n");
	exit(1);
      }
       }
    else if ( strcmp(argv[i],"-g") == 0 ||  strcmp(argv[i],"--g_br") == 0 ) {
      if ( i+1 < argc ) {
	g_br = atof(argv[++i]);
	//fprintf(stdout,"  g_br <%f>\n",g_br);
      }
      else {
	fprintf(stderr,"  Missing argument to -g flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-g_bs") == 0 )  {
      if ( i+1 < argc ) {
	g_bs = atof(argv[++i]);
	//fprintf(stdout,"  g_bs <%f>\n",g_bs);
      }
      else {
	fprintf(stderr,"  Missing argument to -g_bs flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-g_ss") == 0 )  {
      if ( i+1 < argc ) {
	g_ss = atof(argv[++i]);
	//fprintf(stdout,"  g_ss <%f>\n",g_ss);
      }
      else {
	fprintf(stderr,"  Missing argument to -g_ss flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--header") == 0    ) {
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
    else if ( strcmp(argv[i],"-p") == 0 || strcmp(argv[i],"--profile") == 0    ) {
      profile = true;
      //fprintf(stdout,"  Profile\n");
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
    else if ( strcmp(argv[i],"-scp") == 0 ) {
      if ( i+1 < argc ) {
	fname_scp = argv[++i];
	//fprintf(stdout,"  SCP file <%s>\n",fname_scp);
	fname_scp_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -scp flag. \n");
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
    else if ( strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0    ) {
      verbose = true;
      //fprintf(stdout,"  Verbose\n");
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

  if (! (fname_colour_set && fname_od_set && fname_wd_set) ) {
    fprintf(stderr,"Colour file -c, od file -o, wd file -w are mandatory.\n\n");
    fprintf(stderr,"  -c --colourfile  <colourfile>\n");
    fprintf(stderr,"  -g --g_br        <g_br>\n");
    fprintf(stderr,"  -g_bs            <g_bs>\n");
    fprintf(stderr,"  -g_ss            <g_ss>\n");
    fprintf(stderr,"  -h --header      \n");
    fprintf(stderr,"  -o --odfile      <odfile>\n");
    fprintf(stderr,"  -p --profile     \n");
    fprintf(stderr,"  -s --surfile     <surfile>\n");
    fprintf(stderr,"  -scp             <scpfile>\n");
    fprintf(stderr,"  -t --timestep    <timestep>\n");
    fprintf(stderr,"  -v --verbose     \n");
    fprintf(stderr,"  -w --wdfile      <wdfile>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  if (hdfReadAll(fname_colour,&colour,&dims) != 0) {
    exit(1);
  }
  hdfFreeDims(dims);
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
  if (fname_scp_set) {
    hdfFreeDims(dims);
    if (hdfReadAll(fname_scp   ,&scp   ,&dims) != 0) {
      exit(1);
    }
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  // Calculate center positions
  cx = nx/2 - 1;
  cy = ny/2 - 1;
  cz = nz/2 - 1;

  // Calculate r thru phi == 0 fit

#ifndef NOMKL
  r_fit = calc_r_fit(fname_colour);
#endif

  // Calculate r thru linear interpolation

  r_lin = calc_r_lin(colour,cx,cy,nz);

  // Calculate densities

  rho_in = 0.0;
  scp_in = 0.0;
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
	if (fname_scp_set) {
	  scp_in += scp[cz+k][cy+j][cx+i];
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
  scp_in /= pow((2*avg + 1), 3);

  //fprintf(stdout,"P_in: %f\n",pin);

  rho_corner = 0.0;
  scp_corner = 0.0;
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
        if (fname_scp_set) {
          scp_corner += scp[k][j][i];
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
  scp_corner /= pow((2*avg + 1), 3);


  rho_top = 0.0;
  scp_top = 0.0;
  od_top  = 0.0;
  wd_top  = 0.0;
  sur_top = 0.0;
  for(int k=0;k<=2*avg;k++) {
    for(int j=-avg;j<=avg;j++) {
      for(int i=-avg;i<=avg;i++) {
        rho_top += od[k][cy+j][cx+i];
        rho_top += wd[k][cy+j][cx+i];
        if (fname_sur_set) {
          rho_top += sur[k][cy+j][cx+i];
          sur_top += sur[k][cy+j][cx+i];
        }
        if (fname_scp_set) {
          scp_top += scp[k][cy+j][cx+i];
        }
        od_top += od[k][cy+j][cx+i];
        wd_top += wd[k][cy+j][cx+i];
      }
    }
  }
  rho_top /= pow((2*avg + 1), 3);
  p_top    = rho_top / 3.0; // speed of sound ^2
  od_top  /= pow((2*avg + 1), 3);
  wd_top  /= pow((2*avg + 1), 3);
  sur_top /= pow((2*avg + 1), 3);
  scp_top /= pow((2*avg + 1), 3);

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

  r_mass = calc_r_mass(od_in, od_top, od_corner, od_total, nx, ny, nz);


  // Calculate corrected pressures

  // + G cs^2 psi^2 / 2.0 = = g psi^2 /6.0

  // pc_in     = p_in     - (g_br*psi(rho_in    )*psi(rho_in    ))/6.0;
  // pc_top    = p_top    - (g_br*psi(rho_top   )*psi(rho_top   ))/6.0;
  // pc_corner = p_corner - (g_br*psi(rho_corner)*psi(rho_corner))/6.0;

  // +  cs^2 Sum{a != b} g_ab psi(a) psi(b)

  pc_in     = p_in     + 18.0*2.0*(g_br*psi(od_in    )*psi(wd_in    ))/3.0 + 18.0*2.0*(g_bs*psi(od_in    )*psi(sur_in    ))/3.0 + 18.0*2.0*(g_bs*psi(wd_in    )*psi(sur_in    ))/3.0 + 18.0*(g_ss*psi(sur_in    )*psi(sur_in    ))/3.0;
  pc_top    = p_top    + 18.0*2.0*(g_br*psi(od_top   )*psi(wd_top   ))/3.0 + 18.0*2.0*(g_bs*psi(od_top   )*psi(sur_top   ))/3.0 + 18.0*2.0*(g_bs*psi(wd_top   )*psi(sur_top   ))/3.0 + 18.0*(g_ss*psi(sur_top   )*psi(sur_top   ))/3.0;
  pc_corner = p_corner + 18.0*2.0*(g_br*psi(od_corner)*psi(wd_corner))/3.0 + 18.0*2.0*(g_bs*psi(od_corner)*psi(sur_corner))/3.0 + 18.0*2.0*(g_bs*psi(wd_corner)*psi(sur_corner))/3.0 + 18.0*(g_ss*psi(sur_corner)*psi(sur_corner))/3.0;

  double M_d = (od_total - nx*ny*nz*od_corner);

  if (profile) {
    if (header) {
      fprintf(stdout,"#\n");
      fprintf(stdout,"# GENERATED BY: droplet-surfacetension\n");
      fprintf(stdout,"# WARNING: NOSURFACTANT\n");
      fprintf(stdout,"#\n");
      fprintf(stdout,"#  z       od       wd     od+wd       phi\n");
    }
    for(int k=0;k<nz;k++) {
      fprintf(stdout,"%#4.4d %f %f %+f %+f\n",k ,od[k][cy][cx] ,wd[k][cy][cx] ,od[k][cy][cx]+wd[k][cy][cx], od[k][cy][cx]-wd[k][cy][cx]);
    }
  }
  else {
    if (header) {
      fprintf(stdout,"#\n");
      fprintf(stdout,"# GENERATED BY: droplet-surfacetension\n");
      fprintf(stdout,"#\n");
      fprintf(stdout,"#      t        R_fit        R_lin       R_mass         P_in        P_top     P_corner        Pc_in       Pc_top    Pc_corner         g_br         g_bs         g_ss          M_d\n");
      fprintf(stdout,"# ===============================================================================================================================================================================\n");
    }
    //fprintf(stdout,"%#8.8d %f %f %f %f %f %f %f %f %f %f %f %f\n",t,r_fit,r_lin,r_mass,p_in,p_top,p_corner,pc_in,pc_top,pc_corner,g_br,g_bs,g_ss);
    fprintf(stdout,"%#8.8d %E %E %E %E %E %E %E %E %E %E %E %E %E\n",t,r_fit,r_lin,r_mass,p_in,p_top,p_corner,pc_in,pc_top,pc_corner,g_br,g_bs,g_ss,M_d);
    if (verbose) {
      fprintf(stdout,"# %#8.8d wd_total = %E , wd_in = %E , wd_corner = %E \n", t, wd_total, wd_in, wd_corner);
      fprintf(stdout,"# %#8.8d od_total = %E , od_in = %E , od_corner = %E => droplet mass = %E , droplet excess density = %E\n",t, od_total, od_in, od_corner, M_d, (od_in-od_corner));
    }
  }

  //fprintf(stdout,"Freeing colour!\n");
  hdfFree3DArray(colour);
  //fprintf(stdout,"Freeing od!\n");
  hdfFree3DArray(od);
  //fprintf(stdout,"Freeing wd!\n");
  hdfFree3DArray(wd);
  if (fname_scp_set) {
    //fprintf(stdout,"Freeing scp!\n");
    hdfFree3DArray(scp);
  }
  if (fname_sur_set) {
    //fprintf(stdout,"Freeing sur!\n");
    hdfFree3DArray(sur);
  }

  return 0;
}
