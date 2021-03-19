#include <string.h>
#include <math.h>
//#include <limits>

#include "droplet-deformation.h"
#include "hdf5_helper.h"
#include "lbe_version.h"
#include "dsyevj3.h"

using namespace std;

double uncorr_error_quot(const double a, const double s_a, const double b, const double s_b) {
  double f2 = (a/b)*(a/b);
  double s_f2 =  f2 * ( (s_a/a) * (s_a/a) + (s_b/b) * (s_b/b) );
  return sqrt(s_f2);
}

double uncorr_error_product(const double a, const double s_a, const double b, const double s_b) {
  double f2 = (a*b)*(a*b);
  double s_f2 =  f2 * ( (s_a/a) * (s_a/a) + (s_b/b) * (s_b/b) );
  return sqrt(s_f2);
}

int find_extentminz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int k=0;k<dz;k++) {
    for(int j=0;j<dy;j++) {
      for(int i=0;i<dx;i++) {
	if (N[k][j][i] > cutoff) {
	  return k;
	}
      }
    }
  }
  return -1;
}

int find_extentmaxz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int k=dz-1;k>=0;k--) {
    for(int j=0;j<dy;j++) {
      for(int i=0;i<dx;i++) {
	if (N[k][j][i] > cutoff) {
	  return k;
	}
      }
    }
  }
  return -1;
}

double gammadot(float ***vel, const int dx, const int j, const int k) {
  double gd;
  for(int i=1;i<dx;i++) {
    if (vel[k][j][i] < vel[k][j][i-1]) {
      gd = (vel[k][j][i] - vel[k][j][0]) / (double) i;
      return gd;
    }
  }
  gd = (vel[k][j][dx-1]-vel[k][j][0]) / (double) (dx -1);
  return gd;
  //return numeric_limits<double>::quiet_NaN();;
}

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz;

  double dx, dy ,dz;

  // File names
  char *fname_od, *fname_vel;

  // Array pointers
  float ***od, ***vel;

  double M, cutoff;
  // Tensor
  // double Axx, Axy, Axz, Ayy, Ayz, Azz;
  double Rx, Ry, Rz;
  double Ixx, Ixy, Ixz, Iyy, Iyz, Izz;

  // Parameters
  double radius  = 0.0;
  double rho     = 0.0;
  double shear_u = 0.0;
  double sigma   = 0.0;
  double tau     = 0.0;

  double Dradius = 0.0;
  double Drho    = 0.0;
  double Dsigma  = 0.0;

  int t = 0;

  hsize_t *dims;

  bool fname_od_set = false;
  bool fname_vel_set = false;
  //bool dir_set = false;
  bool header = false;

  HDF_VERBOSITY = NONE;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-c") == 0 ) {
      if ( i+1 < argc ) {
	cutoff = atof(argv[++i]);
	//fprintf(stdout,"  Cutoff <%f>\n",cutoff);
      }
      else {
	fprintf(stderr,"  Missing argument to -c flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-Dradius") == 0 ) {
      if ( i+1 < argc ) {
	Dradius = atof(argv[++i]);
	//fprintf(stdout,"  Dradius <%f>\n",Dradius);
      }
      else {
	fprintf(stderr,"  Missing argument to -Dradius flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-Drho") == 0 ) {
      if ( i+1 < argc ) {
	Drho = atof(argv[++i]);
	//fprintf(stdout,"  Drho <%f>\n",Drho);
      }
      else {
	fprintf(stderr,"  Missing argument to -Drho flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-Dsigma") == 0 ) {
      if ( i+1 < argc ) {
	Dsigma = atof(argv[++i]);
	//fprintf(stdout,"  Dsigma <%f>\n",Dsigma);
      }
      else {
	fprintf(stderr,"  Missing argument to -Dsigma flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-h") == 0 ) {
      header = true;
      //fprintf(stdout,"  Header\n");
    }
    else if ( strcmp(argv[i],"-o") == 0 ) {
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
    else if ( strcmp(argv[i],"-radius") == 0 ) {
      if ( i+1 < argc ) {
	radius = atof(argv[++i]);
	//fprintf(stdout,"  Radius <%f>\n",radius);
      }
      else {
	fprintf(stderr,"  Missing argument to -radius flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-rho") == 0 ) {
      if ( i+1 < argc ) {
	rho = atof(argv[++i]);
	//fprintf(stdout,"  Rho <%f>\n",rho);
      }
      else {
	fprintf(stderr,"  Missing argument to -rho flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-sigma") == 0 ) {
      if ( i+1 < argc ) {
	sigma = atof(argv[++i]);
	//fprintf(stdout,"  Sigma <%f>\n",sigma);
      }
      else {
	fprintf(stderr,"  Missing argument to -sigma flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--timestep") == 0 ) {
      if ( i+1 < argc ) {
	t = atoi(argv[++i]);
	//fprintf(stdout,"  Time step <%d>\n",t);
      }
      else {
	fprintf(stderr,"  Missing argument to -t flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-v") == 0 ) {
      if ( i+1 < argc ) {
	fname_vel = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname_vel);
	fname_vel_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -v flag. \n");
	exit(1);
      }
    }
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }

  }

  if (! ( fname_od_set ) ) {
    fprintf(stdout,"droplet-deformation (%s)\n\n", GIT_DESC);
    fprintf(stderr,"od file -o is mandatory.\n\n");
    fprintf(stderr,"  -c               <density cutoff>\n");
    fprintf(stderr,"  -h               \n");
    fprintf(stderr,"  -o               <odfile>\n");
    fprintf(stderr,"  -radius          <initial droplet radius>\n");
    fprintf(stderr,"  -Dradius         <error in initial droplet radius>\n");
    fprintf(stderr,"  -rho             <density>\n");
    fprintf(stderr,"  -Drho            <error in density>\n");
    fprintf(stderr,"  -sigma           <surface tension>\n");
    fprintf(stderr,"  -Dsigma          <error in surface tension>\n");
    fprintf(stderr,"  -t               <timestep>\n");
    fprintf(stderr,"  -v               <velzfile>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  if (fname_vel_set) {
    if (hdfReadAll(fname_vel,&vel,&dims) != 0) {
      exit(1);
    }
    hdfFreeDims(dims);
  }
  if (hdfReadAll(fname_od    ,&od    ,&dims) != 0) {
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  int np_sphere = 0;
  int np_sphere_cap = 0;
  int np = 0;
  double Rp = 0.0;
  double Ro = 0.0;
  // Get attributes
  LB3dAttributes attr = LB3dAttributes(fname_od);
  attr.getDouble("shear_u",&shear_u);
  attr.getDouble("tau_r",&tau);
  attr.getInteger("np_sphere",&np_sphere);
  attr.getInteger("np_sphere_cap",&np_sphere_cap);
  attr.getDouble("R_para",&Rp);
  attr.getDouble("R_orth",&Ro);
  np = np_sphere - 2*np_sphere_cap;

  double gd_imp = 2.0 * shear_u / (nx - 1.0);

  double gd_max = 0.0;
  double gd_avg_drop = 0.0;
  double gd_avg_all = 0.0;
  
  if (fname_vel_set) {
    // Do effective shear calculations

    int extentmin = find_extentminz(od,nx,ny,nz,cutoff);
    int extentmax = find_extentmaxz(od,nx,ny,nz,cutoff);

    // fprintf(stdout,"# Drop extent: z-min = %d , z-max = %d\n",extentmin,extentmax);

    // Found extent

    // if (header) {
    //   fprintf(stdout,"# z    calc     imp\n");
    // }

    int j = (ny/2)-1;

    for (int k=0;k<nz;k++) {
      double gd = gammadot(vel,nx,j,k);
      if (gd > gd_max) gd_max = gd;
      if (k <= extentmax && k >= extentmin) gd_avg_drop += gd;
      gd_avg_all += gd;

      // fprintf(stdout,"%d %lf %lf\n",k,gd,2.0*shear_u/(nx-1.0));
    }

    gd_avg_drop /= (double) (extentmax-extentmin+1);
    gd_avg_all /= nz;

    // fprintf(stdout,"# gd_max = %lf, gd_avg_drop = %lf, gd_avg_all = %lf\n",gd_max,gd_avg_drop,gd_avg_all);
  }

  // Calculate centre of mass

  Rx = 0.0; Ry = 0.0; Rz = 0.0;
  M = 0.0;
  for (int k=0; k < nz; k++) {
    for (int j=0; j < ny; j++) {
      for (int i=0; i < nx; i++) {
	double m;
	m = od[k][j][i];
	if ( m > cutoff) {
	  M  += m;
	  Rx += m*i;
	  Ry += m*j;
	  Rz += m*k;
	}
      }
    }
  }

  Rx /= M;
  Ry /= M;
  Rz /= M;

  // Use CoM as offset

  dx = Rx;
  dy = Ry;
  dz = Rz;

  // Axx = 0.0; Axy = 0.0; Axz = 0.0; Ayy = 0.0; Ayz = 0.0; Azz = 0.0;
  Ixx = 0.0; Ixy = 0.0; Ixz = 0.0; Iyy = 0.0; Iyz = 0.0; Izz = 0.0;
  for (int k=0; k < nz; k++) {
    for (int j=0; j < ny; j++) {
      for (int i=0; i < nx; i++) {
	double x = i - dx;
	double y = j - dy;
	double z = k - dz;

	if ( x > nx/2) x -= nx;
	if ( y > ny/2) y -= ny;
	if ( z > nz/2) z -= nz;

	if ( x < -nx/2) x += nx;
	if ( y < -ny/2) y += ny;
	if ( z < -nz/2) z += nz;

	double m = od[k][j][i];
	if ( m > cutoff) {
	  Ixx += m*(y*y+z*z);
	  Ixy -= m*x*y;
	  Ixz -= m*x*z;
	  Iyy += m*(x*x+z*z);
	  Iyz -= m*y*z;
	  Izz += m*(x*x+y*y);
	}
      }
    }
  }

  double I[3][3]; // Moment of inertia tensor (symmetric)
  double evec[3][3]; // Room for eigenvectors
  double eval[3]; // Room for eigenvalues
 
  I[0][0] = Ixx;
  I[0][1] = Ixy;
  I[0][2] = Ixz;
  I[1][1] = Iyy;
  I[1][2] = Iyz;
  I[2][2] = Izz;

  // These should not be needed for the calculation, actually - cf. dsyevj3.cpp
  I[1][0] = Ixy;
  I[2][0] = Ixz;
  I[2][1] = Iyz;

 // Calculate eigensystem
  if ( dsyevj3(I, evec, eval) != 0) {
    fprintf(stderr,"Solving of eigensystem failed.\n");
  }

  // tan(theta) = x / z
  // Use correct eigenvector (smallest eigenvalue)
  int evpos = 0;
  double ev = eval[evpos];
  double theta = atan(evec[0][evpos] / evec[2][evpos]);

  
  if ( eval[1] < ev ) {
    evpos = 1;
    ev = eval[evpos];
    theta = atan(evec[0][evpos] / evec[2][evpos]);

  }

  if ( eval[2] < ev ) {
    evpos = 2;
    ev = eval[evpos];
    theta = atan(evec[0][evpos] / evec[2][evpos]);
  }

  double e1 = eval[0];
  double e2 = eval[1];
  double e3 = eval[2];

  // Now get the axes out
  double E1 = 5.0*e1/M;
  double E2 = 5.0*e2/M;
  double E3 = 5.0*e3/M;

  double a2 = 0.5*(E2 - E1 + E3);
  double b2 = E3 - a2;
  double c2 = E2 - a2;

  double a = sqrt(a2);
  double b = sqrt(b2);
  double c = sqrt(c2);

  double L = max(a,max(b,c));
  double B  = min(a,min(b,c));
  double D = (L - B)/(L + B);

  // Capillary number?
  double mu = ( 2.0*tau - 1.0 ) * rho / 6.0;
  double Dmu = ( 2.0*tau - 1.0 ) * Drho / 6.0;

  double radius_sigma = radius / sigma;
  double Dradius_sigma = uncorr_error_quot(radius, Dradius, sigma, Dsigma);

  double mu_radius_sigma = mu * radius_sigma;
  double Dmu_radius_sigma = uncorr_error_product(mu, Dmu, radius_sigma, Dradius_sigma);

  double Ca = gd_imp * mu_radius_sigma;
  double DCa = gd_imp * Dmu_radius_sigma;

  double Ca_gdavg = gd_avg_drop * mu_radius_sigma;
  double DCa_gdavg = gd_avg_drop * Dmu_radius_sigma;

  double Ca_gdmax = gd_max * mu_radius_sigma;
  double DCa_gdmax = gd_max * Dmu_radius_sigma;

  double radius2 = radius * radius;
  double Dradius2 = 4*Dradius; /// ???

  double Re = radius2 * ( 6.0 * gd_imp /  ( 2.0*tau - 1.0) ) ;
  double DRe = Dradius2 * ( 6.0 * gd_imp /  ( 2.0*tau - 1.0) ) ;

  double chi = np * ( ( radius - sqrt( radius * radius - Rp * Ro ) ) / (2.0* radius) );

  if (header) {
    fprintf(stdout,"#\n");
    fprintf(stdout,"# GENERATED BY: droplet-deformation (%s)\n", GIT_DESC);
    fprintf(stdout,"#\n");
    fprintf(stdout,"# radius  = %E , Dradius = %E\n", radius, Dradius);
    fprintf(stdout,"# rho     = %E , Drho    = %E\n", rho, Drho);
    fprintf(stdout,"# sigma   = %E , Dsigma  = %E\n", sigma, Dsigma);
    fprintf(stdout,"# shear_u = %f\n", shear_u);
    fprintf(stdout,"# tau     = %f\n", tau);
    fprintf(stdout,"# cutoff  = %f\n", cutoff);
    fprintf(stdout,"#\n");
    fprintf(stdout,"# CoM     = (%f , %f , %f) , M = %f\n",Rx, Ry, Rz, M);
    fprintf(stdout,"# I : { {%f , %f , %f} , {%f , %f , %f} , {%f,%f,%f}}\n",Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz);
    fprintf(stdout,"# Eigenvalues %f %f %f\n",eval[0],eval[1],eval[2]);
    fprintf(stdout,"# Using eigenvector #%d\n",evpos+1);
    fprintf(stdout,"# Eigenvector #1 %f %f %f\n",evec[0][0],evec[1][0],evec[2][0]);
    fprintf(stdout,"# Eigenvector #2 %f %f %f\n",evec[0][1],evec[1][1],evec[2][1]);
    fprintf(stdout,"# Eigenvector #3 %f %f %f\n",evec[0][2],evec[1][2],evec[2][2]);
    fprintf(stdout,"# I Axes  :  %f , %f , %f\n", a, b, c);
    fprintf(stdout,"#? %5.5s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s %12.12s\n", "t","D","theta","Ca","DCa","Ca_gdavg","DCa_gdavg","Ca_gdmax","DCa_gdmax","Re","DRe","gd_imp","gd_avg_drop","gd_max","chi");
    fprintf(stdout,"# ============================================================================================================================================================================================\n");
  }

  fprintf(stdout,"%#8.8d %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E %12.6E\n", t, D, theta, Ca, DCa, Ca_gdavg, DCa_gdavg, Ca_gdmax, DCa_gdmax, Re, DRe, gd_imp, gd_avg_drop, gd_max, chi);

  //fprintf(stdout,"Freeing od!\n");
  hdfFree3DArray(od);
  if (fname_vel_set) {
    //fprintf(stdout,"Freeing velz!\n");
    hdfFree3DArray(vel);
  }


  return 0;
}
