#ifndef NOMKL
#include "droplet_fit.h"
#else
#include "hdf5_helper.h"
#endif

#include <string.h>
#include <math.h>

using namespace std;

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz;

  double x, y, z;
  double dx, dy ,dz;

  // File names
  char *fname_colour, *fname_od;

  // Array pointers
  float ***colour, ***od;

  double D_F = 0.0;

  double M, m, cutoff;
  // Tensor
  // double Axx, Axy, Axz, Ayy, Ayz, Azz;
  double Rx, Ry, Rz;
  double Ixx, Ixy, Ixz, Iyy, Iyz, Izz;

  // Parameters
  double radius  = 0.0;
  double rho     = 0.30;
  double shear_u = 0.050;
  double sigma   = 0.0;
  double tau     = 1.0;

  double Dradius = 0.0;
  double Drho    = 0.0;
  double Dsigma  = 0.0;

  int t = 0;

  hsize_t *dims;

  bool fname_colour_set = false;
  bool fname_od_set = false;
  bool header = false;

  HDF_VERBOSITY = NONE;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--colourfile") == 0 ) {
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
    else if ( strcmp(argv[i],"-cutoff") == 0 ) {
      if ( i+1 < argc ) {
	cutoff = atof(argv[++i]);
	//fprintf(stdout,"  Cutoff <%f>\n",cutoff);
      }
      else {
	fprintf(stderr,"  Missing argument to -cutoff flag. \n");
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
    else if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--header") == 0    ) {
      header = true;
      //fprintf(stdout,"  Header\n");
    }
    else if ( strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--odfile") == 0 ) {
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
    else if ( strcmp(argv[i],"-shear_u") == 0 ) {
      if ( i+1 < argc ) {
	shear_u = atof(argv[++i]);
	//fprintf(stdout,"  Shear_u <%f>\n",shear_u);
      }
      else {
	fprintf(stderr,"  Missing argument to -shear_u flag. \n");
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
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }

  }

  if (! ( fname_colour_set && fname_od_set ) ) {
    fprintf(stderr,"Colour file -c and od file -o are mandatory.\n\n");
    fprintf(stderr,"  -c --colourfile  <colourfile>\n");
    fprintf(stderr,"  -cutoff          <density cutoff>\n");
    fprintf(stderr,"  -o --odfile      <odfile>\n");
    fprintf(stderr,"  -radius          <initial droplet radius>\n");
    fprintf(stderr,"  -rho             <density>\n");
    fprintf(stderr,"  -shear_u         <shear_u>\n");
    fprintf(stderr,"  -sigma           <surface tension>\n");
    fprintf(stderr,"  -t --timestep    <timestep>\n");
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
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  //  LB3dAttributes attr(fname_colour);
  
//   int niter;
//   double g_br;
//   string restore_string;
//   bool MRT;
//   int st;
//   st = attr.getInteger("n_iteration",&niter);
//   if (st == 0) {
//     fprintf(stdout,"\nn_iteration = %d\n",niter);
//   }
//   st = attr.getDouble("g_br",&g_br);
//   if (st == 0) {
//     fprintf(stdout,"\ng_br = %lf\n",g_br);
//   }
//   st = attr.getString("restore_string",&restore_string);
//   if (st == 0) {
//     fprintf(stdout,"\nrestore_string = %s\n",restore_string.c_str());
//   }
//   st = attr.getBool("MRT",&MRT);
//   if (st == 0) {
//     fprintf(stdout,"\nMRT = %d\n",MRT);
//   }

  //fprintf(stdout,"\nAnalyzing file %s...\n",file_in);


  // double R;
  // double fr1;
  // double shear_u;
  // double gdot;
  // double tau;
  // double fr, fb;
  // double sigmaCa;
  
  // attr.getDouble("fr1",&fr1);
  // attr.getInteger("nx",&nx);
  // R = fr1 * nx;
  // if ( attr.getDouble("shear_u",&shear_u)==LB3dAttributes::error) {
  //   exit(1);
  // }
  // gdot = shear_u / nx;
  // attr.getDouble("tau_r",&tau);
  // attr.getDouble("fr",&fr);
  // attr.getDouble("fb",&fb);
  
  // fprintf(stdout,"fr1: %lf, nx: %d, R: %lf, shear_u: %lf, gdot: %lf\n",fr1,nx,R,shear_u,gdot);
  
  // sigmaCa = ( (2.0*tau-1.0)*(fr+fb)*gdot*R ) / 6;

  // fprintf(stdout,"sigmaCa: %lf\n",sigmaCa);

  // Calculate centre of mass

  Rx = 0.0; Ry = 0.0; Rz = 0.0;
  M = 0.0;
  for (int k=0; k < nz; k++) {
    for (int j=0; j < ny; j++) {
      for (int i=0; i < nx; i++) {
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

  // dx = (nx - 1.0)/2.0;
  // dy = (ny - 1.0)/2.0;
  // dz = (nz - 1.0)/2.0;

  // Use CoM as offset

  dx = Rx;
  dy = Ry;
  dz = Rz;

  // for ( double thetad = 0.0; thetad < 20.0 ; thetad += 1.0 ) {
  // double theta = thetad * 2 * M_PI / 360.0;

  // Axx = 0.0; Axy = 0.0; Axz = 0.0; Ayy = 0.0; Ayz = 0.0; Azz = 0.0;
  Ixx = 0.0; Ixy = 0.0; Ixz = 0.0; Iyy = 0.0; Iyz = 0.0; Izz = 0.0;
  for (int k=0; k < nz; k++) {
    for (int j=0; j < ny; j++) {
      for (int i=0; i < nx; i++) {
	x = i - dx;
	y = j - dy;
	z = k - dz;

	if ( x > nx/2) x -= nx;
	if ( y > ny/2) y -= ny;
	if ( z > nz/2) z -= nz;

	if ( x < -nx/2) x += nx;
	if ( y < -ny/2) y += ny;
	if ( z < -nz/2) z += nz;

	// x = x;
	// y = y*cos(theta) - z*sin(theta);
	// z = y*sin(theta) + z*cos(theta);

	// x = x*cos(theta) + z*sin(theta);
	// y = y;
	// z = -x*sin(theta) + z*cos(theta);

	// x = x*cos(theta) - y*sin(theta);
	// y = x*sin(theta) + y*cos(theta);
	// z = z;
			  
	m = od[k][j][i];
	if ( m > cutoff) {
	  // Axx += x*x*m;
	  // Axy += x*y*m;
	  // Axz += x*z*m;
	  // Ayy += y*y*m;
	  // Ayz += y*z*m;
	  // Azz += z*z*m;
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

  // Axx /= M;
  // Axy /= M;
  // Axz /= M;
  // Ayy /= M;
  // Ayz /= M;
  // Azz /= M;

  // Ixx = Axx;
  // Ixy = Axy;
  // Ixz = Axz;
  // Iyy = Ayy;
  // Iyz = Ayz;
  // Izz = Azz;

  // Solve eigenvalues for a 3x3 symmetric matrix
  double Mm = (Ixx + Iyy + Izz)/3.0;

  double Kxx = Ixx - Mm;
  double Kxy = Ixy;
  double Kxz = Ixz;
  double Kyy = Iyy - Mm;
  double Kyz = Iyz;
  double Kzz = Izz - Mm;
  double Kdet = Kxx*(Kyy*Kzz - Kyz*Kyz) - Kxy*(Kxy*Kzz - Kyz*Kxz) + Kxz*(Kxy*Kyz - Kyy*Kxz);

  double p = (Kxx*Kxx + Kyy*Kyy + Kzz*Kzz + 2*Kxy*Kxy + 2*Kxz*Kxz + 2*Kyz*Kyz)/6.0;

  double phi = (1.0/3.0)*acos( (Kdet/2.0) / (pow(p,1.5)));

  // Some edge cases
  if (fabs(Kdet/2.0) >= fabs(pow(p,1.5))) phi = 0;
  if (phi < 0) phi += M_PI/3.0;

  // Yay, eigenvalues
  double e1 = Mm + 2*sqrt(p)*cos(phi);
  double e2 = Mm - sqrt(p)*(cos(phi) + sqrt(3.0)*sin(phi));
  double e3 = Mm - sqrt(p)*(cos(phi) - sqrt(3.0)*sin(phi));

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

  // double a = sqrt( (1.0/2.0) * (E2 + E3 - E1) );
  // double b = sqrt( (1.0/2.0) * (E3 + E1 - E2) );
  // double c = sqrt( (1.0/2.0) * (E1 + E2 - E3) );

  double L_I = max(a,max(b,c));
  double B_I = min(a,min(b,c));
  // double M_I = ( a > B_I && a < L_I ) ? a : ( ( b > B_I && b < L_I ) ? b : c );
  double D_I = (L_I - B_I)/(L_I + B_I);

  //fprintf(stdout,"%f ; %f %f %f ; %f %f %f ; %f %f \n",thetad,a,b,c,L_I,M_I,B_I,D_I,M_I/B_I);
  //}

  // Capillary number?

  //rho = 1.0 - pow(M_E,-rho);
  double Ca = ( ( 2.0*tau - 1.0 ) * rho * shear_u ) / ( 3.0 * nx ) ;
  Ca *= radius/sigma;

  double CaP = ( ( 2.0*tau - 1.0 ) * (rho+Drho) * shear_u ) / ( 3.0 * nx ) ;
  CaP *= (radius+Dradius)/(sigma-Dsigma);

  double CaM = ( ( 2.0*tau - 1.0 ) * (rho-Drho) * shear_u ) / ( 3.0 * nx ) ;
  CaM *= (radius-Dradius)/(sigma+Dsigma);

  double DCa = max(CaP-Ca,Ca-CaM);

  double Re = ( radius * radius * 12.0 * shear_u ) / ( ( 2.0*tau - 1.0) * nx ) ;

#ifdef DBGMSG
  //fprintf(stdout,"%f %f %f %f %f %f %f\n",M,Axx,Axy,Axz,Ayy,Ayz,Azz);
  fprintf(stdout,"R = %f , %f , %f\n",Rx, Ry,Rz);
  fprintf(stdout,"{%f,%f,%f,%f,%f,%f,%f,%f}\n",M,Ixx,Ixy,Ixz,Iyy,Iyz,Izz,Kdet);
  fprintf(stdout,"{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}\n",Ixx,Ixy,Ixz,Ixy,Iyy,Iyz,Ixz,Iyz,Izz);
  fprintf(stdout,"EV: %f, %f, %f\n",e1,e2,e3);
  fprintf(stdout,"Axes: %f, %f, %f\n",a,b,c);
#endif


#ifndef NOMKL

  //fprintf(stdout,"\nAnalyzing X-slice:\n");
  //EllipseFit myfitX = getEllipseFitFromSlice(fname_colour,X);
  //myfitX.displayResult();
  //myfitX.writeMathematicaFormat("myfitX.txt", ny, nz);

  //fprintf(stdout,"\nAnalyzing Y-slice:\n");
  EllipseFit myfitY = getEllipseFitFromSlice(fname_colour,Y);
  //myfitY.displayResult();
  //myfitY.writeMathematicaFormat("myfitY.txt", nz, nx);

  //fprintf(stdout,"\nAnalyzing Z-slice:\n");
  //EllipseFit myfitZ = getEllipseFitFromSlice(fname_colour,Z);
  //myfitZ.displayResult();
  //myfitZ.writeMathematicaFormat("myfitZ.txt", nx, ny);

  // fprintf(stdout,"# %d %f %f %f %f %f %f %f %f %f\n", t,
  //   myfitX.getLength(), myfitX.getBreadth(), myfitX.getRotation(),
  //   myfitY.getLength(), myfitY.getBreadth(), myfitY.getRotation(),
  //   myfitZ.getLength(), myfitZ.getBreadth(), myfitZ.getRotation());

  D_F = (myfitY.getLength()-myfitY.getBreadth()) / (myfitY.getLength()+myfitY.getBreadth());
  
#endif

//  fprintf(stdout,"D_F = %f , D_I = %f, Ca = %f, 35/32 Ca = %f\n", D_F, D_I, Ca, 35.0*Ca/32.0);

    if (header) {
      fprintf(stdout,"#\n");
      fprintf(stdout,"# GENERATED BY: droplet-deformation\n");
      fprintf(stdout,"#\n");
      fprintf(stdout,"# R       = %E , Dradius = %E\n", radius, Dradius);
      fprintf(stdout,"# sigma   = %E , Dsigma  = %E\n", sigma, Dsigma);
      fprintf(stdout,"# shear_u = %f\n", shear_u);
      fprintf(stdout,"# rho     = %f     , Drho    = %E\n", rho, Drho);
      fprintf(stdout,"# cutoff  = %f\n", cutoff);
      fprintf(stdout,"#\n");
      fprintf(stdout,"# CoM     = (%f , %f , %f) , M = %f\n",Rx, Ry, Rz, M);
      fprintf(stdout,"# I : %f , %f , %f , %f , %f , %f\n",Ixx, Ixy, Ixz, Iyy, Iyz, Izz);
      fprintf(stdout,"# EV: %f , %f , %f\n", e1, e2, e3);
      fprintf(stdout,"# I Axes  :  %f , %f , %f\n", a, b, c);
      fprintf(stdout,"#      t          D_F          D_I           Ca          DCa           Re\n");
      fprintf(stdout,"# =======================================================================\n");
    }

    fprintf(stdout,"%#8.8d %E %E %E %E %E\n", t, D_F, D_I, Ca, DCa, Re);

  return 0;
}
