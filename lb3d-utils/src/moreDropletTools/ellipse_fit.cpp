#include "ellipse_fit.h"
#include <cmath>

using namespace std;

/*

Fitting of points to F(x,y) = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0 based on "Numerically stable direct least squares fitting of ellipses"  by R. Halir and J. Flusser
Naming conventions are as in the article.

Calculation of deformation and angles as per http://mathworld.wolfram.com/Ellipse.html
New temporary variables are defined for one-to-one matching of the variable names on that page.

*/

EllipseFit::EllipseFit(const vector< vector<double> > pvec, VerbosityLevel verb) {
  ELLIPSE_FIT_VERBOSITY = verb;
  p = pvec;
}

int EllipseFit::doFit() {
  int vdim   = 1;         // smaller dimension of a column/row [V]ector (1)
  int smdim  = 3;         // dimension of a [S]ub[M]atrix
  int wsmdim = 8*smdim;   // dimension of a [W]ork [S]ub[M]atrix (estimate)
  const char transpose = 'T';
  const char noop = 'N';
  char charV = 'V';
  char charN = 'N';
  int npts = p.size(); // number of datapoints
  double alpha, beta;     // parameters for dgemm calls (C := alpha*(op(A).op(B)) + beta C )
  int ipiv[wsmdim];       // ???
  double work[wsmdim];    // work matrix, temp space
  int info = 0;           // return value of some LAPACK operations

  // ATTENTION, for compatibility with LAPACK we use MATRIX[col][row] --- FORTRAN ORDER!

  if (npts < 1) {
    return 1;
  }

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"N = %d\n",npts);
  }

  if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
    for (int i = 0; i < npts ; i++) {
      fprintf(stdout,"{%f, %f},\n",p[i][0],p[i][1]);
    }
  }

  // Creating D1 from vector of datapoints:
  // Rows of D1 are ( x^2   x*y   y^2 )
  // x[i] = pvec[i][0], y[i] = pvec[i][1]
  double D1[smdim][npts];
  for (int i = 0; i < npts ; i++) {
    D1[0][i] = (p[i][0])*(p[i][0]);
    D1[1][i] = (p[i][0])*(p[i][1]);
    D1[2][i] = (p[i][1])*(p[i][1]);
  }
  // Creating D2 from vector of datapoints:
  // Rows of D2 are ( x     y     1   )
  double D2[smdim][npts];
  for (int i = 0; i < npts ; i++) {
    D2[0][i] = p[i][0];
    D2[1][i] = p[i][1];
    D2[2][i] = 1;
  }

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"D = (D1 | D2):\n");
    for (int i = 0; i < npts ; i++) {
      fprintf(stdout,"\t%f %f %f %f %f %f\n",D1[0][i],D1[1][i],D1[2][i],D2[0][i],D2[1][i],D2[2][i]);
    }
  }

  // Calculate S1 = D1' D1
  alpha = 1.0;
  beta  = 0.0;
  double S1[smdim][smdim];
  dgemm(&transpose, &noop, &smdim, &smdim, &npts, &alpha, D1[0], &npts, D1[0], &npts, &beta, S1[0], &smdim);

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"S1 = D1' D1:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S1[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate S2 = D1' D2
  alpha = 1.0;
  beta  = 0.0;
  double S2[smdim][smdim];
  dgemm(&transpose, &noop, &smdim, &smdim, &npts, &alpha, D1[0], &npts, D2[0], &npts, &beta, S2[0], &smdim);

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"S2 = D1' D2:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S2[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate S2 = D2' D2
  alpha = 1.0;
  beta  = 0.0;
  double S3[smdim][smdim];
  dgemm(&transpose, &noop, &smdim, &smdim, &npts, &alpha, D2[0], &npts, D2[0], &npts, &beta, S3[0], &smdim);

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"S3 = D2' D2:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S3[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // C1 is actually never needed in the calculation, so we immediately construct C1I
  double C1I[smdim][smdim];
  for (int i = 0; i < smdim ; i++ ) {
    for (int j = 0; j < smdim ; j++ ) {
      C1I[j][i] = 0.0;
    }
  }
  C1I[2][0] =  0.5;
  C1I[1][1] = -1.0;
  C1I[0][2] =  0.5;

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"C1I:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",C1I[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate S3I
  // We first have to run dgetrf_ to decompose the matrix in preparation for the dgetri_ call which performs the actual inversion.
  // Also, the inversion operation overwrites the original matrix, so we make a copy of S3 first.
  // See LAPACK documentation for details.
  double S3I[smdim][smdim];
  for (int i = 0; i < smdim ; i++ )
    for (int j = 0; j < smdim ; j++ )
      S3I[j][i] = S3[j][i];

  dgetrf_(&smdim, &smdim, S3I[0], &smdim, ipiv, &info);
  if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
    fprintf(stdout,"DEBUG: dgetrf_ returned info = %d\n",info);
  }
  dgetri_(&smdim, S3I[0], &smdim, ipiv, work, &wsmdim, &info);
  if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
    fprintf(stdout,"DEBUG: dgetri_ returned info = %d\n",info);
  }

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"S3I:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S3I[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate T = - S3I.S2'
  alpha = -1.0;
  beta  =  0.0;
  double T[smdim][smdim];
  dgemm(&noop, &transpose, &smdim, &smdim, &smdim, &alpha, S3I[0], &smdim, S2[0], &smdim, &beta, T[0], &smdim);

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"T = S3I S2':\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",T[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate S1 - S2.S3I.S2' = S1 + S2.T
  // Using dgemm to calculate this all at one means that we overwrite the matrix S1 with this new value.
  alpha = 1.0;
  beta  = 1.0;
  if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
    fprintf(stdout,"DEBUG: Trashing matrix S1...\n");
  }
  dgemm(&noop, &noop, &smdim, &smdim, &smdim, &alpha, S2[0], &smdim, T[0], &smdim, &beta, S1[0], &smdim);
  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"S1 - S2 S3I S2':\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S1[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Calculate M = C1I.( S1 + S2.S3I.S2' ) = C1I.( S1 + S2.T )
  alpha = 1.0;
  beta  = 0.0;
  double M[smdim][smdim];
  dgemm(&noop, &noop, &smdim, &smdim, &smdim, &alpha, C1I[0], &smdim, S1[0], &smdim, &beta, M[0], &smdim);
  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"M = C1I ( S1 - S2 S3I S2' ):\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",S1[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Now we get to solve the eigenvalue problem M.a1 = lamba*a1
  // lamba is split up in real and imaginary component Re and Im respectively
  // These will contain the eigenvalues matching the eigenvector at the same index
  // VL will contain nothing, since we are not interested in the left eigenvalues
  // VR will contain our (right) eigenvalues
  double Re[smdim];
  double Im[smdim];
  double VL[smdim][smdim];
  double VR[smdim][smdim];
  dgeev_(&charN, &charV, &smdim, M[0], &smdim, Re, Im, VL[0], &smdim, VR[0], &smdim, work, &wsmdim,&info);
  if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
    fprintf(stdout,"DEBUG: dgeev_ returned info = %d\n",info);
  }

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"Eigenvalues (Re,Im):\n");
    for (int i = 0; i < smdim; i++) {
      fprintf(stdout,"\t(%f,%f)\n",Re[i],Im[i]);
    }
    fprintf(stdout,"Right eigenvectors:\n");
    for (int j = 0; j < smdim ; j++) {
      fprintf(stdout,"\t");
      for (int i = 0; i < smdim ; i++) {
        fprintf(stdout,"%f ",VR[j][i]);
      }
      fprintf(stdout,"\n");
    }
  }

  // Almost there... we need to find the smallest positive eigenvalue and take the matching right eigenvector as our answer.
  // aCa = a1' C1 a1
  double aCa[smdim];
  for (int i = 0; i < smdim; i++) {
    aCa[i] = 4*VR[i][0]*VR[i][2]-VR[i][1]*VR[i][1]; // 4ac-b^2, directly calculating a1' C1 a1
  }

  if (ELLIPSE_FIT_VERBOSITY >= REPORT) {
    fprintf(stdout,"a1' C1 a1:\n");
    for (int i = 0; i < smdim; i++) {
      fprintf(stdout,"\t%f\n",aCa[i]);
    }
  }

  int caCai = 0;                // [C]hosen [aCa] [I]ndex
  double caCav = aCa[caCai];    // [C]hosen [aCa] [V]alue
  for (int i = 0; i < smdim; i++) {
    if ((aCa[i] > 0.0) && (aCa[i] < caCav || caCav < 0.0)) {
      caCav = aCa[i];
      caCai = i;
      if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
        fprintf(stdout,"DEBUG: New choice of eigenvector: %d (%lf)\n",i,aCa[i]);
      }
    }
  }

  if (caCav < 0.0) {
    if (ELLIPSE_FIT_VERBOSITY >= DEBUG) {
      fprintf(stdout,"DEBUG: Less than zero, error!");
    }
  }

  // We now take a1 := our chosen eigenvector
  double a1[1][smdim];
  a1[0][0] = VR[caCai][0];
  a1[0][1] = VR[caCai][1];
  a1[0][2] = VR[caCai][2];

  // One last call to dgemm to calculate a2 = T.a1
  alpha = 1.0;
  beta  = 0.0;
  double a2[1][smdim];
  dgemm(&noop, &noop, &smdim, &vdim, &smdim, &alpha, T[0], &smdim, a1[0], &smdim, &beta, a2[0], &smdim);

  // And finally we construct a = ( a1' | a2' ) = ( a b c d e f )'
  a = a1[0][0];
  b = a1[0][1];
  c = a1[0][2];
  d = a2[0][0];
  e = a2[0][1];
  f = a2[0][2];

  /// \todo Should add some sanity checks here, to see if it's an actual ellipse?

  // We have gone through the fitting algorithm now, but we also want to extract some other parameters from 
  // the 6 fitted values...

  // Warning: Wolfram uses F(x,y) = a x^2 + 2 b x y + c y^2 + 2 d x + 2 f y + g = 0
  // For clarity, we define these temporary variables in this style
  double wa = 1.0*a;
  double wb = 0.5*b;
  double wc = 1.0*c;
  double wd = 0.5*d;
  double wf = 0.5*e;
  double wg = 1.0*f;

  // Center position of the ellipse
  x0 = (wc*wd - wb*wf) / (wb*wb - wa*wc);
  y0 = (wa*wf - wb*wd) / (wb*wb - wa*wc);

  // Semi-minor and semi-major axes
  L = sqrt(
        2.0* ( wa*wf*wf + wc*wd*wd + wg*wb*wb - 2.0*wb*wd*wf - wa*wc*wg) /
        ( (wb*wb - wa*wc) * (  1.0 * sqrt( (wa-wc)*(wa-wc) + 4.0*wb*wb ) - (wa+wc) ) )
      );
  B = sqrt(
        2.0* ( wa*wf*wf + wc*wd*wd + wg*wb*wb - 2.0*wb*wd*wf - wa*wc*wg) /
        ( (wb*wb - wa*wc) * ( -1.0 * sqrt( (wa-wc)*(wa-wc) + 4.0*wb*wb) - (wa+wc) ) )
      );

  phi = 0.0;

  // Swap axes and rotate ellipse if required.
  if ( B > L ) {
    double Ltmp = L;
    L = B;
    B = Ltmp;
    phi = 0.5*M_PI;
  }
  else {
    phi = 0.0;
  }

  // Calculate deformation parameter
  D = ( L - B ) / ( L + B );

  // Calculate rotation angle
  if ( abs(b) < EPSILON ) {
    if ( wa < wc ) {
      phi += 0.0;
    }
    else {
      phi += 0.5*M_PI;
    }
  }
  else {
    if ( wa < wc ) {
      phi += 0.5*atan( 2.0*wb / (wa - wc) ); // arccot(x) = arctan(1/x)
    }
    else {
      phi += 0.5*M_PI + 0.5*atan( 2.0*wb / (wa - wc) ); // arccot(x) = arctan(1/x)
    }
  }

  if (ELLIPSE_FIT_VERBOSITY >= RESULT) {
    displayResult();
  }

  return 0;
}

void EllipseFit::writeMathematicaFormat(char *filename, const int xmax, const int ymax) {
  // This will write a Mathematica-compatible .txt file which can be evaulated to check the fit.
  // We need to supply a plot range, so xmax and ymax are mandatory arguments.
  // These are normally the "other" elements of the system size, i.e. not the slice direction.
  int npts = p.size();
  FILE *fp;
  fp = fopen (filename,"w");
  fprintf(fp,"a = %f;\n",a);
  fprintf(fp,"b = %f;\n",b);
  fprintf(fp,"c = %f;\n",c);
  fprintf(fp,"d = %f;\n",d);
  fprintf(fp,"e = %f;\n",e);
  fprintf(fp,"f = %f;\n",f);
  fprintf(fp,"l = {");
  for (int i = 0; i < npts ; i++) {
    if (i == npts - 1 ) {
      fprintf(fp,"{%f, %f, 0}" ,p[i][0],p[i][1]);
    }
    else {
      fprintf(fp,"{%f, %f, 0},",p[i][0],p[i][1]);
    }
  }
  fprintf(fp,"};\n");
  fprintf(fp,"Plot3D[a*x^2 + b*x*y + c*y^2 + d*x + e*y + f, {x, 0, %d}, {y, 0, %d}];\n", xmax, ymax);
  fprintf(fp,"Plot3D[0, {x, 0, %d}, {y, 0, %d}];\n", xmax, ymax);
  fprintf(fp,"ListPointPlot3D[l];\n");
  fprintf(fp,"Show[%%, %%%%, %%%%%%, PlotRange -> Automatic]\n");
  fclose(fp);
}


void EllipseFit::displayResult() {
  // This just dumps the final result to stdout, in human-readable format.
  int npts = p.size(); // number of datapoints
  fprintf(stdout,"\n==========================================================================\n");
  fprintf(stdout,"\nFitted ellipse through %d points.\n",npts);
  fprintf(stdout,"\nEllipse parameters of F(x,y) = a x^2 + b x y + c y^2 + d x + e y + f = 0:\n \ta = %f\n \tb = %f\n \tc = %f\n \td = %f\n \te = %f\n \tf = %f\n",a,b,c,d,e,f);
  fprintf(stdout,"\nCenter coordinates of ellipse (x0, y0):\n \t(%f, %f)\n",x0,y0);
  fprintf(stdout,"\nLength and breadth of ellipse (L, B):\n \t(%f, %f)\n",L,B);
  fprintf(stdout,"\nDeformation parameter D:\n \t%f\n",D);
  fprintf(stdout,"\nCounterclockwise angle of rotation from the x-axis\n  to the major axis of the ellipse phi:\n \t%f (%f degrees)\n",phi,180.0*phi/M_PI);
  fprintf(stdout,"\n==========================================================================\n"); 
}

/*

Accessor functions, yay!

*/

double EllipseFit::getDeformation() {
  return D;
}

double EllipseFit::getRotation() {
  return phi;
}

double EllipseFit::getLength() {
  return L;
}

double EllipseFit::getBreadth() {
  return B;
}

vector<double> EllipseFit::getCentre() {
  vector<double> centre;
  centre.push_back(x0);
  centre.push_back(y0);
  return centre;
}

vector<double> EllipseFit::getFunctionParameters() {
  vector<double> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  params.push_back(d);
  params.push_back(e);
  params.push_back(f);
  return params;
}

