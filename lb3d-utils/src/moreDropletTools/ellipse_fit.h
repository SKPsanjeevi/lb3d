#include <stdio.h>
#include <vector>
#include <iostream>

#ifndef VERBOSITY_LEVEL_H
#define VERBOSITY_LEVEL_H
#include "verbosity_level.h"
#endif

#ifndef CARDINALDIRECTION_H
#define CARDINALDIRECTION_H
#include "cardinaldirection.h"
#endif

#include "mkl.h"

using namespace std;

// Used for checking if float == 0, sorta
const double EPSILON = 0.00000000000001;

class EllipseFit {
  VerbosityLevel ELLIPSE_FIT_VERBOSITY;
  vector< vector<double> > p;   // Datapoints are stored here
  double a, b, c, d, e, f;      // Fitting parameters
  double x0, y0, L, B, D, phi;  // Extended fitting parameters

  public:
    // Constructor
    EllipseFit(const vector< vector<double> >, VerbosityLevel);

    // Calculation
    int doFit();

    // Accessor functions
    vector<double> getFunctionParameters();
    vector<double> getCentre();
    double getDeformation();
    double getRotation();
    double getLength();
    double getBreadth();

    // Miscellaneous output
    void writeMathematicaFormat(char *filename, const int xmax, const int ymax);
    void displayResult();
};

