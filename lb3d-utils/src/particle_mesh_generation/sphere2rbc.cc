#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

#define SQ(x) ((x) * (x))

using namespace std;

/************************************************************
************ FUNCTION DECLARATIONS **************************
************************************************************/

double xy(double);
double z(double, double);

/************************************************************
************ VARIABLE DECLARATIONS **************************
************************************************************/

/// parameters for RBC shape

double C_0 = 0.207;
double C_1 = 2.0;
double C_2 = -1.123;
double R = 1;

/// other parameters

int num_lines = 0; // number of lines in input file
string dummy_string; // dummy string for reading data
int vertex; // vertex index
double X, Y, Z; // coordinates from input file
double XX, YY, ZZ; // coordinates for output file

/************************************************************
************ MAIN FUNCTION **********************************
************************************************************/

int main() {

  /// open input and output files

  ifstream file_i("sphere.msh");
  ofstream file_o("rbc.msh");

  /// get number of lines from input file

  while(getline(file_i, dummy_string)) {
    num_lines += 1;
  }

  cout << num_lines << " lines in the input file" << endl << endl;

  file_i.close();

  /// write output file

  file_i.open("sphere.msh");

  for(int i = 0; i < num_lines; ++i) {

    /// read coordinates from input file

    file_i >> vertex >> X >> Y >> Z;

    /// compute coordinates for output file

    //XX = xy(X);
    //YY = xy(Y);
    ZZ = z(X, Y);

    Z < 0 ? ZZ *= -1 : ZZ = ZZ;

    /// write coordinates to output file

    file_o << setprecision(15);
    file_o << vertex << " " << X * R << " " << Y * R << " ";

    if(isnan(ZZ) || abs(ZZ) < 1e-4) {
      file_o << 0 << "\n";
    }
    else {
      file_o << ZZ * R << "\n";
    }
  }

  file_i.close();

  return 0;
}

/************************************************************
************ MODIFY COORDINATES *****************************
************************************************************/

double xy(double x) {

  if(x > 0) {
    return (0.01 * x * x + 0.99 * x);
  }
  else {
    return (-0.01 * x * x + 0.99 * x);
  }
}

/***********************************************************/

double z(double x, double y) {

  return (0.5 * sqrt(1. - (SQ(x) + SQ(y)))
      * (C_0 + C_1 * (SQ(x) + SQ(y)) + C_2 * SQ(SQ(x) + SQ(y))));
}
