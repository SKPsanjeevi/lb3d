#ifndef CARDINALDIRECTION_H
#define CARDINALDIRECTION_H
#include "cardinaldirection.h"
#endif

using namespace std;

int find_extentminz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
int find_extentmaxz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
int find_topz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
int find_botz(float*** N,const int dx, const int dy, const int dz, const double cutoff);
double gammadot(float***N,const int dx, const int j, const int dz, const int ox, const int oy, const int oz);
