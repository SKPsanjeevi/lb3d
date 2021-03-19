#include "GL_CQuadIO.h"

int main( int argc, char **argv )
{
 GL_CQuadIO CQuadIO("plane.qd","w");

 int n;
 fscanf(stdin,"%i",&n);
 CQuadIO.writenumberofquads(n);
 GL_CQuad CQuad;

for(int i=0;i<n;i++){
 double x,y,z;
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.set1(x,y,z);
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.set2(x,y,z);
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.set3(x,y,z);
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.set4(x,y,z);
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.setnorm(x,y,z);
 fscanf(stdin,"%le%le%le",&x,&y,&z);
 CQuad.setrgb((char)(256*x),(char)(256*y),(char)(256*z));

 CQuadIO.writequad(CQuad);
}
/*
1
0 -1 -1   0 -1 1   0 1 1   0 1 -1   0 0 1   0 0 0

*/

 return 0;
}
