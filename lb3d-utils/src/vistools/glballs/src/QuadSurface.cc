#include <cmath>

#include "T_CVector.h"
#include "GL_CQuadIO.h"
#include <vector>

char readBlock(FILE *InputStream,std::vector<T_CVector > &Vec,std::vector<T_CVector > &Col)
{
 // skip empty lines at the beginning
 char c;
 do {
  c=fgetc(InputStream);
 } while(c=='\n' || c==' ');
 ungetc(c,InputStream);

 Vec.clear();
 Col.clear();
 do {
  c=fgetc(InputStream);
  if(c!='\n' && c!=EOF){
   ungetc(c,InputStream);  
   T_R3_t x,y,z;
   char tmp;
   fscanf(InputStream,"%"T_R3_SCANF"%"T_R3_SCANF"%"T_R3_SCANF"%c",&x,&y,&z,&tmp);
   Vec.push_back(T_CVector(x,y,z));
   if(tmp=='#'){
    fscanf(InputStream,"%"T_R3_SCANF"%"T_R3_SCANF"%"T_R3_SCANF"%c",&x,&y,&z,&tmp);
    if(x<0)x=0;
    if(y<0)y=0;
    if(z<0)z=0;
    if(x>1.)x=1.;
    if(y>1.)y=1.;
    if(z>1.)z=1.;
    Col.push_back(T_CVector(x,y,z));
   } else {
    Col.push_back(T_CVector(0.,0.,0.));
   }
  }
 } while(c!='\n' && c!=EOF); 
 return c;
}

void writeTriangle(GL_CQuadIO &CQuadIO,T_CVector A,T_CVector B,T_CVector C,T_CVector Col)
{
 T_CVector BA=B-A;
 T_CVector CA=C-A;
 T_CVector N=BA.cross(CA); N.normalize();
 GL_CQuad CQuad;

 CQuad.set1(A);
 CQuad.set2(A);
 CQuad.set3(B);
 CQuad.set4(C);
 CQuad.setnorm(N);
 CQuad.setrgb((unsigned char)(256.*Col.getX()),(unsigned char)(256.*Col.getY()),(unsigned char)(256.*Col.getZ()));
 CQuadIO.writequad(CQuad);
}


int main( int argc, char **argv )
{
 if(argc!=2){
  return 1;
 }
 char *pString=(char *)malloc(sizeof(argv[1])+4);
 pString[0]=0;
 strcat(pString,argv[1]);
 strcat(pString,".qd");
 GL_CQuadIO CQuadIO(pString,"w");
 FILE *InputStream;
 InputStream=fopen(argv[1],"r");

 CQuadIO.writenumberofquads(0);

 std::vector<T_CVector> VecOld,VecNew;
 std::vector<T_CVector> ColOld,ColNew;

 readBlock(InputStream,VecNew,ColNew);

 int it=0;
 int imax=VecNew.size();
 char c;
 do {
  VecOld.swap(VecNew);
  ColOld.swap(ColNew);
  c=readBlock(InputStream,VecNew,ColNew);
  for(int i=0;i<imax-1;i++){
   T_CVector M=(VecNew[i]+VecNew[i+1]+VecOld[i]+VecOld[i+1])*0.25;
   T_CVector MCol=(ColNew[i]+ColNew[i+1]+ColOld[i]+ColOld[i+1])*0.25;
   writeTriangle(CQuadIO,VecNew[i+0],M,VecNew[i+1],(ColNew[i+0]+MCol+ColNew[i+1])*0.333333);
   writeTriangle(CQuadIO,VecNew[i+1],M,VecOld[i+1],(ColNew[i+1]+MCol+ColNew[i+1])*0.333333);
   writeTriangle(CQuadIO,VecOld[i+1],M,VecOld[i+0],(ColOld[i+1]+MCol+ColOld[i+0])*0.333333);
   writeTriangle(CQuadIO,VecOld[i+0],M,VecNew[i+0],(ColOld[i+0]+MCol+ColNew[i+0])*0.333333);
   it+=4;
  }
 }while(c!=EOF);

 CQuadIO.rewritenumberofquads(it);

 return 0;
}
