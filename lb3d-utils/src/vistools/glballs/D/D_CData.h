#if !defined(D_CDATA_H)
#define D_CDATA_H
#include "T.h"
#include "MD.h"

class D_CData {
 public:
  T_R3_t x,y,z,vx,vy,vz,r;
  T_R3_t qw,qx,qy,qz,wx,wy,wz;
  T_Id_t id;
  T_ParticleClassExtern_t type;
  MD_Color_t color;
};
typedef class D_CData D_Data_t;

class D_CHeader {
 public:
 MD_Particle_Counter_t Number;
 T_Time_t Time;
 T_R3_t xMin,yMin,zMin,xMax,yMax,zMax;
};
typedef class D_CHeader D_Header_t;

class D_CDataLine {
 public:
  T_R3_t x1,y1,z1,x2,y2,z2;
  MD_Color_t color;
};
typedef class D_CDataLine D_DataLine_t;

class D_CHeaderLine {
 public:
 MD_Particle_Counter_t Number;
};
typedef class D_CHeaderLine D_HeaderLine_t;

#endif
