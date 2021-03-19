#include "T_math.h"
#include <stdio.h>

int main(int argc,char *argv[]){
 double t=0.;
 double dt=0.1;
 double x,y,z,scale;
 double Cos,Sin;
 double theta,dtheta;

 theta=0.;
 dtheta=0.05;
 x=0.;
 y=1.;
 z=0.;

 scale=1./sqrt(x*x+y*y+z*z);
 x*=scale; y*=scale; z*=scale;
 
 int i;
 int iiimax=1;
for(i=0;i<1000;++i){
 fprintf(stdout,"%i %f -2 -2 -2 2 2 2\n",6+iiimax,t);
 Cos=std::cos(0.5*theta);
 Sin=std::sin(0.5*theta);
 fprintf(stdout,"1. 0. 0. 0. 0. 0. 0.4 %f %f %f %f 0 0 0 1 1 150\n",Cos,Sin*x,Sin*y,Sin*z);
 fprintf(stdout,"%f %f 0. %f %f 0. 0.2 %f 0 0 %f 0 0 0 2 1 150\n",std::cos(t),std::sin(t),-std::sin(t),std::cos(t),std::cos(0.5*t),std::sin(0.5*t));
 fprintf(stdout,"0. %f 0. 0. %f 0. 0.2 0 1 0 0 0 0 0 3 1 50\n",2.*std::sin(t),2*std::cos(t));
 fprintf(stdout,"0. 0.5 0. 0. 0. 0. 0.2 0 1 0 0 0 0 0 4 2 50\n");
 fprintf(stdout,"0. %f %f 0. %f %f 0.2 0 1 0 0 0 0 0 5 2 10\n",std::cos(t),std::sin(t),-std::sin(t),std::cos(t));
 fprintf(stdout,"0. 0. %f 0. 0. %f 0.2 0 1 0 0 0 0 0 6 2 10240\n",2.*std::sin(t),2*std::cos(t));
// fprintf(stdout,"#0. 0. %f 0. 0. %f 0.2 0 1 0 0 0 0 0 6 2 10240\n",2.*std::sin(t),2*std::cos(t));
for(int iii=0;iii<iiimax;iii++){
 fprintf(stdout,"0. 0. %f 0. 0. %f 0.2 0 1 0 0 0 0 0 6 2 10240\n",2.*std::sin(t),2*std::cos(t));
}
 fprintf(stdout,"#L 2\n");
 fprintf(stdout,"# 0. 0. %f 0. 0. %f %i\n",-2.*std::sin(t),-2*std::cos(t),(int)(256*std::sin(t)+256*256*std::cos(t)));
 fprintf(stdout,"# 0. %f %f 0. %f %f %i\n",std::cos(t),std::sin(t),-std::sin(t),std::cos(t),256*256*255+256*255+255);
 fprintf(stdout,"#\n");
 fflush(stdout);
//sleep(10);
 t+=dt;
 theta+=dtheta;
}
 return 0;
}
