#include <iostream>
#include<iostream>
#include<fstream>
#include<string>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>


using namespace std;

int main()
{
    //calculates the contact angle of a particle adsorbed at the interface
    const int Gx=128;
    const int Gy=128;
    const int Gz=128;
    int gx=0;
    int gy=0;
    int gz=0;
    const int V=0;
    const int Nz=40;
    string datzeile;
    int Ni=0;
    bool anfang=false;
    bool inzahl=false;
    int nz=0;
    string Zahlstring[Nz];
    string ZahlstringCol[Gx][Gy][Gz];
    char Zahlchar;
    string ZstringCol;
    double Col[Gx][Gy][Gz];
    double Gf[Gx][Gz];
    bool posi=false;
    bool phase1=false;
    int lp=0;
    int en=0;
    double my=0;
    double cy=0;
    const string lb_gr_out_file="output";
    const string lbID="0861629126";
    const string tstep="00010000";
    string filemd="";
    string filecol="";
    string filemdpart1="";
    string filecolpart1="";
    string filepart2="";
    string filemdpart1a="md-cfg_";
    string filemdpart1b="_t";
    string filepartm="-";
    string fileend=".asc";
    string filecolpart1a="colour_";
    const int XNr=0;
    const int YNr=1;
    const int ZNr=2;
    const int OXNr=9;
    const int OYNr=10;
    const int OZNr=11;
    double X=0;
    double Y=0;
    double Z=0;
    double Ox=0;
    double Oy=0;
    double Oz=0;
    double D=0;
    double Dx=0;
    double Dy=0;
    double Dz=0;
    double Dp=0;
    double Do=0;
    bool inparticle=false;
    bool inexp=false;
    string Conv="";
    string Convz="";
    string Conve="";
    int Ne=0;
    const double Rp=10;
    const double Ro=10;
    double dGf=0;
    int Ndgf=0;
    double md=Rp+4;
    double theta=0;
    double dtheta=0;
    double thetag=0;
    double theta2=0;
    double theta2g=0;
    double myp=0;



    filemdpart1=filemdpart1a+lb_gr_out_file+filemdpart1b;
    filepart2=lbID+fileend;
    filecolpart1=filecolpart1a+lb_gr_out_file+filemdpart1b;
    filecol=filecolpart1+tstep+filepartm+filepart2;
    filemd=filemdpart1+tstep+filepartm+filepart2;
    //cout<<filecol<<"    "<<filemd<<endl;
    ifstream colfile(filecol.c_str());
    if (colfile.good())
    {
        gy=0;
        gz=V;
        colfile.seekg(0L, ios::beg);//cout<<"test"<<endl;
        while (! colfile.eof())//while(gx<Gx)
	{
	    getline(colfile,datzeile);//cout<<"test"<<endl;
	    Ni=int(datzeile.length())+1;//cout<<Ni<<endl;
            anfang=false;
            inzahl=false;
            nz=0;
            for(int i=0; i<Nz; i++)
	    {
	        Zahlstring[i]="";
            }
            for(int i=0; i<Ni; i++)
	    {
                Zahlchar=datzeile[i];
                if(Zahlchar==':')anfang=true;
                if(anfang)
                {
                    if(inzahl)if(Zahlchar!=',')Zahlstring[nz]+=datzeile[i];
                    if(Zahlchar==' ')inzahl=true;
                    //if(gx==Gx)inzahl=false;
                    //if(gx==Gx-1)if(gy==Gy-1)if(gx==Gx-1)inzahl=false;
                    if(Zahlchar==',')
   	            {
                        inzahl=false;
                        nz+=1;
                    }
                    //if(Zahlchar=='v')inzahl=false;
	        }
                
            }
            for(int i=0; i<(nz); i++)
            {
	        if(gy<Gy)
		{
                    if(gx<Gx)ZahlstringCol[gx][gy][gz]=Zahlstring[i];
                    //if(gx<Gx)cout<<gx<<"    "<<gy<<"    "<<gz<<"    "<<ZahlstringCol[gx][gy][gz]<<endl;
                }
                if(gy==Gy)cout<<"O: gy="<<Gy<<endl;
                gz+=1;
                if(gz==Gz)
	        {
	            gz=0;
		    gy+=1;
                }
                if(gy==Gy)
	        {
	            gy=0;
		    gx+=1;
                }
            }
        }
    }//cout<<"test"<<endl;
    for(int ix=0; ix<(Gx); ix++)
    {
        for(int iy=0; iy<(Gy); iy++)
	{
            for(int iz=0; iz<(Gz); iz++)
	    {
	        ZstringCol=ZahlstringCol[ix][iy][iz];
                std::stringstream strCol(ZstringCol);
                strCol >> Col[iz][iy][ix];
            }
        }
    }
    ifstream mdfile(filemd.c_str());
    if (mdfile.good())
    {
        mdfile.seekg(0L, ios::beg);
        {
            getline(mdfile,datzeile);
            Ni=int(datzeile.length())+1;
            //cout<<Ni<<" test"<<endl;
            inzahl=false;
            nz=0;
            for(int i=0; i<Nz; i++)
            {
                Zahlstring[i]="";
            }
            for(int i=0; i<Ni; i++)
            {
                Zahlchar=datzeile[i];
                if(Zahlchar=='+')inzahl=true;
                if(Zahlchar=='0')inzahl=true;
                if(Zahlchar=='-')inzahl=true;
                if(inzahl)if(Zahlchar!=' ')Zahlstring[nz]+=datzeile[i];
                //cout<<Zahlchar<<" test"<<endl;
                if(Zahlchar==' ')
   	        {
                    inzahl=false;
                    nz+=1;
                }
            }
        }//cout<<Zahlstring[XNr]<<endl;
        Conv=Zahlstring[XNr];
        Ni=int(Zahlstring[XNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";//cout<<"test"<<endl;
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
            {
                inzahl=false;
                inexp=true;
            }
        }//cout<<X<<endl;
        std::stringstream strConvex(Conve);
        strConvex >> Ne;
        std::stringstream strConvzx(Convz);
        strConvzx >> X;
        //cout<<X<<endl;
        if(Ne==1)X=X*10;
        if(Ne==2)X=X*100;
        if(Ne==3)X=X*1000;
        if(Ne==4)X=X*10000;
        if(Ne==-1)X=X*0.1;
        if(Ne==-2)X=X*0.01;
        if(Ne==-3)X=X*0.001;
        if(Ne==-4)X=X*0.0001;
        if(Ne==-5)X=X*0.00001;
        if(Ne==-6)X=X*0.000001;
        if(Ne==-7)X=X*0.0000001;
        if(Ne==-8)X=X*0.00000001;
        if(Ne==-9)X=X*0.000000001;
        if(Ne==-10)X=X*0.0000000001;
        if(Ne==-11)X=X*0.00000000001;
        if(Ne==-12)X=X*0.000000000001;
        if(Ne==-13)X=X*0.0000000000001;
        if(Ne==-14)X=X*0.00000000000001;
        if(Ne==-15)X=X*0.000000000000001;
        if(Ne==-16)X=X*0.0000000000000001;
        if(Ne<-16)X=0;
        //cout<<X<<endl;
        Conv=Zahlstring[YNr];
        Ni=int(Zahlstring[YNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
            {
                inzahl=false;
                inexp=true;
            }
        }
        std::stringstream strConvey(Conve);
        strConvey >> Ne;
        std::stringstream strConvzy(Convz);
        strConvzy >> Y;
        if(Ne==1)Y=Y*10;
        if(Ne==2)Y=Y*100;
        if(Ne==3)Y=Y*1000;
        if(Ne==4)Y=Y*10000;
        if(Ne==-1)Y=Y*0.1;
        if(Ne==-2)Y=Y*0.01;
        if(Ne==-3)Y=Y*0.001;
        if(Ne==-4)Y=Y*0.0001;
        if(Ne==-5)Y=Y*0.00001;
        if(Ne==-6)Y=Y*0.000001;
        if(Ne==-7)Y=Y*0.0000001;
        if(Ne==-8)Y=Y*0.00000001;
        if(Ne==-9)Y=Y*0.000000001;
        if(Ne==-10)Y=Y*0.0000000001;
        if(Ne==-11)Y=Y*0.00000000001;
        if(Ne==-12)Y=Y*0.000000000001;
        if(Ne==-13)Y=Y*0.0000000000001;
        if(Ne==-14)Y=Y*0.00000000000001;
        if(Ne==-15)Y=Y*0.000000000000001;
        if(Ne==-16)Y=Y*0.0000000000000001;
        if(Ne<-16)Y=0;
        Conv=Zahlstring[ZNr];
        Ni=int(Zahlstring[ZNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
            {
                inzahl=false;
                inexp=true;
            }
        }
        std::stringstream strConvez(Conve);
        strConvez >> Ne;
        std::stringstream strConvzz(Convz);
        strConvzz >> Z;
        if(Ne==1)Z=Z*10;
        if(Ne==2)Z=Z*100;
        if(Ne==3)Z=Z*1000;
        if(Ne==4)Z=Z*10000;
        if(Ne==-1)Z=Z*0.1;
        if(Ne==-2)Z=Z*0.01;
        if(Ne==-3)Z=Z*0.001;
        if(Ne==-4)Z=Z*0.0001;
        if(Ne==-5)Z=Z*0.00001;
        if(Ne==-6)Z=Z*0.000001;
        if(Ne==-7)Z=Z*0.0000001;
        if(Ne==-8)Z=Z*0.00000001;
        if(Ne==-9)Z=Z*0.000000001;
        if(Ne==-10)Z=Z*0.0000000001;
        if(Ne==-11)Z=Z*0.00000000001;
        if(Ne==-12)Z=Z*0.000000000001;
        if(Ne==-13)Z=Z*0.0000000000001;
        if(Ne==-14)Z=Z*0.00000000000001;
        if(Ne==-15)Z=Z*0.000000000000001;
        if(Ne==-16)Z=Z*0.0000000000000001;
        if(Ne<-16)Z=0;
        Conv=Zahlstring[OXNr];
        Ni=int(Zahlstring[OXNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
            {
                inzahl=false;
                inexp=true;
            }
        }//cout<<"test"<<endl;
        std::stringstream strConveox(Conve);
        strConveox >> Ne;
        std::stringstream strConvzox(Convz);
        strConvzox >> Ox;
        if(Ne==1)Ox=Ox*10;
        if(Ne==2)Ox=Ox*100;
        if(Ne==3)Ox=Ox*1000;
        if(Ne==4)Ox=Ox*10000;
        if(Ne==-1)Ox=Ox*0.1;
        if(Ne==-2)Ox=Ox*0.01;
        if(Ne==-3)Ox=Ox*0.001;
        if(Ne==-4)Ox=Ox*0.0001;
        if(Ne==-5)Ox=Ox*0.00001;
        if(Ne==-6)Ox=Ox*0.000001;
        if(Ne==-7)Ox=Ox*0.0000001;
        if(Ne==-8)Ox=Ox*0.00000001;
        if(Ne==-9)Ox=Ox*0.000000001;
        if(Ne==-10)Ox=Ox*0.0000000001;
        if(Ne==-11)Ox=Ox*0.00000000001;
        if(Ne==-12)Ox=Ox*0.000000000001;
        if(Ne==-13)Ox=Ox*0.0000000000001;
        if(Ne==-14)Ox=Ox*0.00000000000001;
        if(Ne==-15)Ox=Ox*0.000000000000001;
        if(Ne==-16)Ox=Ox*0.0000000000000001;
        if(Ne<-16)Ox=0;
        Conv=Zahlstring[OYNr];
        Ni=int(Zahlstring[OYNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
	    {
                inzahl=false;
                inexp=true;
            }
        }
        std::stringstream strConveoy(Conve);
        strConveoy >> Ne;
        std::stringstream strConvzoy(Convz);
        strConvzoy >> Oy;
        if(Ne==1)Oy=Oy*10;
        if(Ne==2)Oy=Oy*100;
        if(Ne==3)Oy=Oy*1000;
        if(Ne==4)Oy=Oy*10000;
        if(Ne==-1)Oy=Oy*0.1;
        if(Ne==-2)Oy=Oy*0.01;
        if(Ne==-3)Oy=Oy*0.001;
        if(Ne==-4)Oy=Oy*0.0001;
        if(Ne==-5)Oy=Oy*0.00001;
        if(Ne==-6)Oy=Oy*0.000001;
        if(Ne==-7)Oy=Oy*0.0000001;
        if(Ne==-8)Oy=Oy*0.00000001;
        if(Ne==-9)Oy=Oy*0.000000001;
        if(Ne==-10)Oy=Oy*0.0000000001;
        if(Ne==-11)Oy=Oy*0.00000000001;
        if(Ne==-12)Oy=Oy*0.000000000001;
        if(Ne==-13)Oy=Oy*0.0000000000001;
        if(Ne==-14)Oy=Oy*0.00000000000001;
        if(Ne==-15)Oy=Oy*0.000000000000001;
        if(Ne==-16)Oy=Oy*0.0000000000000001;
        if(Ne<-16)Oy=0;
        Conv=Zahlstring[OZNr];
        Ni=int(Zahlstring[OZNr].length())+1;
        inzahl=true;
        inexp=false;
        Convz="";
        Conve="";
        for(int i=0; i<Ni; i++)
        {
            Zahlchar=Conv[i];
            if(Zahlchar!='E')
            {
                if(inzahl)Convz+=Conv[i];
                if(inexp)Conve+=Conv[i];
            }
            if(Zahlchar=='E')
	    {
                inzahl=false;
                inexp=true;
            }
        }
        std::stringstream strConveoz(Conve);
        strConveoz >> Ne;
        std::stringstream strConvzoz(Convz);
        strConvzoz >> Oz;
        if(Ne==1)Oz=Oz*10;
        if(Ne==2)Oz=Oz*100;
        if(Ne==3)Oz=Oz*1000;
        if(Ne==4)Oz=Oz*10000;
        if(Ne==-1)Oz=Oz*0.1;
        if(Ne==-2)Oz=Oz*0.01;
        if(Ne==-3)Oz=Oz*0.001;
        if(Ne==-4)Oz=Oz*0.0001;
        if(Ne==-5)Oz=Oz*0.00001;
        if(Ne==-6)Oz=Oz*0.000001;
        if(Ne==-7)Oz=Oz*0.0000001;
        if(Ne==-8)Oz=Oz*0.00000001;
        if(Ne==-9)Oz=Oz*0.000000001;
        if(Ne==-10)Oz=Oz*0.0000000001;
        if(Ne==-11)Oz=Oz*0.00000000001;
        if(Ne==-12)Oz=Oz*0.000000000001;
        if(Ne==-13)Oz=Oz*0.0000000000001;
        if(Ne==-14)Oz=Oz*0.00000000000001;
        if(Ne==-15)Oz=Oz*0.000000000000001;
        if(Ne==-16)Oz=Oz*0.0000000000000001;
        if(Ne<-16)Oz=0;
    }//cout<<"test"<<endl;



    for(int iy=0; iy<(Gy); iy++)
    {
        for(int iz=0; iz<(Gz); iz++)
        {
	    posi=true;
            phase1=true;
            inparticle=false;
	    if(Col[2][iy][iz]<0)posi=false;
	    for(int ix=2; ix<(Gx-2); ix++)
	    {
	        Dx=ix-X;
                if(Dx<0)Dx=-Dx;
	        Dy=iy-Y;
                if(Dy<0)Dy=-Dy;
	        Dz=iz-Z;
                if(Dz<0)Dz=-Dz;
                D=sqrt(Dx*Dx+Dy*Dy+Dz*Dz);
                Dp=Dx*Ox+Dy*Oy+Dz*Oz;
                Do=sqrt(D*D-Dp*Dp);
                if(((Dp*Dp)/(Rp*Rp)+(Do*Do)/(Ro*Ro))<=1)inparticle=true;
                if(posi)
                {
                    if(Col[ix][iy][iz]<0)
                    {
                        if(phase1)en=ix;
                        phase1=false;
		    }
		    if(phase1)lp=ix;
		}
	        if(Col[ix][iy][iz]<0)
		{
		    if(posi)en=ix;
                    posi=false;
                }
	        if(posi)lp=ix;
            }
            my=Col[en][iy][iz]-Col[lp][iy][iz];
            cy=Col[lp][iy][iz]-my*lp;
            Gf[iy][iz]=-cy/my;
            if(inparticle)Gf[iy][iz]=0;
        }
    }
    //for(int iy=2; iy<(Gy-2); iy++)cout<<iy<<"    "<<Col[4][iy][4]<<endl;


    dGf=0;
    Ndgf=0;
    md=Rp+4;
    if(Ro>Rp)md=Ro+4;
    for(int iy=0; iy<(Gy); iy++)
    {
        for(int iz=0; iz<(Gz); iz++)
        {
	    Dx=iy-X;
            if(Dx<0)Dx=-Dx;
            Dz=iz-Z;
            if(Dz<0)Dz=-Dz;
            if(Dx>md)
	    {
	        dGf+=Gf[iy][iz];
                Ndgf+=1;
            }
            else if(Dz>md)
	    {
	        dGf+=Gf[iy][iz];
                Ndgf+=1;
            }
        }
    }
    dGf=dGf/Ndgf;
    cout<<"dGf="<<dGf<<"    Ndgf="<<Ndgf<<endl;
    cout<<"Y="<<Y<<endl;
    Dy=dGf-Y+1;
    cout<<"Dy="<<Dy<<endl;
    if(Dy<0)Dy=-Dy;
    cout<<"Dy="<<Dy<<endl;
    //theta=asin(Dy/Ro);
    myp=-((Ro*Dy)/(Rp*Rp))/sqrt(1-((Dy*Dy)/(Rp*Rp)));
    dtheta=atan(myp);
    theta=(M_PI/2.)+dtheta;
    thetag=theta*180/M_PI;
    cout<<"Angle (radian) theta="<<theta<<endl;
    cout<<"angle (in degree): theta="<<thetag<<endl;
    theta2=M_PI-theta;
    theta2g=180-thetag;
    cout<<"angle II (radian) theta2="<<theta2<<endl;
    cout<<"angle II (in degree): theta2="<<theta2g<<endl;






    return 0;
} 
