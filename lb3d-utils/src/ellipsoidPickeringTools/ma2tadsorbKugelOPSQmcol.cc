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
    const int Nz=40;
    const int Nta=0;
    const int Nt=10000;
    //const int Nta=0;
    //const int Nt=200;
    //const int Nt=0;
    const int NT=24;
    const int dt=100;
    string datzeile;
    int Ni=0;
    bool inzahl=false;
    bool inexp=false;
    int nz=0;
    string Zahlstring[Nz][NT+1];
    string Zahlstringcol[Nz];
    char Zahlchar;
    string datei1st="md-cfg_absorb01_t00000000-0705201734.asc";
    const int TNr=15;
    const int XNr=0;
    const int YNr=1;
    const int ZNr=2;
    const int OXNr=9;
    const int OYNr=10;
    const int OZNr=11;
    string lb_gr_out_file="droplet";
    string lbID="1780644712";
    //const string lbID="0825974025";
    //const string lb_gr_out_file="";
    //const string lbID="";
    string fileend=".asc";
    string filepart1="";
    string filecolpart1="";
    string filepart1a="md-cfg_";
    string filecolpart1a="colour_";
    string filepart1b="_t";
    string filepart2="";
    string outfilestr="";
    string outcolfilestr="";
    double X[NT+1];
    double Y[NT+1];
    double Z[NT+1];
    double Ox[NT+1];
    double Oy[NT+1];
    double Oz[NT+1];
    int Ne=0;
    string Conv="";
    string Convz="";
    string Conve="";
    double YG[Nt+1];
    double YGi[Nt+1];
    const double Gi=32;
    double G=32;
    const double Rp=12;
    const double Ro=6;
    double dGf[Nt+1];
    int Ndgf=0;
    double md=Rp+4;
    const int Gx=64;
    const int Gy=64;
    const int Gz=64;
    int gx=0;
    int gy=0;
    int gz=0;
    const int V=0;
    bool anfang=false;
    bool inbereich=false;
    string ZahlstringCol[Gx][Gy][Gz];
    string ZstringCol;
    double Col[Gx][Gy][Gz];
    double Gf[Gx][Gz];
    bool posi=false;
    int lp=0;
    int en=0;
    double my=0;
    double cy=0;
    double Dx=0;
    double Dy=0;
    double Dz=0;
    double S[Nt+1];
    double Q[Nt+1];
    double Qre[Nt+1];
    double Qn1[Nt+1];
    double Qn2[Nt+1];
    double Qn3[Nt+1];
    double Qn4[Nt+1];
    double Qn5[Nt+1];
    double Qmax[Nt+1];
    //double theta[NT+1];
    //double phi[NT+1];
    double thetaa[NT+1];
    double phia[NT+1];
    double thetar[NT+1];
    double phir[NT+1];
    double thetaK[NT+1];
    double phiK[NT+1];
    double phire[NT+1];
    double phineu[NT+1];
    double phineu1[NT+1];
    double phineu2[NT+1];
    double phineu3[NT+1];
    double phineu4[NT+1];
    double phineu5[NT+1];
    const int Nth=1000;
    //const int Nth=10;
    const int Nph=2000;
    //const int Nph=20;
    const int Nph2=20;
    const int Nph3=10;
    const int Nph4=2000;
    //const int Nph4=2;
    //const int Nph4=20;
    double rho[Nth+1][Nph];
    double rhore[Nth+1][Nph];
    double rhoneu1[Nth+1][Nph];
    double rhoneu2[Nth+1][Nph];
    double rhoneu3[Nth+1][Nph];
    double rhoneu4[Nth+1][Nph];
    double rhoneu5[Nth+1][Nph];
    //int rho[Nth+1][Nph];??
    int rphi[Nph2];
    int rphi2[Nph3];
    int rphimax;
    int nrphimax;
    int nrphimaxm;
    int nrphimaxp;
    double avphi=0;
    int navphi=0;
    const int rphic=NT/8;
    const int rphic2=NT/16;
    const double dth=M_PI/Nth;
    const double dph=2*M_PI/Nph;
    const double dph2=2*M_PI/Nph2;
    const double dph3=2*M_PI/Nph3;
    const double dph4=2*M_PI/Nph4;
    double rhodom;
    double rhoredom;
    double rhoneu1dom;
    double rhoneu2dom;
    double rhoneu3dom;
    double rhoneu4dom;
    double rhoneu5dom;
    double Sint;
    double Qint;
    double dphiQmax[Nt+1];
    //const double Cphi4=2.19911485751;
    //const double Cphi4=2.82743338823;
    const double Cphi4=0.6*M_PI;
    double ColM=0;
    double Colrsx=0;
    double Colrsy=0;
    double Colrsz=0;
    int ColrsxI=0;
    int ColrsyI=0;
    int ColrszI=0;
    double rKPx[NT+1];
    double rKPy[NT+1];
    double rKPz[NT+1];
    double rKP[NT+1];
    double eKPx[NT+1];
    double eKPy[NT+1];
    double eKPz[NT+1];
    double eKthx[NT+1];
    double eKthy[NT+1];
    double eKthz[NT+1];
    double OaeKR=0;
    double OaeKth=0;
    double NDx=0;
    double NDy=0;
    double NDz=0;
    double DR=0;
    double NDmax=0;
    double DDef=0;
    int NTTO=0;
    const double dDR=4;
    bool anGf[NT+1];
    double phin[Nt+1];




    /*ifstream fileidc("id.conf");
    if (fileidc.good())
    {
	int i=0;
        while (! fileidc.eof())
        {
            getline(fileidc,datzeile);
            if(i==0)lbID=datzeile;
            i++;
        }
    }
    ifstream fileonc("outname.conf");
    if (fileonc.good())
    {
	int i=0;
        while (! fileonc.eof())
        {
            getline(fileonc,datzeile);
            if(i==0)lb_gr_out_file=datzeile;
            i++;
        }
    }*/
    lb_gr_out_file="droplet";
    lbID="1780644712";


    filepart1=filepart1a+lb_gr_out_file+filepart1b;
    filepart2=lbID+fileend;
    filecolpart1=filecolpart1a+lb_gr_out_file+filepart1b;
    fstream outOP;
    //outOP.open("OPSQ", ios::out);
    outOP.open("OPSQm", ios::out);
    fstream outStatus;
    outStatus.open("StatusOP", ios::out);
    for(int t=Nta; t<(Nt+1); t++)
    {
        if(t%dt==0)
	  {
            string st;
            stringstream sst;
            string Nl;
            if (t<10)
	    {
	        Nl="0000000";
            }
            else if(t<100)
	    {
                Nl="000000";
            }
            else if(t<1000)
	    {
                Nl="00000";
            }
            else if(t<10000)
	    {
                Nl="0000";
            }
            else if(t<100000)
	    {
                Nl="000";
            }
            if (sst << t)
            {
                std::string ssst(sst.str());
                st=ssst;
            }
            outfilestr=filepart1+Nl+st+"-"+filepart2;
            outcolfilestr=filecolpart1+Nl+st+"-"+filepart2;
            char *outfile;
            Ni=int(outfilestr.length())+1;
            for(int i=0; i<Ni; i++)outfile+=outfilestr[i];
            ifstream file1st(outfilestr.c_str());
            //cout<<outfilestr<<endl;
            if (file1st.good())
            {
	        int nT=0;
                file1st.seekg(0L, ios::beg);
                while (! file1st.eof())
                {
                    getline(file1st,datzeile);
                    Ni=int(datzeile.length())+1;
                    inzahl=false;
                    nz=0;
                    for(int i=0; i<Nz; i++)
    	            {
	                Zahlstring[i][nT]="";
                    }
                    for(int i=0; i<Ni; i++)
	            {
                        Zahlchar=datzeile[i];
                        if(Zahlchar=='+')inzahl=true;
                        if(Zahlchar=='0')inzahl=true;
                        if(Zahlchar=='-')inzahl=true;
                        if(inzahl)if(Zahlchar!=' ')Zahlstring[nz][nT]+=datzeile[i];
                        if(Zahlchar==' ')
   	                {
                            inzahl=false;
                            nz+=1;
                        }
                        //nT+=1;
                    }
                    nT+=1;
                }
                //evtl folgendes in eigene Funktion auslagern
 



                for(int nt=0; nt<(NT); nt++)
		{
                    Conv=Zahlstring[XNr][nt];
                    Ni=int(Zahlstring[XNr][nt].length())+1;
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
                    std::stringstream strConvxe(Conve);
                    strConvxe >> Ne;
                    std::stringstream strConvxz(Convz);
                    strConvxz >> X[nt];
                    if(Ne==1)X[nt]=X[nt]*10;
                    if(Ne==2)X[nt]=X[nt]*100;
                    if(Ne==3)X[nt]=X[nt]*1000;
                    if(Ne==4)X[nt]=X[nt]*10000;
                    if(Ne==-1)X[nt]=X[nt]*0.1;
                    if(Ne==-2)X[nt]=X[nt]*0.01;
                    if(Ne==-3)X[nt]=X[nt]*0.001;
                    if(Ne==-4)X[nt]=X[nt]*0.0001;
                    if(Ne==-5)X[nt]=X[nt]*0.00001;
                    if(Ne==-6)X[nt]=X[nt]*0.000001;
                    if(Ne==-7)X[nt]=X[nt]*0.0000001;
                    if(Ne==-8)X[nt]=X[nt]*0.00000001;
                    if(Ne==-9)X[nt]=X[nt]*0.000000001;
                    if(Ne==-10)X[nt]=X[nt]*0.0000000001;
                    if(Ne==-11)X[nt]=X[nt]*0.00000000001;
                    if(Ne==-12)X[nt]=X[nt]*0.000000000001;
                    if(Ne==-13)X[nt]=X[nt]*0.0000000000001;
                    if(Ne==-14)X[nt]=X[nt]*0.00000000000001;
                    if(Ne==-15)X[nt]=X[nt]*0.000000000000001;
                    if(Ne==-16)X[nt]=X[nt]*0.0000000000000001;
                    if(Ne<-16)X[nt]=0;
                    





                    Conv=Zahlstring[YNr][nt];
                    Ni=int(Zahlstring[YNr][nt].length())+1;
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
                    std::stringstream strConvye(Conve);
                    strConvye >> Ne;
                    std::stringstream strConvyz(Convz);
                    strConvyz >> Y[nt];
                    if(Ne==1)Y[nt]=Y[nt]*10;
                    if(Ne==2)Y[nt]=Y[nt]*100;
                    if(Ne==3)Y[nt]=Y[nt]*1000;
                    if(Ne==4)Y[nt]=Y[nt]*10000;
                    if(Ne==-1)Y[nt]=Y[nt]*0.1;
                    if(Ne==-2)Y[nt]=Y[nt]*0.01;
                    if(Ne==-3)Y[nt]=Y[nt]*0.001;
                    if(Ne==-4)Y[nt]=Y[nt]*0.0001;
                    if(Ne==-5)Y[nt]=Y[nt]*0.00001;
                    if(Ne==-6)Y[nt]=Y[nt]*0.000001;
                    if(Ne==-7)Y[nt]=Y[nt]*0.0000001;
                    if(Ne==-8)Y[nt]=Y[nt]*0.00000001;
                    if(Ne==-9)Y[nt]=Y[nt]*0.000000001;
                    if(Ne==-10)Y[nt]=Y[nt]*0.0000000001;
                    if(Ne==-11)Y[nt]=Y[nt]*0.00000000001;
                    if(Ne==-12)Y[nt]=Y[nt]*0.000000000001;
                    if(Ne==-13)Y[nt]=Y[nt]*0.0000000000001;
                    if(Ne==-14)Y[nt]=Y[nt]*0.00000000000001;
                    if(Ne==-15)Y[nt]=Y[nt]*0.000000000000001;
                    if(Ne==-16)Y[nt]=Y[nt]*0.0000000000000001;
                    if(Ne<-16)Y[nt]=0;




                    Conv=Zahlstring[ZNr][nt];
                    Ni=int(Zahlstring[ZNr][nt].length())+1;
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
                    std::stringstream strConvze(Conve);
                    strConvze >> Ne;
                    std::stringstream strConvzz(Convz);
                    strConvzz >> Z[nt];
                    if(Ne==1)Z[nt]=Z[nt]*10;
                    if(Ne==2)Z[nt]=Z[nt]*100;
                    if(Ne==3)Z[nt]=Z[nt]*1000;
                    if(Ne==4)Z[nt]=Z[nt]*10000;
                    if(Ne==-1)Z[nt]=Z[nt]*0.1;
                    if(Ne==-2)Z[nt]=Z[nt]*0.01;
                    if(Ne==-3)Z[nt]=Z[nt]*0.001;
                    if(Ne==-4)Z[nt]=Z[nt]*0.0001;
                    if(Ne==-5)Z[nt]=Z[nt]*0.00001;
                    if(Ne==-6)Z[nt]=Z[nt]*0.000001;
                    if(Ne==-7)Z[nt]=Z[nt]*0.0000001;
                    if(Ne==-8)Z[nt]=Z[nt]*0.00000001;
                    if(Ne==-9)Z[nt]=Z[nt]*0.000000001;
                    if(Ne==-10)Z[nt]=Z[nt]*0.0000000001;
                    if(Ne==-11)Z[nt]=Z[nt]*0.00000000001;
                    if(Ne==-12)Z[nt]=Z[nt]*0.000000000001;
                    if(Ne==-13)Z[nt]=Z[nt]*0.0000000000001;
                    if(Ne==-14)Z[nt]=Z[nt]*0.00000000000001;
                    if(Ne==-15)Z[nt]=Z[nt]*0.000000000000001;
                    if(Ne==-16)Z[nt]=Z[nt]*0.0000000000000001;
                    if(Ne<-16)Z[nt]=0;




                    Conv=Zahlstring[OXNr][nt];
                    Ni=int(Zahlstring[OXNr][nt].length())+1;
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
                    std::stringstream strConvex(Conve);
                    strConvex >> Ne;
                    std::stringstream strConvzx(Convz);
                    strConvzx >> Ox[nt];
                    if(Ne==1)Ox[nt]=Ox[nt]*10;
                    if(Ne==2)Ox[nt]=Ox[nt]*100;
                    if(Ne==3)Ox[nt]=Ox[nt]*1000;
                    if(Ne==4)Ox[nt]=Ox[nt]*10000;
                    if(Ne==-1)Ox[nt]=Ox[nt]*0.1;
                    if(Ne==-2)Ox[nt]=Ox[nt]*0.01;
                    if(Ne==-3)Ox[nt]=Ox[nt]*0.001;
                    if(Ne==-4)Ox[nt]=Ox[nt]*0.0001;
                    if(Ne==-5)Ox[nt]=Ox[nt]*0.00001;
                    if(Ne==-6)Ox[nt]=Ox[nt]*0.000001;
                    if(Ne==-7)Ox[nt]=Ox[nt]*0.0000001;
                    if(Ne==-8)Ox[nt]=Ox[nt]*0.00000001;
                    if(Ne==-9)Ox[nt]=Ox[nt]*0.000000001;
                    if(Ne==-10)Ox[nt]=Ox[nt]*0.0000000001;
                    if(Ne==-11)Ox[nt]=Ox[nt]*0.00000000001;
                    if(Ne==-12)Ox[nt]=Ox[nt]*0.000000000001;
                    if(Ne==-13)Ox[nt]=Ox[nt]*0.0000000000001;
                    if(Ne==-14)Ox[nt]=Ox[nt]*0.00000000000001;
                    if(Ne==-15)Ox[nt]=Ox[nt]*0.000000000000001;
                    if(Ne==-16)Ox[nt]=Ox[nt]*0.0000000000000001;
                    if(Ne<-16)Ox[nt]=0;
                    Conv=Zahlstring[OYNr][nt];
                    Ni=int(Zahlstring[OYNr][nt].length())+1;
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
                    strConvzy >> Oy[nt];
                    if(Ne==1)Oy[nt]=Oy[nt]*10;
                    if(Ne==2)Oy[nt]=Oy[nt]*100;
                    if(Ne==3)Oy[nt]=Oy[nt]*1000;
                    if(Ne==4)Oy[nt]=Oy[nt]*10000;
                    if(Ne==-1)Oy[nt]=Oy[nt]*0.1;
                    if(Ne==-2)Oy[nt]=Oy[nt]*0.01;
                    if(Ne==-3)Oy[nt]=Oy[nt]*0.001;
                    if(Ne==-4)Oy[nt]=Oy[nt]*0.0001;
                    if(Ne==-5)Oy[nt]=Oy[nt]*0.00001;
                    if(Ne==-6)Oy[nt]=Oy[nt]*0.000001;
                    if(Ne==-7)Oy[nt]=Oy[nt]*0.0000001;
                    if(Ne==-8)Oy[nt]=Oy[nt]*0.00000001;
                    if(Ne==-9)Oy[nt]=Oy[nt]*0.000000001;
                    if(Ne==-10)Oy[nt]=Oy[nt]*0.0000000001;
                    if(Ne==-11)Oy[nt]=Oy[nt]*0.00000000001;
                    if(Ne==-12)Oy[nt]=Oy[nt]*0.000000000001;
                    if(Ne==-13)Oy[nt]=Oy[nt]*0.0000000000001;
                    if(Ne==-14)Oy[nt]=Oy[nt]*0.00000000000001;
                    if(Ne==-15)Oy[nt]=Oy[nt]*0.000000000000001;
                    if(Ne==-16)Oy[nt]=Oy[nt]*0.0000000000000001;
                    if(Ne<-16)Oy[nt]=0;
                    Conv=Zahlstring[OZNr][nt];
                    Ni=int(Zahlstring[OZNr][nt].length())+1;
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
                    std::stringstream strConvzoz(Convz);
                    strConvzoz >> Oz[nt];
                    if(Ne==1)Oz[nt]=Oz[nt]*10;
                    if(Ne==2)Oz[nt]=Oz[nt]*100;
                    if(Ne==3)Oz[nt]=Oz[nt]*1000;
                    if(Ne==4)Oz[nt]=Oz[nt]*10000;
                    if(Ne==-1)Oz[nt]=Oz[nt]*0.1;
                    if(Ne==-2)Oz[nt]=Oz[nt]*0.01;
                    if(Ne==-3)Oz[nt]=Oz[nt]*0.001;
                    if(Ne==-4)Oz[nt]=Oz[nt]*0.0001;
                    if(Ne==-5)Oz[nt]=Oz[nt]*0.00001;
                    if(Ne==-6)Oz[nt]=Oz[nt]*0.000001;
                    if(Ne==-7)Oz[nt]=Oz[nt]*0.0000001;
                    if(Ne==-8)Oz[nt]=Oz[nt]*0.00000001;
                    if(Ne==-9)Oz[nt]=Oz[nt]*0.000000001;
                    if(Ne==-10)Oz[nt]=Oz[nt]*0.0000000001;
                    if(Ne==-11)Oz[nt]=Oz[nt]*0.00000000001;
                    if(Ne==-12)Oz[nt]=Oz[nt]*0.000000000001;
                    if(Ne==-13)Oz[nt]=Oz[nt]*0.0000000000001;
                    if(Ne==-14)Oz[nt]=Oz[nt]*0.00000000000001;
                    if(Ne==-15)Oz[nt]=Oz[nt]*0.000000000000001;
                    if(Ne==-16)Oz[nt]=Oz[nt]*0.0000000000000001;
                    if(Ne<-16)Oz[nt]=0;
                    thetaa[nt]=acos(Oz[nt]);
                    phia[nt]=asin(Oy[nt]/sin(thetaa[nt]));
                    if(thetaa[nt]>M_PI/2)thetaa[nt]=M_PI-thetaa[nt];
                    phia[nt]+=M_PI;
                    //cout<<t<<"    "<<nt<<"    "<<Oz[nt]<<"    "<<thetaa[nt]<<"    "<<phia[nt]<<endl;
                }
                
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rho[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
                    if(anGf[nT])
                    {
  	 	        if(thetaa[nT]>=ith*dth)if(thetaa[nT]<(ith+1)*dth)if(phia[nT]>=iph*dph)if(phia[nT]<(iph+1)*dph)
	  	        {
	  	            rho[ith][iph]+=1;
                        }
                    }
                }
                //if(rho[ith][iph]!=0)cout<<t<<"    "<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            ifstream colfile(outcolfilestr.c_str());
            if (colfile.good())
            {
	        gx=0;
	        gy=0;
                gz=V;
                colfile.seekg(0L, ios::beg);
                inbereich=true;
                while (! colfile.eof())
		{
                    getline(colfile,datzeile);
                    Ni=int(datzeile.length())+1;
                    anfang=false;
                    inzahl=false;
                    nz=0;
                    for(int i=0; i<Nz; i++)
	            {
	                Zahlstringcol[i]="";
                    }
                    for(int i=0; i<Ni; i++)
	            {
                        Zahlchar=datzeile[i];
                        if(Zahlchar==':')anfang=true;
                        if(Zahlchar=='U')anfang=false;
                        if(Zahlchar=='L')anfang=false;
                        if(Zahlchar=='I')anfang=false;
                        if(Zahlchar=='i')anfang=false;
                        if(Zahlchar=='M')anfang=false;
                        if(Zahlchar=='N')anfang=false;
                        if(Zahlchar=='E')anfang=false;
                        if(Zahlchar=='e')anfang=false;
                        if(Zahlchar=='m')inbereich=false;
                        if(anfang)if(inbereich)
                        {
                            if(inzahl)if(Zahlchar!=',')Zahlstringcol[nz]+=datzeile[i];
                            if(Zahlchar==' ')inzahl=true;
                            if(Zahlchar==',')
   	                    {
                                inzahl=false;
                                nz+=1;
                            }
	                }
                    }
                    for(int i=0; i<(nz); i++)
		    {
	                if(gy<Gy)
		        {
                            ZahlstringCol[gx][gy][gz]=Zahlstringcol[i];
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
            }

            for(int ix=0; ix<(Gx); ix++)
            {
                for(int iy=0; iy<(Gy); iy++)
                {
                    for(int iz=0; iz<(Gz); iz++)
                    {
	                ZstringCol=ZahlstringCol[ix][iy][iz];
                        std::stringstream strCol(ZstringCol);
                        strCol >> Col[ix][iy][iz];
                    }
                }
            }
            for(int ix=0; ix<(Gx); ix++)
            {
                for(int iz=0; iz<(Gz); iz++)
                {
	            posi=true;
                    for(int iy=0; iy<(Gy); iy++)
                    {
	                if(Col[ix][iy][iz]<0)
                        {
		            if(posi)en=iy;
                            posi=false;
                        }
	                if(posi)lp=iy;
                    }
                    my=Col[ix][en][iz]-Col[ix][lp][iz];
                    cy=Col[ix][lp][iz]-my*lp;
                    Gf[ix][iz]=-cy/my;
                }
            }
            cout<<"Dateien ausgelesen"<<endl;
            outStatus<<"Dateien ausgelesen"<<endl;

            //cout<<Col[0][0][0]<<"    "<<Col[32][32][32]<<"    "<<outcolfilestr<<endl;
            ColM=0;
            Colrsx=0;
            Colrsy=0;
            Colrsz=0;
            for(int ix=0; ix<(Gx); ix++)
            {
                for(int iz=0; iz<(Gz); iz++)
                {
                    for(int iy=0; iy<(Gy); iy++)
                    {
	                if(Col[ix][iy][iz]>0)
                        {
                            ColM+=1;
                            Colrsx+=ix;
                            Colrsy+=iy;
                            Colrsz+=iz;
                        }
                    }
                }
            }
            Colrsx=Colrsx/ColM;
            Colrsy=Colrsy/ColM;
            Colrsz=Colrsz/ColM;
            //cout<<t<<"    "<<Colrsx<<"    "<<Colrsy<<"    "<<Colrsz<<endl;
            ColrsxI=int(Colrsx+0.5);
            ColrsyI=int(Colrsy+0.5);
            ColrszI=int(Colrsz+0.5);
            //cout<<ColrsxI<<"    "<<ColrsyI<<"    "<<ColrszI<<endl;
            NDx=0;
            NDy=0;
            NDz=0;
            for(int ix=0; ix<(Gx); ix++)
	    {
                if(Col[ix][ColrsyI][ColrszI]>0)NDx+=1;
            }
            for(int iy=0; iy<(Gy); iy++)
	    {
                if(Col[ColrsxI][iy][ColrszI]>0)NDy+=1;
            }
            for(int iz=0; iz<(Gz); iz++)
	    {
                if(Col[ColrsxI][ColrsyI][iz]>0)NDz+=1;
            }
            DR=(NDx+NDy+NDz)/6;
            NDmax=NDx;
            if(NDmax<NDy)NDmax=NDy;
            if(NDmax<NDz)NDmax=NDz;
            DDef=(((NDmax-NDx)/(NDmax+NDx))+((NDmax-NDy)/(NDmax+NDy))+((NDmax-NDz)/(NDmax+NDz)))/2.;
            NTTO=0;
            for(int nT=0; nT<NT; nT++)
	    {
                rKPx[nT]=X[nT]-Colrsx;
                rKPy[nT]=Y[nT]-Colrsy;
                rKPz[nT]=Z[nT]-Colrsz;
                rKP[nT]=sqrt(rKPx[nT]*rKPx[nT]+rKPy[nT]*rKPy[nT]+rKPz[nT]*rKPz[nT]);
                anGf[nT]=false;
                //if((DR-dDR)>rKP[nT])cout<<nT<<"    kleiner    "<<DR<<"    "<<rKP[nT]<<endl;
                //if((DR+dDR)<rKP[nT])cout<<nT<<"    größer    "<<DR<<"    "<<rKP[nT]<<endl;
                //cout<<t<<"    "<<nT<<"    "<<DR<<"    "<<rKP[nT]<<endl;
                if((DR-dDR)<rKP[nT])if((DR+dDR)>rKP[nT])
		{
                    eKPx[nT]=rKPx[nT]/rKP[nT];
                    eKPy[nT]=rKPy[nT]/rKP[nT];
                    eKPz[nT]=rKPz[nT]/rKP[nT];
                    OaeKR=Ox[nT]*eKPx[nT]+Oy[nT]*eKPy[nT]+Oz[nT]*eKPz[nT];
                    thetar[nT]=acos(OaeKR);
                    thetaK[nT]=acos(eKPz[nT]);
                    phiK[nT]=asin(eKPy[nT]/sin(thetaK[nT]));
                    if(thetaK[nT]>M_PI/2)thetaK[nT]=M_PI-thetaK[nT];
                    phiK[nT]+=M_PI;
                    //cout<<nT<<"    "<<thetaK[nT]<<"    "<<phiK[nT]<<endl;
                    //cout<<t<<"    "<<nT<<"    "<<DR<<"    "<<rKP[nT]<<endl;
                    eKthx[nT]=cos(thetaK[nT])*cos(phiK[nT]);
                    eKthy[nT]=cos(thetaK[nT])*sin(phiK[nT]);
                    eKthz[nT]=sin(thetaK[nT]);
                    OaeKth=Ox[nT]*eKthx[nT]+Oy[nT]*eKthy[nT]+Oz[nT]*eKthz[nT];
                    phir[nT]=acos(OaeKth);
                    NTTO+=1;
                    anGf[nT]=true;
                    //cout<<t<<"    "<<nT<<"    "<<thetar[nT]<<"    "<<phir[nT]<<endl;
                }
                //eKPx[nT]=rKPx[nT]/rKP[nT];
                //eKPy[nT]=rKPy[nT]/rKP[nT];
                //eKPz[nT]=rKPz[nT]/rKP[nT];
                //OaeKR=Ox[nT]*eKPx[nT]+Oy[nT]*eKPy[nT]+Oz[nT]*eKPz[nT];
                //thetar[nT]=acos(OaeKR);
                //thetaK[nT]=acos(eKPz[nT]);
                //phiK[nT]=asin(eKPy[nT]/sin(thetaK[nT]));
                //if(thetaK[nT]>M_PI/2)thetaK[nT]=M_PI-thetaK[nT];
                //phiK[nT]+=M_PI;
                ////cout<<nT<<"    "<<thetaK[nT]<<"    "<<phiK[nT]<<endl;
                //eKthx[nT]=cos(thetaK[nT])*cos(phiK[nT]);
                //eKthy[nT]=cos(thetaK[nT])*sin(phiK[nT]);
                //eKthz[nT]=sin(thetaK[nT]);
                //OaeKth=Ox[nT]*eKthx[nT]+Oy[nT]*eKthy[nT]+Oz[nT]*eKthz[nT];
                //phir[nT]=acos(OaeKth);
                //anGf[nT]=true;

	    }


            for(int nT=0; nT<NT; nT++)
	    {
	        phire[nT]=phir[nT]-phin[t];
                if(phire[nT]<0)phire[nT]+=2*M_PI;
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rho[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
		    if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phir[nT]>=iph*dph)if(phir[nT]<(iph+1)*dph)
		    {
		        rho[ith][iph]+=1;
                    }
                }
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhore[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
		    if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phire[nT]>=iph*dph)if(phire[nT]<(iph+1)*dph)
		    {
		        rhore[ith][iph]+=1;
                    }
                }
            }
            rhodom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhodom+=sin(ith*dth)*rho[ith][iph];
            }
            rhodom=rhodom*M_PI*2*M_PI/((Nth+1)*Nph);
            rhoredom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoredom+=sin(ith*dth)*rhore[ith][iph];
            }
            rhoredom=rhoredom*M_PI*2*M_PI/((Nth+1)*Nph);
            Sint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Sint+=sin(ith*dth)*rho[ith][iph]
                      *(1.0/2.0)*(3*cos(ith*dth)*cos(ith*dth)-1);
            }
            Sint=Sint*M_PI*2*M_PI/((Nth+1)*Nph);
            S[t]=Sint/rhodom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Qint+=sin(ith*dth)*rho[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Q[t]=Qint/rhodom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Qint+=sin(ith*dth)*rhore[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Qre[t]=Qint/rhoredom;



            /*            
            for(int iph=0; iph<Nph2; iph++)rphi[iph]=0;
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
                for(int iph=0; iph<Nph2; iph++)
 		    {
		        if(phir[nT]>=iph*dph2)if(phir[nT]<(iph+1)*dph2)
  		        {
                            rphi[iph]+=1;
                        }
                    }
                }
            }
            for(int iph=0; iph<Nph2; iph++)
	    {
	      //cout<<iph<<"    "<<rphi[iph]<<endl;
            }
            rphimax=0;
            nrphimax=0;
            for(int iph=0; iph<Nph2; iph++)
	    {
                if(rphi[iph]>rphic)if(rphi[iph]>rphimax)
		{
		    rphimax=rphi[iph];
		    nrphimax=iph;
                }
	    }
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
  	            phineu1[nT]=phir[nT]-nrphimax*dph2;
                    if(phineu1[nT]<0)phineu1[nT]+=2*M_PI;
                }
            }
            
            avphi=0;
            navphi=0;
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
                    avphi+=phir[nT];
                    navphi+=1;
                }
	    }
            avphi=avphi/navphi;
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
   	            phineu2[nT]=phir[nT]-avphi;
                    if(phineu2[nT]<0)phineu2[nT]+=2*M_PI;
                }
            }

            for(int iph=0; iph<Nph3; iph++)rphi2[iph]=0;
            for(int nT=0; nT<NT; nT++)
	    {
                if(anGf[nT])for(int iph=0; iph<Nph3; iph++)
		{
		    if(phir[nT]>=iph*dph3)if(phir[nT]<(iph+1)*dph3)
		    {
                        rphi2[iph]+=1;
                    }
                }
            }

            rphimax=0;
            nrphimax=0;
            for(int iph=0; iph<Nph3; iph++)
	    {
                if(rphi[iph]>rphic2)if(rphi[iph]>rphimax)
		{
		    rphimax=rphi2[iph];
		    nrphimax=iph;
                }
	    }
            nrphimaxp=nrphimax+1;
            nrphimaxm=nrphimax-1;
            if(nrphimaxm<0){nrphimaxm+=Nph2;}
            if(nrphimaxp=Nph2){nrphimaxp=0;}
            avphi=0;
            navphi=0;
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
  	            if(phir[nT]>=nrphimax*dph3)if(phir[nT]<(nrphimax+1)*dph3)
		    {
                        avphi+=phir[nT];
                        navphi+=1;
                    }
	            if(phir[nT]>=nrphimaxm*dph3)if(phir[nT]<(nrphimaxm+1)*dph3)
	   	    {
                        avphi+=phir[nT];
                        navphi+=1;
                    }
	            if(phir[nT]>=nrphimaxp*dph3)if(phir[nT]<(nrphimaxp+1)*dph3)
		    {
                        avphi+=phir[nT];
                        navphi+=1;
                    }
                }
	    }
            avphi=avphi/navphi;
            //cout<<avphi<<endl;
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
  	            phineu3[nT]=phir[nT]-avphi;
                    if(phineu3[nT]<0)phineu3[nT]+=2*M_PI;
                }
            }
            for(int nT=0; nT<NT; nT++)
	    {
	        if(anGf[nT])
		{
	            phineu4[nT]=phir[nT]-Cphi4;
  	            //phineu4[nT]=phi[nT]-avphi;
                    if(phineu4[nT]<0)phineu4[nT]+=2*M_PI;
                    if(phineu4[nT]>=(2*M_PI))phineu4[nT]=phineu4[nT]-2*M_PI;
                }
            }
            Qmax[t]=0;
	    dphiQmax[t]=0;
            for(int iph=0; iph<Nph4; iph++)
	    {
	        double dphi=dph4*iph;
                int ntest=0;
                //for(int nT=0; nT<NT; nT++)
   	        //{
                    //phi[nT]=0;
                    //theta[nT]=M_PI/2;
                //}
                for(int nT=0; nT<NT; nT++)
   	        {
  	            if(anGf[nT])
   		    {
  	                phineu5[nT]=phir[nT]-dphi;
                        if(phineu5[nT]<0)phineu5[nT]+=2*M_PI;
                        if(phineu5[nT]>=(2*M_PI))phineu5[nT]=phir[nT]-2*M_PI;
                    }
                }
                for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
  	        {
	            rhoneu5[ith][iph]=0;
                    for(int nT=0; nT<NT; nT++)
  		    {
  	                if(anGf[nT])
   		        {
		            if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phineu5[nT]>=iph*dph)if(phineu5[nT]<(iph+1)*dph)
 		            {
		                rhoneu5[ith][iph]+=1;
                                 ntest+=1;
                            }
                        }
		    }
                }
                //rho[Nth/2][0]=1;
                rhoneu5dom=0;
 	        for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
 	        {
	            rhoneu5dom+=sin(ith*dth)*rhoneu5[ith][iph];
                    //cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
                }
                rhoneu5dom=rhoneu5dom*M_PI*2*M_PI/((Nth+1)*Nph);
                Qint=0;
   	        for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	        {
	            Qint+=sin(ith*dth)*rhoneu5[ith][iph]
                          *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
                }
                Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
                Qn5[t]=Qint/rhoneu5dom;
                //cout<<iph<<"    "<<dphi<<"    "<<Qn5[t]<<"    "<<ntest<<endl;
                //cout<<iph<<"    "<<dphi<<"    "<<Qn5[t]<<"    "<<rhoneu5dom<<endl;
                if(Qn5[t]>Qmax[t])
                {
                    Qmax[t]=Qn5[t];
                    dphiQmax[t]=dphi;
                }

	    }





            




            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rho[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
                    if(anGf[nT])
                    {
  	 	        if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phir[nT]>=iph*dph)if(phir[nT]<(iph+1)*dph)
	  	        {
	  	            rho[ith][iph]+=1;
                        }
                    }
                }
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu1[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
                    if(anGf[nT])
                    {
   		        if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phineu1[nT]>=iph*dph)if(phineu1[nT]<(iph+1)*dph)
		        {
		            rhoneu1[ith][iph]+=1;
                        }
                    }
                }
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu2[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
		    if(anGf[nT])if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phineu2[nT]>=iph*dph)if(phineu2[nT]<(iph+1)*dph)
		    {
		        rhoneu2[ith][iph]+=1;
                    }
                }
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu3[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
		    if(anGf[nT])if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phineu3[nT]>=iph*dph)if(phineu3[nT]<(iph+1)*dph)
		    {
		        rhoneu3[ith][iph]+=1;
                    }
                }
            }
            for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu4[ith][iph]=0;
                for(int nT=0; nT<NT; nT++)
		{
		    if(anGf[nT])if(thetar[nT]>=ith*dth)if(thetar[nT]<(ith+1)*dth)if(phineu4[nT]>=iph*dph)if(phineu4[nT]<(iph+1)*dph)
		    {
		        rhoneu4[ith][iph]+=1;
                    }
                }
            }
            rhodom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhodom+=sin(ith*dth)*rho[ith][iph];
                //if(rho[ith][iph]!=0)cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            rhodom=rhodom*M_PI*2*M_PI/((Nth+1)*Nph);
            //cout<<rhodom<<endl;
            rhoneu1dom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu1dom+=sin(ith*dth)*rhoneu1[ith][iph];
                //cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            rhoneu1dom=rhoneu1dom*M_PI*2*M_PI/((Nth+1)*Nph);
            rhoneu2dom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu2dom+=sin(ith*dth)*rhoneu2[ith][iph];
                //cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            rhoneu2dom=rhoneu2dom*M_PI*2*M_PI/((Nth+1)*Nph);
            rhoneu3dom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu3dom+=sin(ith*dth)*rhoneu3[ith][iph];
                //cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            rhoneu3dom=rhoneu3dom*M_PI*2*M_PI/((Nth+1)*Nph);
            rhoneu4dom=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        rhoneu4dom+=sin(ith*dth)*rhoneu4[ith][iph];
                //cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            rhoneu4dom=rhoneu4dom*M_PI*2*M_PI/((Nth+1)*Nph);
            //cout<<rhodom<<endl;
            Sint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Sint+=sin(ith*dth)*rho[ith][iph]
                      *(1.0/2.0)*(3*cos(ith*dth)*cos(ith*dth)-1);
            }
            Sint=Sint*M_PI*2*M_PI/((Nth+1)*Nph);
            S[t]=Sint/rhodom;
            //cout<<t<<"    "<<S[t]<<"    "<<Sint<<"    "<<rhodom<<endl;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        //Qint+=sin(ith*dth)*rho[ith][iph]
                      //*(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(iph*dph));
	        //Qint+=(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(iph*dph));
                //Qint+=sin(ith*dth)*rho[ith][iph];
	        Qint+=sin(ith*dth)*rho[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
	        //Qint+=sin(ith*dth)*rho[ith][iph]
                      //*(3.0/2.0)*(sin(ith*dth)*sin(ith*dth));
	        //Qint+=sin(ith*dth)*rho[ith][iph]
                      //*(3.0/2.0)*(cos(2*iph*dph));
                //if(t==10000)if(ith==2)cout<<iph<<"    "<<cos(2*iph*dph)<<endl;
                //if(t==10000)if(ith==2)cout<<iph<<"    "<<rho[ith][iph]<<endl;
                //if(t==10000)if(rho[ith][iph]!=0)cout<<ith<<"    "<<iph<<"    "<<rho[ith][iph]<<endl;
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Q[t]=Qint/rhodom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++) 
	    {
	        Qint+=sin(ith*dth)*rhoneu1[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Qn1[t]=Qint/rhoneu1dom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Qint+=sin(ith*dth)*rhoneu2[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Qn2[t]=Qint/rhoneu2dom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Qint+=sin(ith*dth)*rhoneu3[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Qn3[t]=Qint/rhoneu3dom;
            Qint=0;
	    for(int ith=0; ith<(Nth+1); ith++)for(int iph=0; iph<Nph; iph++)
	    {
	        Qint+=sin(ith*dth)*rhoneu4[ith][iph]
                      *(3.0/2.0)*(sin(ith*dth)*sin(ith*dth)*cos(2*iph*dph));
            }
            Qint=Qint*M_PI*2*M_PI/((Nth+1)*Nph);
            Qn4[t]=Qint/rhoneu4dom;
            //cout<<t<<"    "<<Q[t]<<"    "<<Qint<<"    "<<rhodom<<endl;
            
            //for(int nT=0; nT<NT; nT++)
	    //{

            //}
            //cout<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qn1[t]<<"    "<<Qn2[t]<<"    "<<Qn3[t]<<"    "<<Qn4[t]<<"    "<<Qmax[t]<<"    "<<dphiQmax[t]<<endl;
            cout<<t<<"    "<<S[t]<<"    "<<Qmax[t]<<"    "<<NTTO<<"    "<<DDef<<endl;
            outOP<<t<<"    "<<S[t]<<"    "<<Qmax[t]<<"    "<<NTTO<<"    "<<DDef<<endl;
            */
            //cout<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<endl;
            //outOP<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<endl;
            cout<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<NTTO<<"    "<<DDef<<endl;
            outOP<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<NTTO<<"    "<<DDef<<endl;
            







 
        }

    }


    for(int t=Nta; t<(Nt+1); t++)
    {
        //if(t%dt==0)cout<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qn1[t]<<"    "<<Qn2[t]<<"    "<<Qn3[t]<<"    "<<Qn4[t]<<"    "<<Qmax[t]<<"    "<<dphiQmax[t]<<endl;
    }


    return 0;
} 
