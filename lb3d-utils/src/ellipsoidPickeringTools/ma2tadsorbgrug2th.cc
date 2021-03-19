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
    const int Nt=330000;
    const int NT=1600;
    const int dt=10000;
    string datzeile;
    int Ni=0;
    bool inzahl=false;
    bool inexp=false;
    int nz=0;
    string Zahlstring[Nz][NT+1];
    string Zahlstringcol[Nz];
    char Zahlchar;
    string datei1st="md-cfg_absorb01_t00000000-0705201734.asc";
    //const int TNr=15;
    //const int XNr=0;
    //const int YNr=1;
    //const int ZNr=2;
    //const int OXNr=9;
    //const int OYNr=10;
    //const int OZNr=11;
    const int TNr=16;
    const int XNr=1;
    const int YNr=2;
    const int ZNr=3;
    const int OXNr=10;
    const int OYNr=11;
    const int OZNr=12;
    const string lb_gr_out_file="lauf08";
    const string lbID="";
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
    const double Rp=8;
    const double Ro=4;
    double dGf[Nt+1];
    int Ndgf=0;
    double md=Rp+4;
    const int Gx=512;
    const int Gy=80;
    const int Gz=512;
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
    double theta[NT+1];
    double phi[NT+1];
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
    double Qxx=0;
    double Qxy=0;
    double Qyx=0;
    double Qyy=0;
    double L=0;
    double vex=0;
    double vey=0;
    double vba=0;
    double phin[Nt+1];
    int Ngr[NT+1][Gx/2+2];
    double gr[NT+1][Gx/2+2];
    double Gr[Gx/2+2];
    double Grn[Gx/2+2];
    double g2th[NT+1][Gx/2+2];
    double G2th[Gx/2+2];
    double Dr=0;
    int iDr=0;
    double dNgr[Gx/2+2];
    double Grm[Gx/4+2];
    double Grmn[Gx/4+2];
    double Grm2[Gx/8+2];
    double Grm2n[Gx/8+2];
    double Grm3[Gx/2+2];
    double Grm3n[Gx/2+2];
    double Grm4[Gx/8+2];
    double Grm4n[Gx/8+2];
    double Grm5[Gx/2+2];
    double Grm5n[Gx/2+2];
    int NdGrb=0;
    double dGrb=0;
    double G2thm[Gx/4+2];
    double G2thm2[Gx/2+2];
    const double dphic=M_PI/10;
    int Ngrp[NT+1][Gx/2+2];
    double grp[NT+1][Gx/2+2];
    double Grp[Gx/2+2];
    double Grpn[Gx/2+2];
    int NdGrpb=0;
    double dGrpb=0;
    double dNgrp[Gx/2+2];
    string stringGra="";
    string stringGrm5a="";
    string stringG2tha="";
    string stringG2thm2a="";


    filepart1=filepart1a+lb_gr_out_file+filepart1b;
    filepart2=lbID+fileend;
    filecolpart1=filecolpart1a+lb_gr_out_file+filepart1b;
    fstream outOP;
    outOP.open("OPSQm", ios::out);
    fstream outGR;
    outGR.open("GR", ios::out);
    fstream outGr;
    outGr.open("Gr", ios::out);
    fstream outGrm;
    outGrm.open("Grm", ios::out);
    fstream outGrm2;
    outGrm2.open("Grm2", ios::out);
    fstream outGrm3;
    outGrm3.open("Grm3", ios::out);
    fstream outGrm4;
    outGrm4.open("Grm4", ios::out);
    fstream outGrm5;
    outGrm5.open("Grm5", ios::out);
    fstream outG2th;
    outG2th.open("G2th", ios::out);
    fstream outG2thm;
    outG2thm.open("G2thm", ios::out);
    fstream outG2thm2;
    outG2thm2.open("G2thm2", ios::out);
    fstream outGrp;
    outGrp.open("Grp", ios::out);
    fstream outTest;
    outTest.open("testOP", ios::out);
    fstream outStatus;
    outStatus.open("StatusOP", ios::out);
    fstream outLongStatus;
    outLongStatus.open("LongStatusOP", ios::out);
    for(int t=Nta; t<(Nt+1); t++)
    {
        if(t%dt==0)
	{
            cout<<"t="<<t<<"    t/dt="<<t/dt<<endl;
            outStatus<<"t="<<t<<"    t/dt="<<t/dt<<endl;
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
            else if(t<1000000)
	    {
                Nl="00";
            }
            if (sst << t)
            {
                std::string ssst(sst.str());
                st=ssst;
            }
            stringGra="Gra_"+Nl+st;
            stringGrm5a="Grm5a_"+Nl+st;
            stringG2tha="G2tha_"+Nl+st;
            stringG2thm2a="G2thm2a_"+Nl+st;
            fstream outGra;
            outGra.open(stringGra.c_str() , ios::out);
            fstream outGrm5a;
            outGrm5a.open(stringGrm5a.c_str() , ios::out);
            fstream outG2tha;
            outG2tha.open(stringG2tha.c_str() , ios::out);
            fstream outG2thm2a;
            outG2thm2a.open(stringG2thm2a.c_str() , ios::out);
            outfilestr=filepart1+Nl+st+filepart2;
            outcolfilestr=filecolpart1+Nl+st+filepart2;
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
                    if(Ne<-16)X[nt]=0;//cout<<nt<<"    "<<X[nt]<<endl;
                    





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
                    theta[nt]=acos(Oy[nt]);
                    phi[nt]=asin(Ox[nt]/sin(theta[nt]));
                    if(Ox[nt]/sin(theta[nt])<-1)phi[nt]=asin(-1);
                    if(Ox[nt]/sin(theta[nt])>1)phi[nt]=asin(1);
                    if(theta[nt]>M_PI/2)theta[nt]=M_PI-theta[nt];
                    phi[nt]+=M_PI;
                }
                
            }



            cout<<"Dateien ausgelesen"<<endl;
            outStatus<<"Dateien ausgelesen"<<endl;


            for(int nT=0; nT<NT; nT++) 
            {
                for(int i=0; i<Gx/2; i++)
		{
                    Ngr[nT][i]=0;
                    gr[nT][i]=0;
                    Gr[i]=0;
                    dNgr[i]=0;
                    g2th[nT][i]=0;
                    G2th[i]=0;
                }
            }
            for(int nT=0; nT<NT; nT++)
            {
                for(int nt=0; nt<NT; nt++)if(nt!=nT)
		{
                    Dx=X[nt]-X[nT];
                    if(Dx<0)Dx=-Dx;
                    if(Dx>Gx/2)Dx=Gx-Dx;
                    Dz=Z[nt]-Z[nT];
                    if(Dz<0)Dz=-Dz;
                    if(Dz>Gz/2)Dz=Gz-Dz;
                    Dr=sqrt(Dx*Dx+Dz*Dz);
                    iDr=int(Dr+0.5);
                    //Ngr[nT][iDr]+=1;
                    //if(nT==0)cout<<nt<<"    "<<X[nt]<<"    "<<X[nT]<<"    "<<Dx<<"    "<<Dr<<"    "<<iDr<<endl;
                    //if(nT==0)cout<<nt<<"    "<<iDr<<"    "<<Ngr[nT][iDr]<<endl;
                    //g2thsum[nT][iDr]+=cos(2*(phi[nT]-phi[nt]));
                    //g2th[nT][iDr]+=cos(2*(phi[nT]-phi[nt]));
                    //if(nT==0)if(iDr<20)cout<<nt<<"    "<<iDr<<"    "<<Dr<<endl;
                    //if(iDr<6)cout<<nt<<"    "<<nT<<"    "<<iDr<<"    "<<Dr<<endl;
                    //if(Dr<6)cout<<nt<<"    "<<nT<<"    "<<iDr<<"    "<<Dr<<endl;
                    //if(Dr>Gx/2-2)cout<<nt<<"    "<<nT<<"    "<<Dr<<"    "<<iDr<<"    "<<Gx/2<<endl;
                    //if(Dr>Gx/2-2)cout<<nt<<"    "<<nT<<"    "<<Dr<<"    "<<Dx<<"    "<<Dz<<endl;
                    if(Dr<Gx/2)
		    {
                        double dphi=0;
                        Ngr[nT][iDr]+=1;
                        //g2th[nT][iDr]+=cos(2*(phi[nT]-phi[nt]));
                        dphi=phi[nT]-phi[nt];
                        if(dphi<0)dphi=-dphi;
                        if(dphi>M_PI/2)dphi=M_PI-dphi;
                        g2th[nT][iDr]+=cos(2*dphi);
                        if(dphi<dphic)Ngrp[nT][iDr]+=1;
                    }
                }
                //if(gr[0][1]>0)cout<<nT<<"    "<<Ngr[nT][1]<<"    "<<gr[nT][1]<<"    "<<Gr[1]<<endl;
                for(int i=1; i<Gx/2; i++)
		{
                    gr[nT][i]=Ngr[nT][i]/(2*M_PI*i);
                    //if(nT==0)cout<<i<<"    "<<Ngr[nT][i]<<endl;
                    Gr[i]+=gr[nT][i];
                    //g2th[nT][iDr]=g2thsum[nT][iDr]/Ngr[nT][iDr];
                    //g2th[nT][i]=g2th[nT][i]/Ngr[nT][i];
                    if(Ngr[nT][i]>0)g2th[nT][i]=g2th[nT][i]/Ngr[nT][i];
                    dNgr[i]+=Ngr[nT][i];
                    if(Ngr[nT][i]>0)G2th[i]+=g2th[nT][i];
                    grp[nT][i]=Ngrp[nT][i]/(2*M_PI*i);
                    Grp[i]+=grp[nT][i];
                }
            }
            //cout<<Ngr[0][1]<<"    "<<gr[0][1]<<"    "<<Gr[1]<<endl;
            NdGrb=0;
            NdGrpb=0;
            dGrpb=0;
            for(int i=1; i<Gx/2; i++)
	    {
                Gr[i]=Gr[i]/NT;
                dNgr[i]=dNgr[i]/NT;
                G2th[i]=G2th[i]/NT;
                Grp[i]=Grp[i]/NT;
                if(i>Gx/4)
		{
                    dGrb+=Gr[i];
                    NdGrb+=1;
                    dGrpb+=Grp[i];
                    NdGrpb+=1;
                }
            }
            dGrb=dGrb/NdGrb;
            dGrpb=dGrpb/NdGrpb;
            Grm[0]=0;
            for(int i=1; i<Gx/4; i++)
	    {
                Grm[i]=(Gr[2*i-1]+2*Gr[2*i]+Gr[2*i+1])/4;
            }
            Grm2[0]=0;
            for(int i=1; i<Gx/8; i++)
	    {
                Grm2[i]=(Gr[4*i-2]+2*Gr[4*i-1]+2*Gr[4*i]+2*Gr[4*i+1]+Gr[4*i+2])/8;
            }
            Grm3[0]=0;
            for(int i=1; i<Gx/2; i++)
	    {
                Grm3[i]=(Gr[i-2]+2*Gr[i-1]+4*Gr[i]+2*Gr[i+1]+Gr[i+2])/10;
            }
            Grm4[0]=0;
            for(int i=1; i<Gx/8; i++)
	    {
                Grm4[i]=(Gr[4*i-3]+2*Gr[4*i-2]+3*Gr[4*i-1]+4*Gr[4*i]+3*Gr[4*i+1]+2*Gr[4*i+2]+Gr[4*i+3])/16;
            }
            Grm5[0]=0;
            for(int i=1; i<Gx/2; i++)
	    {
                Grm5[i]=(Gr[i-3]+2*Gr[i-2]+3*Gr[i-1]+4*Gr[i]+3*Gr[i+1]+2*Gr[i+2]+Gr[i+3])/16;
                if(Gr[i]==0)Grm5[i]=0;
            }
            G2thm[0]=0;
            for(int i=1; i<Gx/4; i++)
	    {
                G2thm[i]=(G2th[2*i-1]+2*G2th[2*i]+G2th[2*i+1])/4;
            }
            G2thm2[0]=0;
            for(int i=1; i<Gx/2; i++)
	    {
                G2thm2[i]=(G2th[i-3]+2*G2th[i-2]+3*G2th[i-1]+4*G2th[i]+5*G2th[i+1]+2*G2th[i+2]+G2th[i+3])/16;
                if(G2th[i]==0)G2thm2[i]=0;
            }
            for(int i=1; i<Gx/2; i++)
	    {
                //Gr[i]=Gr[i]/NT;
                //dNgr[i]=dNgr[i]/NT;
                //G2th[i]=G2th[i]/NT;
                Grn[i]=Gr[i]/dGrb;
                Grpn[i]=Grp[i]/dGrpb;
                if(i%2==0)Grmn[i/2]=Grm[i/2]/dGrb;
                if(i%4==0)Grm2n[i/4]=Grm2[i/4]/dGrb;
                Grm3n[i]=Grm3[i]/dGrb;
                if(i%4==0)Grm4n[i/4]=Grm4[i/4]/dGrb;
                Grm5n[i]=Grm5[i]/dGrb;
                //if(i%2==0)cout<<i<<endl;
                //if(i%2==0)Grm[i/2]=(Gr[i-1]+2*Gr[i]+Gr[i+1])/4;
                cout<<t<<"    "<<i<<"    "<<Gr[i]<<"    "<<Grn[i]<<"    "<<gr[0][i]<<"    "<<gr[1][i]<<"    "<<Gr[i]*NT<<"    "<<dNgr[i]<<endl;
                //if(i%2==0)cout<<t<<"    "<<i<<"    "<<Grm[i/2]<<endl;
                outGr<<t<<"    "<<i<<"    "<<Gr[i]<<"    "<<Grn[i]<<"    "<<gr[0][i]<<"    "<<gr[1][i]<<"    "<<Gr[i]*NT<<"    "<<dNgr[i]<<endl;
                outGra<<t<<"    "<<i<<"    "<<Gr[i]<<"    "<<Grn[i]<<"    "<<gr[0][i]<<"    "<<gr[1][i]<<"    "<<Gr[i]*NT<<"    "<<dNgr[i]<<endl;
                //if(i%2==0)outGrm<<t<<"    "<<i<<"    "<<Grm[i/2]<<endl;
                if(i%2==0)outGrm<<t<<"    "<<i<<"    "<<Grm[i/2]<<"    "<<Grmn[i/2]<<endl;
                if(i%4==0)outGrm2<<t<<"    "<<i<<"    "<<Grm2[i/4]<<"    "<<Grm2n[i/4]<<endl;
                outGrm3<<t<<"    "<<i<<"    "<<Grm3[i]<<"    "<<Grm3n[i]<<endl;
                if(i%4==0)outGrm4<<t<<"    "<<i<<"    "<<Grm4[i/4]<<"    "<<Grm4n[i/4]<<endl;
                outGrm5<<t<<"    "<<i<<"    "<<Grm5[i]<<"    "<<Grm5n[i]<<endl;
                outGrm5a<<t<<"    "<<i<<"    "<<Grm5[i]<<"    "<<Grm5n[i]<<endl;
                //cout<<"test"<<endl;
                outG2th<<t<<"    "<<i<<"    "<<G2th[i]<<"    "<<g2th[0][i]<<"    "<<g2th[1][i]<<endl;
                outG2tha<<t<<"    "<<i<<"    "<<G2th[i]<<"    "<<g2th[0][i]<<"    "<<g2th[1][i]<<endl;
                if(i%2==0)outG2thm<<t<<"    "<<i<<"    "<<G2thm[i/2]<<endl;
                outG2thm2<<t<<"    "<<i<<"    "<<G2thm2[i]<<endl;
                outG2thm2a<<t<<"    "<<i<<"    "<<G2thm2[i]<<endl;
                //cout<<"test"<<endl;
                outGrp<<t<<"    "<<i<<"    "<<Grp[i]<<"    "<<Grpn[i]<<"    "<<grp[0][i]<<"    "<<grp[1][i]<<"    "<<Grp[i]*NT<<"    "<<dNgrp[i]<<endl;
            }





            //cout<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<endl;
            //outOP<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qre[t]<<"    "<<phin[t]<<endl;






 
        }

    }

    //fstream outOP;
    //outOP.open("OPSQ", ios::out);
    for(int t=Nta; t<(Nt+1); t++)
    {
        //if(t%dt==0)outOP<<t<<"    "<<S[t]<<"    "<<Q[t]<<"    "<<Qn1[t]<<"    "<<Qn2[t]<<"    "<<Qn3[t]<<"    "<<Qn4[t]<<"    "<<Qmax[t]<<"    "<<dphiQmax[t]<<endl;
    }


    return 0;
} 
