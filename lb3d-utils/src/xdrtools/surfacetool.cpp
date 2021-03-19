// code by Fabian Doerfler, MPI-MF Stuttgart, Sept. 2008
//
//     creates xdr output in order to specify surface patterns to be used within 
//     LBE simulations of wetting and fluid dynamics on chemically patterned substrates
//
//     for further information and instructions read <surfacetool.readme>
//

//#include<cstdlib>
//#include<cstring>
#include <iostream>
#include<iomanip>
#include<fstream>
#include<strstream>

#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <rpc/types.h>
#include <math.h>

#define max 121

using namespace std;

int main(int argc, char* argv[])
{
  if (argc!=2) {
    cout << "\nsyntax:\n\t" << argv[0] << " <inputfile>" << "\n \n";
    exit(-1);
  }

  int nx, ny, nz;
  int thick, hwidth;
  int bigrad, smallrad, radquad;   // bigrad: big radius; smallrad: small radius
  int midx, midy, midz;
  float vac=0;
  float color[3]={0,0,vac};   // color[0]: ring; color[1]: vicinity

  XDR xdrs;
  FILE* fp;

  ifstream ifs;
  ofstream ofs;
  
  char txtfile[max];
  char xdrfile[max];
  char* infile=argv[1];

  ostrstream strtxt(txtfile,120);
  ostrstream strxdr(xdrfile,120);

  strtxt << infile << ".txt" << '\0';
  strxdr << infile << ".xdr" << '\0';
  
  ifs.open(infile, ios::in);
  if(!infile) {
    cout << "\n" << "could not open " << infile;
    exit(-1);
  }
  ofs.open(strtxt.str(), ios::out);
  fp=fopen(xdrfile, "w");
  
  // get parameters from input file:
  ifs >> nx >> ny >> nz >> bigrad >> smallrad >> midx >> midy >> midz >> color[0] >> color[1] >> thick;
  cout << "\n";
  cout << "lattice: " << nx << ", " << ny << ", " << nz << "\n";
  cout << "radii: " << bigrad << ", " << smallrad << "\n";
  cout << "center: " << midx << ", " << midy << ", " << midz << "\n";
  cout << "rockcolors: " << color[0] << ", " << color[1] <<'\n';
  cout << "thickness: " << thick << "\n \n";

  // allocate memory for 3-d array //////////////////////////////
  float*** p;  
  p=new float**[nx];
  if(p==0) {
    cout << "memory allocation for <float** p[nx]> failed";
    exit(1);
  }
  p[0]=new float*[nx*ny];
  if(p[0]==0) {
    cout << "memory allocation for <float* p[nx][ny]> failed";
    exit(1);
  }
  p[0][0]=new float[nx*ny*nz];
  if(p[0][0]==0) {
    cout << "memory allocation for <float p[nx ][ny][nz]> failed";
    exit(1);
  }

  // set the pointers 
  for(int j=1;j<ny;j++){       
    p[0][j]=p[0][j-1]+nz;      
  }                            
   for(int i=1;i<nx;i++){
    p[i]=p[i-1]+ny;            
    p[i][0]=p[i-1][0]+ny*nz;   
    for(int j=1;j<ny;j++){
      p[i][j]=p[i][j-1]+nz;    
    
    }
  }

  // set everything to zero                               
  for(int k=0; k<nz; k++) {    
    for(int j=0; j<ny; j++) {
      for(int i=0; i<nx; i++) {
	p[i][j][k]=0;
      }
    }
  }

  // create the rockcolor pattern
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      for(int k=0; k<nz; k++) {
	if(j<thick || j>=ny-thick) { // bottom and top substrate layers with y-normal
	  radquad=(i-midx)*(i-midx)+(k-midz)*(k-midz);
	  hwidth=(bigrad-smallrad)/2;
	  if(radquad>=bigrad*bigrad && abs(midx-i)<=hwidth) { // straight channel
	    p[i][j][k]=color[0];
	  }
	  else{
	    p[i][j][k]=color[1];
	  }
	  if(radquad<bigrad*bigrad && radquad>=smallrad*smallrad){ // ring-channel
	    p[i][j][k]=color[0];
	  }
	}
	else {
	  p[i][j][k]=color[2];
	}
      }
    }
  }

  // create xdr output and text output
  xdrstdio_create(&xdrs, fp, XDR_ENCODE);
  for(int k=0;k<nz;++k){
    for(int j=0;j<ny;++j){
      for(int i=0;i<nx;++i){
	xdr_float(&xdrs,&p[i][j][k]);                                   // xdr output
	ofs << i << " " << j << " " << k << " " << p[i][j][k] << "\n";  // text output
      }
    }
  }
  
  xdr_destroy(&xdrs);
  fclose(fp);
  delete p[0][0];
  delete p[0];
  delete p;
  
}



