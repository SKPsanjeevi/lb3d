/* 
 *  PROGRAM BY Frank Raischel, UNI STUTTGART
 *  developed April 2007
 *
 *  reads 1d rock files 
 *  writes xdr,vtk and ascii file as output file and creates wall around the sample 
 *  don't allocate array, set -DDEBUG for text output
 *  working version
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <math.h>
#include <string.h>



#define MAXSTRING 160



int main (int argc, char *argv[])
{
	FILE	    *f0,*f1,*f2, *VTK, *DBG;
	XDR	    xdrs;
	char        *fname;
	char        extxdr[] = ".xdr",ext[] = ".txt", extvtk[]=".vtk";
	char        xdrfile[MAXSTRING],txtfile[MAXSTRING], vtkfile[MAXSTRING], debugfile[MAXSTRING];
	char        command[MAXSTRING], syscmd[MAXSTRING];
	char        unziponefile[MAXSTRING];
	char        usage[MAXSTRING];
	char        undecided[5]="";
	int         nx,ny,nz,nxx,nyy,nzz,chanx,chany,chanz, mh,rdn,hgt,seed;
	float       ***pos;
	int       *tx1,*tx2,*ty1,*ty2,*tz1,*tz2 ;/*3d channels */
	int       *width; /* 'Width' of the channel*/
	float	   *coulor;
	float	  mhgt;
	int     i,j,k,tmp, countsites;
	char line[180]="";
	float im, imc;
	int cx,cy,cz;
	int zipflag=0;
	int u2m=1,u2p=0;//default
	int encase=0;

	sprintf(usage, "Error: Syntax: %s <filename> <(3xint) input rock size> <(3xint) output size > <(char) u2m/u2p\n",argv[0]);
	if (argc<8||argc>10) {
		fprintf(stderr,usage);
		fprintf(stderr, "Error: Less than 7 or more than 9 arguments given!\n");
		return -1;
	}

	fname=argv[1];
	if(NULL!=(strstr(fname, ".gz"))){
		zipflag=1;
	}
	
	sscanf(argv[2], "%d", &nxx);// input sizes
	sscanf(argv[3], "%d", &nyy);//"" 
	sscanf(argv[4], "%d", &nzz);//""
	sscanf(argv[5], "%d", &nx);//output sizes
	sscanf(argv[6], "%d", &ny);//""
	sscanf(argv[7], "%d", &nz);//""
	sprintf(undecided, "u2m");
	if(argc==9){
		if(NULL!=strstr(argv[8], "u2m")){//convert undecided sites (-1) to matrix
			u2m=1;u2p=0;
			sprintf(undecided, "u2m");
			fprintf(stdout, "Option 'u2m' is set.\n");
		}
		else if(NULL!=strstr(argv[8], "u2p")){//convert undecided sites (-1) to pores
			u2p=1;u2m=0;
			sprintf(undecided, "u2p");
			fprintf(stdout, "Option 'u2p' is set.\n");
		}
		else{
			fprintf(stderr,usage);
			fprintf(stderr, "Give option 'u2p' or 'u2m'\n");
			return -1;
		}
	}
		

	//set channel widths in x,y and  fluid zone in z direction: 	


 	if(nx>nxx){
		chanx=(nx-nxx)/2; 
	} 
 	else{
		chanx=1;
	}//	nxx=nx-2;	} 
 	if(ny>nyy){
	 	chany=(ny-nyy)/2; 
	} 
 	else{
		chany=1;
	}//	nyy=ny-2;	} 
 	if(nz>nzz){
	 	chanz=(nz-nzz)/2; 
	} 
 	else{
		chanz=8;
	}//	nzz=nz-16;	} 

	if(nx<=nxx && ny<=nyy && nz>nzz){
		encase=2;
	}
	else if(nx<=nxx && ny<=nyy && nz<=nzz){
		encase=1;
	}
	else{
		encase=0;
	}
	printf("encase=%d\n", encase);
	printf("channels=%d,%d,%d\n", chanx,chany,chanz);



	printf("Output system size: %d %d %d  \n",nx, ny, nz  );
	//srand(seed);
	
	//sprintf(txtfile,"%s%s",fname,ext);

	/* Allocate tensor (i.e. 3d array) and set pointers */

	
/* 	pos=(float ***)malloc(nx*sizeof(float **)); */
/* 	if(!pos){fprintf(stderr,"Error: (1) can't allocate memory \n"); */
/* 		exit(1);} */
/* 	pos[0]=(float **)malloc(nx*ny*sizeof(float *)); */
/* 	if(!pos[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(-1);} */
/* 	pos[0][0]=(float *)calloc(nx*ny*nz,sizeof(float)); */
/* 	if(!pos[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(-1);} */
/* 	/\* printf("pt 1 \n"); *\/ */
	
/*  	for(j=1;j<ny;j++){ */
/* 		pos[0][j]=pos[0][j-1]+nz; */
/* 	}  */
	
/*  	for(i=1;i<nx;i++){  */
/*  		pos[i]=pos[i-1]+ny;  */
/*  		pos[i][0]=pos[i-1][0]+ny*nz;  */
/*  		for(j=1;j<ny;j++){  */
/*  			pos[i][j]=pos[i][j-1]+nz;  */
/*  		}  */
/*  	}   */
		// 	printf("pt 2 \n ");   

	
	sprintf(txtfile,"%s_%i_%i_%i_r7_%s%s",fname,nx,ny,nz,undecided,ext);
	sprintf(xdrfile,"%s_%i_%i_%i_r7_%s%s",fname,nx,ny,nz,undecided,extxdr);
	sprintf(vtkfile,"%s_%i_%i_%i_r7_%s%s",fname,nx,ny,nz,undecided,extvtk);
	sprintf(debugfile,"%s_%i_%i_%i_r7_%s_debug%s",fname,nx,ny,nz,undecided,ext);

	/* Create xdrfile, stop if it already exists/overwrite? */
	sprintf(command,"touch %s ",xdrfile);
	system(command);
	if(NULL==(f2=fopen(xdrfile,"w"))){
		perror("Unable to open output file");
		return -1;
	}        
	printf("Open XDR file %s \n",xdrfile);
	if(NULL==(VTK=fopen(vtkfile,"w"))){
		perror("Unable to open output file");
		return -1;
	}        
	printf("Open VTK file %s \n", vtkfile);
	
#ifdef DEBUG 
	printf(" Open ASCII file %s \n",txtfile );
	f1 = fopen(txtfile,"w");
	printf("Open DEBUG file %s \n",debugfile);
	if(NULL==(DBG=fopen(debugfile,"w"))){
		fprintf(stderr, "Unable to open output file %s \n",debugfile );
		return -1;
	}
#endif 
		
/* 	/\* In the beginning set all points to 0   *\/  */

/*  for(cz=0;cz<nz;cz++){ */
/*            for(cy=0;cy<ny;cy++){ */
/*             for(cx=0;cx<nx;cx++){ */
/* 							if(cx<chanx||cx>=nxx+chanx||cy<chany||cy>=nyy+chany ){ */
/* 								/\*only the channel walls are nonzero:*\/ */
/* 								pos[cx][cy][cz] = 5.0;  */
/* 							} */
/* 							else{ */
/* 								pos[cx][cy][cz] = 0.0;  */
/* 							} */
/* 						} */
/* 					 } */
/*  } */
	/* printf("nach der 0\n" ); */

                                       
		 if(!zipflag){
			 if( NULL==(f0 = fopen(fname,"r")) ){
				 fprintf(stderr, "Cannot open file %s !\n", fname);
				 exit(-1);
			 } 
		 }
		 else{
			 memset(unziponefile, 0, 160);
			 sprintf(unziponefile, "gunzip -c %s", fname);
			 if( NULL==(f0=popen(unziponefile, "r")) ){
				 fprintf(stderr, "Cannot open pipe %s !\n", unziponefile);
				 exit(-1);
			 }
		 }
		 
		 countsites=0;
		 xdrstdio_create(&xdrs,f2,XDR_ENCODE); 

		 if(encase==0){
			 for(cz=0;cz<nz;cz++){ 
				 for(cy=0;cy<ny;cy++){ 
					 for(cx=0;cx<nx;cx++){ 
						 if( (cx>=chanx)&&(cx<nx-chanx)&&(cy>=chany)&&(cy<ny-chany)&&(cz>=chanz)&&(cz<nz-chanz) ){						 						 
							 fgets(line, 180, f0);
							 countsites++;
							 sscanf(line, "%f", &im);
							 if(im==0){//pore
								 im=0.;
							 }
							 else if(im==125){//matrix
								 im=5.;
							 }
							 else if(im>0 && im<125 ){//undecided, is u2m or p2m set?
								 if(u2m){
									 im=5.;
								 }
								 else if(u2p){
									 im=0.;
							 }
								 else{}
							 }
							 else{
								 fprintf(stderr, "Error: Undefined value %f at position %i\n", im, countsites);
								 exit(-1);
							 }
						 }
						 //and now,  channels and fluid:
						 //						 imc=cx;cx=cz;
						 if(  (cx<chanx)||(cx>=nx-chanx)||(cy<chany)||(cy>=ny-chany) ){
							 im=5.;
						 }
						 else if( (cz<chanz)||(cz>=nz-chanz)  ){
							 im=0.;
						 }
						 
						 //pos[cx][cy][cz]=im;
						 xdr_float(&xdrs,&im); 
#ifdef DEBUG 
						 fprintf(f1," %d %d %d %f \n", i,j,k,im);
						 if(pos[i][j][k]!=0.){
							 fprintf(DBG," %d %d %d %f \n", i,j,k,im);
						 }
#endif
						 
					 }
				 }
			 }
		 }//if(encase)
		 else if(encase==1){
			 for(cz=0;cz<nz;cz++){ 
				 for(cy=0;cy<ny;cy++){ 
					 for(cx=0;cx<nx;cx++){ 
						 
						 //for(cy=nyy-1;cy>=0;cy--){ 
						 fgets(line, 180, f0);
						 countsites++;
						 sscanf(line, "%f", &im);
						 if(im==0){//pore
							 im=0.;
							 }
						 else if(im==125){//matrix
							 im=5.;
						 }
						 else if(im>0 && im <125){//undecided, is u2m or p2m set?
							 if(u2m){
								 im=5.;
							 }
							 else if(u2p){
								 im=0.;
							 }
							 else{}
						 }
						 else{
							 fprintf(stderr, "Error: Undefined value %f at position %i\n", im, countsites);
							 exit(-1);
						 }
					 //and now, restore channels and fluid:
						 //						 imc=cx;cx=cz;
						 if(  (cx<chanx)||(cx>=nxx-chanx)||(cy<chany)||(cy>=nyy-chany) ){
							 im=5.;
						 }
						 else if( (cz<chanz)||(cz>=nzz-chanz)  ){
							 im=0.;
						 }
						 //pos[cx][cy][cz]=im;
						 xdr_float(&xdrs,&im); 
#ifdef DEBUG 
						 fprintf(f1," %d %d %d %f \n", i,j,k,im);
						 if(pos[i][j][k]!=0.){
							 fprintf(DBG," %d %d %d %f \n", i,j,k,im);
						 }
#endif
						 
					 }
				 }
			 }
		 }//else(encase)

		 else if(encase==2){
				 //				 for(cx=nxx-1;cx>=0;cx--){ 
			 for(cz=0;cz<nz;cz++){ 
				 for(cy=0;cy<ny;cy++){ 
					 for(cx=0;cx<nx;cx++){ 

						 //for(cy=nyy-1;cy>=0;cy--){ 
						 if( (cz>=chanz)&&(cz<nz-chanz) ){						 
							 fgets(line, 180, f0);
							 countsites++;
							 sscanf(line, "%f", &im);
							 if(im==0){//pore
								 im=0.;
							 }
							 else if(im==125){//matrix
								 im=5.;
							 }
							 else if(im>0 && im<125){//undecided, is u2m or p2m set?
								 if(u2m){
									 im=5.;
								 }
								 else if(u2p){
								 im=0.;
								 }
								 else{}
							 }
							 else{
								 fprintf(stderr, "Error: Undefined value %f at position %i\n", im, countsites);
								 exit(-1);
							 }
						 }
						 //and now, restore channels and fluid:
						 //						 imc=cx;cx=cz;
						 if(  (cx<chanx)||(cx>=nx-chanx)||(cy<chany)||(cy>=ny-chany) ){
							 im=5.;
						 }
						 else if( (cz<chanz)||(cz>=nz-chanz)  ){
							 im=0.;
						 }
						 //pos[cx][cy][cz]=im;
						 xdr_float(&xdrs,&im); 
#ifdef DEBUG 
						 fprintf(f1," %d %d %d %f \n", i,j,k,im);
						 if(pos[i][j][k]!=0.){
							 fprintf(DBG," %d %d %d %f \n", i,j,k,im);
						 }
#endif
						 
					 }
				 }
			 }
		 }//else(encase)
		 
		 printf("number of sites copied:%i\n", countsites);


		 if(!zipflag){
			 fclose(f0);
		 }
		 else{
			 pclose(f0);
		 }

		 xdr_destroy(&xdrs); 
		 fclose(f2); 
#ifdef DEBUG 
		 fclose(f1);
		 fclose(DBG); 
#endif
  
	/*put it into the files */
/* 		 for(k=0;k<nz;++k){ */
/* 			 for(j=0;j<ny;++j){ */
/* 				 for(i=0;i<nx;++i){ */
/* 					 /\* Check output of xdr file*\/ */
/* 				 }      */
/* 			 } */
/* 		 } */
		 
	
	//write vtk file
	fprintf(VTK, "# vtk DataFile Version 2.0\n");
	fprintf(VTK, "test\nBINARY\nDATASET STRUCTURED_POINTS\n");
	fprintf(VTK, "DIMENSIONS %d %d %d\n",nx,ny,nz);
	fprintf(VTK, "ORIGIN 0 0 0\nSPACING 1 1 1\n");
	fprintf(VTK, "POINT_DATA %d\n", ((int) (nx*ny*nz) ));
	fprintf(VTK, "SCALARS OutArray float 1\nLOOKUP_TABLE default\n");
	fclose(VTK);

	sprintf(syscmd, "cat %s >> %s", xdrfile, vtkfile);
	system(syscmd);

	return 0;



         
	     }

