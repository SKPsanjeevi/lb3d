#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpc/rpc.h>
#include <hdf5.h>
#include<string.h>

#define V_LIMIT (double) 0.1
#define MUSQARE2MDARCY(X)  ((double) 1000.0*X)/(double)0.986923


#undef SGL
#ifdef SGL
typedef  float REAL;
#else
typedef  double  REAL;
#endif

typedef struct{
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;
  REAL ***h5data;
} T_H5Array;

typedef struct{
  /*parameters*/
  REAL a;             /*lattice constant = resolution in mu*/
  REAL bforce;        /*accelerating body force*/
  unsigned int forcestart,forceend;
  REAL tau;           /*ralaxation time -> viscosity*/
  int  rockstart;     /*start index of the first z plane. planes are counted from 1,2,...,N*/
  int  rockend;       /* end index of last plan within porous sample , again sounted from 1,2,...,N*/
  unsigned int x,y,z; 
  unsigned int wallx;
  unsigned int wally;
  char vname[30];
  char odname[30];
} T_Params;

T_H5Array  hdfread(char *filename,unsigned int zplane) {
  hid_t       file_id, dataset_id,type_id;  /* identifiers */
  hid_t       dataspace;
  herr_t      status;
  hsize_t     dims[3];
  int         i,j,status_n,nz,ny,nx;
  REAL ***data=NULL;

  T_H5Array h5array;
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT); 
  /* Open an existing dataset. */
  dataset_id = H5Dopen1(file_id,"OutArray"); 
  /* Get dimensions */
  dataspace= H5Dget_space(dataset_id);
  status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
  printf("Dims from .h5 file %s :(x,y,z)  %lu x %lu x %lu \n",filename,(unsigned long)(dims[2]),(unsigned long)(dims[1]),(unsigned long)(dims[0]) );    
  nz=dims[0];
  ny=dims[1];
  nx=dims[2];
#ifdef SGL
  type_id=H5T_NATIVE_FLOAT;
#else
  type_id=H5T_NATIVE_DOUBLE;
#endif
  /*data=h5array.h5data;*/
  data=(double ***)malloc(nz*sizeof(double **));
  if(!data){fprintf(stderr,"Error: (1) can't allocate memory \n");
  exit(1);}
  data[0]=(double **)malloc(nz*ny*sizeof(double *));
  if(!data[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
  data[0][0]=(double *)malloc(nz*ny*nx*sizeof(double));
  if(!data[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}

  for(j=1;j<ny;j++){data[0][j]=data[0][j-1]+nx;}
  for(i=1;i<nz;i++){
    data[i]=data[i-1]+ny;
    data[i][0]=data[i-1][0]+ny*nx;
    for(j=1;j<ny;j++){
      data[i][j]=data[i][j-1]+nx;
    }
  } 
  
  status = H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   &data[0][0][0]);
  status = H5Dclose(dataset_id); 
  status = H5Fclose(file_id);
 
  h5array.nx=nx;
  h5array.ny=ny;
  h5array.nz=nz;
  h5array.h5data=data;
  return h5array;
}

void init_param(char * inputfile,char * simtag,REAL a,\
		unsigned int rs,unsigned int re,\
		unsigned int wallx, unsigned int wally,\
		T_Params *params){
  FILE *IF=NULL;
  char in[256];
  char tmpfile[30];

  if(rs>=re){
    printf("ERROR: rockstart>=rockend \n");
    abort();
  };
    
  params->a=a;
  params->rockstart=rs;
  params->rockend=re;
  params->wallx=wallx;
  params->wally=wally;
  printf("resolution a: %le\n",params->a);
  printf("rockstart: %u rockend: %u\n",params->rockstart,params->rockend);
  printf("Wall stength x:%d y:%d \n",wallx,wally);
  
  strncpy(tmpfile,"vel_v_",6);
  strncpy(tmpfile+6,simtag,18);
  strcpy(tmpfile+24,".h5");
  strncpy(params->vname,tmpfile,30);
  //printf(" vname: %s\n",params->vname);

  strncpy(tmpfile,"od_v_",5);
  strncpy(tmpfile+5,simtag,18);
  strcpy(tmpfile+23,".h5");
  strncpy(params->odname,tmpfile,30);
  //printf(" odname: %s\n",params->odname);

  IF=fopen(inputfile,"r");
    while(fgets(in,256,IF)){
      //printf("%s",in);
      if(!strncmp(in,"g_accn ",7)){
	//printf("found accn %s",in);
	sscanf(in,"g_accn = %lf",&(params->bforce));
	printf("body force is: %16le \n",params->bforce);
      }

      if(!strncmp(in,"g_accn_min ",11)){
	//printf("found g_accn_min %s",in);
	sscanf(in,"g_accn_min = %u",&(params->forcestart));
	printf("forcestart is: %u \n",params->forcestart);
      }

      if(!strncmp(in,"g_accn_max ",11)){
	//printf("found g_accn_max %s",in);
	sscanf(in,"g_accn_max = %u",&(params->forceend));
	printf("forceend is: %u \n",params->forceend);
      }

      if(!strncmp(in,"tau_b ",6)){
	//printf("found tau_b %s",in);
	sscanf(in,"tau_b = %lf",&(params->tau));
	printf("tau is: %16le \n",params->tau);
      }
      if(!strncmp(in,"nx ",3)){
	//printf("found nx %s",in);
	sscanf(in,"nx = %u",&(params->x));
	printf("nx is: %u \n",params->x);
      }
      if(!strncmp(in,"ny ",3)){
	//printf("found ny %s",in);
	sscanf(in,"ny = %u",&(params->y));
	printf("ny is: %u \n",params->y);
      }
      if(!strncmp(in,"nz ",3)){
	//printf("found nx %s",in);
	sscanf(in,"nz = %u",&(params->z));
	printf("nz is: %u \n",params->z);
      }
    }
};

int main(int argc, char* argv[]){
  unsigned int x=0,y=0,z=0,nx=0,ny=0,nz=0;
  REAL plane_vel=0.0, plane_rho=0.0,plane_mf=0.0;
  REAL rho_in=0.0, rho_out=0.0, dp=0.0, phi=0.0;
  REAL sum_zvel=0.0,sum_zz=0.0,min_z=V_LIMIT, max_z=0.0,avg_zvel=0.0,var_z=0.0;
  REAL sum_od=0.0,perm=0.0, vis=0.0, L=0.0,incount=0.0,outcount=0.0; ;
  unsigned long int limit_count=0,pore_count=0,rock_count=0, plane_count=0;

  hid_t FILE_OD,FILE_VEL;
  hid_t DATASET_OD,DATASET_VEL;
  hid_t DATASPACE_OD,DATASPACE_VEL;
  hid_t MEM_OD,MEM_VEL;
  hid_t type_id;

  herr_t status;
  hsize_t dim_od[3],dim_vel[3];
  hsize_t memdim[2],file_count[3],mem_count[2];
  hsize_t file_offset[3],mem_offset[2];

  T_Params params;

  char *input_file,*simtag;
  REAL *od=NULL,*vel=NULL;  


  if(argc!=8){
    printf("Usage: parmcal <input_file>  <simtag(=t123456-0123456789)> <resolution> <samplestart> <sampleend> <wallx> <wally>\n");
    exit(1);
  };
  
  input_file=argv[1];
  simtag=argv[2];
  init_param(input_file,simtag,atof(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),&params);
   
  /*Setting up file and dataspace handles*/
  
  FILE_OD       = H5Fopen(params.odname, H5F_ACC_RDONLY, H5P_DEFAULT);
  DATASET_OD    = H5Dopen1(FILE_OD,"OutArray");
  DATASPACE_OD  = H5Dget_space(DATASET_OD);
  status        = H5Sget_simple_extent_dims(DATASPACE_OD,dim_od,NULL);
  
  FILE_VEL      = H5Fopen(params.vname, H5F_ACC_RDONLY, H5P_DEFAULT);
  DATASET_VEL   = H5Dopen1(FILE_VEL,"OutArray");
  DATASPACE_VEL = H5Dget_space(DATASET_VEL);
  status        = H5Sget_simple_extent_dims(DATASPACE_VEL,dim_vel,NULL);

  /*get dimensions nx,ny,nz*/
  if( (dim_od[0]==dim_vel[0]) &&
      (dim_od[1]==dim_vel[1]) &&
      (dim_od[2]==dim_vel[2]) )
    {
      nz=dim_od[0];
      nx=dim_od[1];
      ny=dim_od[2];
      
    }else{
    printf("ERROR: dims of the datastructures in the  (HDF) velocity and density file do not agree \n");
    abort();
  };
 
  if( (nx!=params.x) || (ny!=params.y) || (nz!=params.z) ){
    printf("ERROR: One of the coordinates in the .h5 file is different from\
            the information defined in the input file:x=%d==%d y=%d==%d z=%d==%d",\
	   nx,params.x,ny,params.y,nz,params.z);
    abort();
  };


  /* Allocate  memory  structures for one z-plane  vel[][] and od[][]
   ans settin um HDF structures*/
  if(! (od=malloc(nx*ny*sizeof(REAL)))){
    printf("ERROR: unable to allocate memory for density data of size %d x %d \n",ny,nz);
    abort();
  };

  if(! (vel=malloc(nx*ny*sizeof(REAL)))){
    printf("ERROR: unable to allocate memory for velocity data of size %d x %d \n",ny,nz);
    abort();
  };

  memset(od,0,nx*ny*sizeof(REAL));
  memset(vel,0,nx*ny*sizeof(REAL));

  memdim[0]=nx;
  memdim[1]=ny;

  MEM_OD  = H5Screate_simple(2,memdim, NULL);
  MEM_VEL = H5Screate_simple(2,memdim, NULL);

 printf("---------------------------------------------------------------------------------------------------------\n");
  printf("-- Lattice constant in mu. Lattice time defined to 1 sec. All output compatible with these definitions --\n");
  printf("-- Precision compiled with real values of length in byte: %lu. Display in 16e.                         --\n", sizeof(REAL) );
  printf("---------------------------------------------------------------------------------------------------------\n");
  printf("[1]: z plane index starting with 1\n");
  printf("[2]: Average porespace velocity[l.u.] on plane z\n");
  printf("[3]: Average porespace density[l.u.]  on plane z\n");
  printf("[4]: Average porespace massflow[l.u.]=avg(rho(x,y)*v(x,y))  on plane z\n");
  printf("[5]: Plane porosity (porespace/x*y, walls not included) \n");
  printf("     [1]        [2]              [3]              [4]              [5]\n");

  for(z=0;z<nz;z++){
    
    // printf("starting new plane .................... %d\n",z);
    /* starting a new z plane */
    plane_count=0.0; 
    plane_vel=0.0;
    plane_rho=0.0;    
 
    /* Setting up hyperslab in file */
    file_offset[0] = z;
    file_offset[1] = 0;
    file_offset[2] = 0;
 
    file_count[0] = 1;
    file_count[1] = nx;
    file_count[2] = ny;

    status= H5Sselect_hyperslab(DATASPACE_OD, H5S_SELECT_SET,file_offset,NULL,file_count,NULL);
    status= H5Sselect_hyperslab(DATASPACE_VEL, H5S_SELECT_SET,file_offset,NULL,file_count,NULL);

    /* Setting up hyperslab in dataspace */
    mem_offset[0] = 0;
    mem_offset[1] = 0;
 
    mem_count[0] = nx;
    mem_count[1] = ny;
    
    status= H5Sselect_hyperslab(MEM_OD, H5S_SELECT_SET, mem_offset,NULL,mem_count,NULL);
    status= H5Sselect_hyperslab(MEM_OD, H5S_SELECT_SET,	mem_offset,NULL,mem_count,NULL);
    
    /* readin data from file */
#ifdef SGL
  type_id=H5T_NATIVE_FLOAT;
#else
  type_id=H5T_NATIVE_DOUBLE;
#endif

  status=H5Dread (DATASET_OD, type_id, MEM_OD, DATASPACE_OD, H5P_DEFAULT, od);
  status=H5Dread (DATASET_VEL, type_id, MEM_VEL, DATASPACE_VEL, H5P_DEFAULT, vel);

  for(y=0;y<ny;y++){
    for(x=0;x<nx;x++){
      unsigned long int plane_pos=(x*nx)+y;
        if(od[plane_pos]!=0.0){
	/*Fluid*/
	if(z==params.rockstart-1){
	  rho_in+=od[plane_pos]; incount++;
	}
	
	if(z==params.rockend-1){
	  rho_out+=od[plane_pos]; outcount++;
	}
	/* printf("%d %d %d: %16e ", x,y,z,od[z][x][y]);
	   if(x==nx-1)
	   printf("\n");
	*/
	if(vel[plane_pos]>=V_LIMIT){
	  limit_count++;
	}
	
	plane_count++; /* counts fluid foxels on a plane. also in bath*/
	
	plane_vel+=vel[plane_pos];
	plane_rho+=od[plane_pos];
	plane_mf=+(vel[plane_pos]*od[plane_pos]); 
	
	if( (z>=params.rockstart-1) && (z<=params.rockend-1)){
	  pore_count++; /*pore voxels. not reset fpor each plane. only in the sample*/
	  sum_zvel+=vel[plane_pos];
	  sum_zz+=(vel[plane_pos]*vel[plane_pos]);
	  sum_od+=od[plane_pos];
	  if(vel[plane_pos] < min_z)
	    min_z=vel[plane_pos];
	  if(vel[plane_pos] > max_z)
	    max_z=vel[plane_pos];	
	}
      }else{
	/*Matrix*/
	if(      (z >= params.rockstart-1) && (z <= params.rockend-1)  
		 && (y >= params.wally) && (y <= ny-params.wally-1)	
		 && (x >= params.wallx) && (x <= nx-params.wallx-1)
		 )
	  { 
	    //printf("%d %d %d: %16e \n", x,y,z,od[z][x][y]);
	    //if(od[z][x][y]!=0) abort();
	    rock_count++;
	    
	  }
      }
    } /* end of x */
  }/* end of y */
  
  plane_vel/=plane_count;
  plane_rho/=plane_count;
  plane_mf/=plane_count;
  
  if(plane_vel==0) {
    printf("ATTENTION !!!: average velocity  on plane %d is zero. Sample possibly not connected \n",z);
  }
  printf(" %6d %16le %16le %16le %16le\n",z+1,plane_vel,plane_rho,plane_mf, (REAL) plane_count/((nx-2*params.wallx)*(ny-2*params.wally)  )); /*plane statistics*/
  
  }/* end of z*/
    

    H5Dclose(DATASET_OD);
    H5Sclose(DATASPACE_OD);
    H5Sclose(MEM_OD);
    H5Fclose(FILE_OD);

    H5Dclose(DATASET_VEL);
    H5Sclose(DATASPACE_VEL);
    H5Sclose(MEM_VEL);
    H5Fclose(FILE_VEL);

   /*
    IN = first sample cross section
    OUT= last sample cross section

    note:!! for the pressure gradient not the cross sections BEFORE/AFTER the sample are used.
         BUT the FIRST/LAST cross section  within  the sample.
   */

  rho_in/=incount;                                     /* mean inlet Density  */
  rho_out/=outcount;                                   /* mean inlet pressure */
  //printf("-------- %16le  %16le \n",rho_in/3.0,rho_out/3.0);
  L=(params.rockend-params.rockstart);                 /* in lu */   
  phi= (REAL) pore_count/(rock_count+pore_count);      /* dimensionless */
  //avg_zvel=params.a*sum_zvel/pore_count;             /* mu/s  */
  avg_zvel=sum_zvel/(pore_count+rock_count);
  sum_zz=sum_zz/(pore_count+rock_count);
  var_z=sum_zz-avg_zvel*avg_zvel;                       /* variance !WARNING! not unbiased variance as in many comp. tools. */
  dp=(rho_out-rho_in)/(3.0*L);                          /* global pressure gradient Pa/mu */
  vis=(2.0*params.tau-1.0)/6.0*(sum_od/pore_count);     /* dyn. vis. in Pa*s=kg/(m*s)     */
  perm=-vis*avg_zvel/dp*(params.a)*(params.a);          /* in mu^2 */
  
  printf("Full sample Nof. pores: %d  Nof. matrix: %d  Phi: %16le \n",  (int) pore_count,(int) rock_count,phi);
  printf("Full sample z-velocity average[l.u.]         :  %16le\n",avg_zvel);
  printf("Full sample z-velocity variance[l.u.]        :  %16le\n",var_z);
  printf("Full sample (min-max) z-velocity[l.u.]       : (%16le , %16e)\n",min_z,max_z);
  printf("Number/Perc. of nodes with z-vel. abv. limit :  %lu / %16le \n", limit_count, (REAL)limit_count/(REAL)(pore_count+rock_count));
  printf("Full sample Dyn.viscosity [l.u.]             :  %16le\n",vis);
  printf("Full sample global pressure gradient [l.u.]  :  %16le\n",dp);
  printf("Full sample permeability[mu^2,mD]            :  %16le , %16le\n",perm,MUSQARE2MDARCY(perm));

  exit(0);
}

