/*  3-D VERSION OF THE LIKE
 *  PROGRAM BY J. CHIN
 *  
 *  This program handles double and single precision data files and only
 *  functions '3' and '4' are working at the moment.
 */


#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <math.h>
#include <string.h>

#define BUFSIZE         1024 /* So shoot me. */
#define TIMECHARSTR     "6" 
#define TIMECHARSTRLEN  6
#define FNAMEBUFSIZE    512

/* (char *)find_path(); */

int main (int argc, char *argv[])
{
   FILE	     *f, *f2;
   int	     i;
   XDR	     xdrs, xdrs2;
   char      *fname, *fname2, *aux, *file, *path;
   
   int       nx,ny,nz,type,action;

#ifdef SGL
   float    buf,buf_b,buf_c,buf_all,buf2;
   float    phimax,phimin,phisum,phisum2,phisum_abs,phisum_abs2;
   float    phimax_b,phimin_b,phisum_b,phisum2_b;
   float    phimax_c,phimin_c,phisum_c,phisum2_c;
   float    phimax_all,phimin_all,phisum_all,phisum2_all;
   float    sd, sd_b, sd_c, sd_all, sd_abs;
   float    max, dens, norm2, sq_sum;
   float    speed,speed_b,speed_c,speedavg,speedavg2,stderr_speed;
   float    speedmax = 0.0;
#else
   double    buf,buf_b,buf_c,buf_all,buf2;
   double    phimax,phimin,phisum,phisum2,phisum_abs,phisum_abs2;
   double    phimax_b,phimin_b,phisum_b,phisum2_b;
   double    phimax_c,phimin_c,phisum_c,phisum2_c;
   double    phimax_all,phimin_all,phisum_all,phisum2_all;
   double    sd, sd_b, sd_c, sd_all, sd_abs;
   double    max, dens, norm2, sq_sum;
   double    speed,speed_b,speed_c,speedavg,speedavg2,stderr_speed;
   double    speedmax = 0.0;
#endif
   
   fname2 = (char *)malloc(FNAMEBUFSIZE * sizeof(char));
   file   = (char *)malloc(FNAMEBUFSIZE * sizeof(char)); 
   path   = (char *)malloc(FNAMEBUFSIZE * sizeof(char));
   aux    = (char *)malloc(FNAMEBUFSIZE * sizeof(char));
   strcpy(fname2,"\0");
   strcpy(file,"\0");
   strcpy(path,"\0");
   strcpy(aux,"\0");


   if (argc!=7) {
      fprintf(
	 stderr,"\nSyntax: "
	 "\t%s <filename> <nx> <ny> <nz> <type> <action>         \n\n"
	 "\twhere:                                               \n\n"
	 "\tfilename   is the name of the data file, including     \n"
	 "\t           full path.                                  \n"
	 "\ttype==1, 2, or 3 depending on the no. of fields per \n"
	 "\t           space point.                                \n"
	 "\taction==1  give full statistics on each field and      \n"
	 "\t           their sum; valid for any <type>.            \n" 
	 "\taction==2  give max{abs(field(x)), x in lattice},      \n"
	 "\t           valid for <type>=1 only.                    \n" 
	 "\taction==3  give 'avg{abs(field(x))} stderr(avg)',      \n"
	 "\t           valid for <type>=1 only.                    \n" 
	 "\taction==4  give 'avg{|(g1,g2,g3)|} stderr(avg)',       \n"
	 "\t           where |.| is a vector norm, and             \n"
	 "\t           g? = arr?(x)/owd(x). Here arr?(x) is        \n"
	 "\t           the field ? of file arr*.xdr and owd(x) is  \n"
	 "\t           field1+field2 of file owd*.xdr, both taken  \n"
	 "\t           at the same lattice site x. Valid for       \n"
	 "\t           <filename>=arr*.xdr (<type>=3) only.        \n"
	 "\taction==5  give 'max{|(g1,g2,g3)|}',                   \n"
	 "\t           same symbols than action==4. Valid for      \n"
	 "\t           <filename>=arr*.xdr (<type>=3) only.      \n\n"
	 "\t IMPORTANT:                                          \n\n" 
	 "\t    The only options giving the correct std. err.      \n"
	 "\t    of the avg. is type=3, action=4. The rest give     \n"
	 "\t    the std.dev. of the sample distribution instead    \n"
	 "\t    WHICH IS UTTERLY WRONG!!                         \n\n"
	 "\t PLEASE DO NOT NAME ANY DIRECTORY IN THE RELEVANT      \n"
	 "\t FILE'S TREE AS 'arr_', AS THIS IS THE PATTERN THIS    \n"
	 "\t CODE USES TO SEARCH FOR THE FILENAME FROM THE         \n"
	 "\t FULL PATH FILENAME.                                 \n\n"
	 "\t    The code needs to be edited each time it's applied \n"
	 "\t    to a different field, namely, changing the field's \n"
	 "\t    name in the strstr() function call.\n\n"
	 ,argv[0]);
      return -1;
   }
   
   fname=argv[1];
   nx=atoi(argv[2]);
   ny=atoi(argv[3]);
   nz=atoi(argv[4]);
   type=atoi(argv[5]);
   action=atoi(argv[6]);

/* DEBUG */
   printf("\nParameters loaded");
/* END */

/* Determine filename without path: */
   file=strstr(fname,"dir_");
   
/* DEBUG */
   printf("\nfname=%s"
	  "\nfile=%s",fname,file);
/* END */

/* Extract path without filename. 
// DEBUG
//   fprintf(stderr,"\npath=%s\n"
//	   "strlen(fname)-strlen(file)-1=%d\n",
//	   path,strlen(fname)-strlen(file)-1);
// END
*/
   getchar();
   strncat(path,fname,strlen(fname)-strlen(file)-1);



/*   exit(0); */

   if(action==4 || action==5) {
/*    Next lines obtain owd*.xdr file name
//    from the arr*.xdr file name.
*/

      strncat(aux,file+3,strlen(file)-3);
      sprintf(fname2,"%s/owd%s",path,aux);
/*      printf("\naux  = %s",aux);
//      printf("\nfname2 = %s",fname2);
//      exit(0);
*/
   }
   

   fprintf(stderr,"\nOpening data file %s\n",fname);	 
   if (NULL==(f=fopen(fname,"r"))) {
      fprintf(stderr,"... Unable to open data file.");
      return -1;
   }
   
   if(action==4 || action==5) {
      fprintf(stderr,"\nOpening data file %s\n",fname2);
      if (NULL==(f2=fopen(fname2,"r"))) {
	 fprintf(stderr,"... Unable to open owd file.");
	 return -1;
      }
   }
   
   xdrstdio_create(&xdrs,f,XDR_DECODE);
   
   if(action==4 || action==5) 
      xdrstdio_create(&xdrs2,f2,XDR_DECODE);
   
/*      Initialise accumulators */

#ifdef SGL
   xdr_float(&xdrs,&buf);
#else
   xdr_double(&xdrs,&buf);
#endif
   phisum=phimax=phimin=buf;
   phisum2=buf*buf;
   phisum_abs=fabs(buf);
   phisum_abs2=fabs(buf)*fabs(buf);
   if(action==4 || action==5) {
#ifdef SGL
      xdr_float(&xdrs2,&buf2);
#else
      xdr_double(&xdrs2,&buf2);
#endif
      dens=buf2;
#ifdef SGL
      xdr_float(&xdrs2,&buf2);
#else
      xdr_double(&xdrs2,&buf2);
#endif
      dens+=buf2;
      if(dens!=0.0)
	 speed=buf/dens;
      else
	 fprintf(stderr,"\n**Zero density at first record**");
   }
   if(type==2) {
#ifdef SGL
      xdr_float(&xdrs,&buf_b); 
#else
      xdr_double(&xdrs,&buf_b); 
#endif
      phisum_b=phimax_b=phimin_b=buf_b;
      phisum2_b=buf_b*buf_b;
      
      buf_all=buf + buf_b;
      phisum_all=phimax_all=phimin_all=buf_all;
      phisum2_all=buf_all*buf_all;
   }
   else if(type==3) {
#ifdef SGL
      xdr_float(&xdrs,&buf_b); 
#else
      xdr_double(&xdrs,&buf_b); 
#endif
      phisum_b=phimax_b=phimin_b=buf_b;
      phisum2_b=buf_b*buf_b;
      
      if(action==4 || action==5) 
	 speed_b=buf_b/dens;	     
      
#ifdef SGL
      xdr_float(&xdrs,&buf_c); 
#else
      xdr_double(&xdrs,&buf_c); 
#endif
      phisum_c=phimax_c=phimin_c=buf_c;
      phisum2_c=buf_c*buf_c;
      
      if(action==4 || action==5)
	 speed_c=buf_c/dens;	     
      
      buf_all=buf + buf_b + buf_c;
      phisum_all=phimax_all=phimin_all=buf_all;
      phisum2_all=buf_all*buf_all;
   }
   if(action==4) {
      speedavg2=speed*speed+speed_b*speed_b+speed_c*speed_c;
      speedavg=sqrt(speedavg2);
   }
   if(action==5) {
      speedmax = speed*speed+speed_b*speed_b+speed_c*speed_c;
   }
   
/*      Loop over data */

   for (i=0;i<(nx*ny*nz-1);i++) {
#ifdef SGL
      xdr_float(&xdrs,&buf);
#else
      xdr_double(&xdrs,&buf);
#endif
      phisum+=buf;
      phisum2+=buf*buf;
      phisum_abs+=fabs(buf);
      phisum_abs2+=fabs(buf)*fabs(buf);
      if (buf>phimax) { phimax=buf;}
      if (buf<phimin) { phimin=buf;}
      
      if(action==4 || action==5) {
#ifdef SGL
	 xdr_float(&xdrs2,&buf2);
#else
	 xdr_double(&xdrs2,&buf2);
#endif
	 dens=buf2;
#ifdef SGL
	 xdr_float(&xdrs2,&buf2);
#else
	 xdr_double(&xdrs2,&buf2);
#endif
	 dens+=buf2;
	 if(dens!=0.0)
	    speed=buf/dens;	     
	 else
	    fprintf(stderr,"\n**Zero density at rec=%d**",i); 
      }
      
      if(type==2) {
#ifdef SGL
	 xdr_float(&xdrs,&buf_b); 
#else
	 xdr_double(&xdrs,&buf_b); 
#endif
	 phisum_b+=buf_b;
	 phisum2_b+=buf_b*buf_b;
	 if (buf_b>phimax_b) { phimax_b=buf_b;}
	 if (buf_b<phimin_b) { phimin_b=buf_b;}
	 
	 buf_all=buf + buf_b;
	 phisum_all+=buf_all;
	 phisum2_all+=buf_all*buf_all;
	 if (buf_all>phimax_all) { phimax_all=buf_all;}
	 if (buf_all<phimin_all) { phimin_all=buf_all;}
	 
      }
      if(type==3) {
#ifdef SGL
	 xdr_float(&xdrs,&buf_b); 
#else
	 xdr_double(&xdrs,&buf_b); 
#endif
	 phisum_b+=buf_b;
	 phisum2_b+=buf_b*buf_b;
	 if (buf_b>phimax_b) { phimax_b=buf_b;}
	 if (buf_b<phimin_b) { phimin_b=buf_b;}
	 
	 if(action==4 || action==5) 
	    speed_b=buf_b/dens;	     
	 
#ifdef SGL
	 xdr_float(&xdrs,&buf_c); 
#else
	 xdr_double(&xdrs,&buf_c); 
#endif
	 phisum_c+=buf_c;
	 phisum2_c+=buf_c*buf_c;
	 if (buf_c>phimax_c) { phimax_c=buf_c;}
	 if (buf_c<phimin_c) { phimin_c=buf_c;}
	 
	 if(action==4 || action==5) 
	    speed_c=buf_c/dens;	     
	 
	 buf_all=buf + buf_b + buf_c;
	 phisum_all+=buf_all;
	 phisum2_all+=buf_all*buf_all;
	 if (buf_all>phimax_all) { phimax_all=buf_all;}
	 if (buf_all<phimin_all) { phimin_all=buf_all;}
      }
      
      if(action==4) {
	 norm2=speed*speed+speed_b*speed_b+speed_c*speed_c;
	 speedavg2+=norm2;
	 speedavg+=sqrt(norm2);
      }
      if((action==5) && (speed*speed+speed_b*speed_b+speed_c*speed_c > speedmax))
	 speedmax = speed*speed+speed_b*speed_b+speed_c*speed_c;
   }
   
   xdr_destroy(&xdrs);
   if(action==4 || action==5)
      xdr_destroy(&xdrs2);
   
#ifdef SGL
   phisum/=(float)(nx*ny*nz);
   phisum2/=(float)(nx*ny*nz);
   phisum_abs/=(float)(nx*ny*nz);
   phisum_abs2/=(float)(nx*ny*nz);
#else
   phisum/=(double)(nx*ny*nz);
   phisum2/=(double)(nx*ny*nz);
   phisum_abs/=(double)(nx*ny*nz);
   phisum_abs2/=(double)(nx*ny*nz);
#endif

   sd=sqrt(phisum2-phisum*phisum);
   sd_abs=sqrt(phisum_abs2-phisum_abs*phisum_abs);

   if(action==4) {
#ifdef SGL
      speedavg2 /= (float)(nx*ny*nz-1);
      sq_sum = speedavg * speedavg / (float)(nx*ny*nz*(nx*ny*nz-1));
      stderr_speed=sqrt((speedavg2 - sq_sum)/nx*ny*nz);
      speedavg/=(float)(nx*ny*nz);
#else
      speedavg2 /= (double)(nx*ny*nz-1);
      sq_sum = speedavg * speedavg / (double)(nx*ny*nz*(nx*ny*nz-1));
      stderr_speed=sqrt((speedavg2 - sq_sum)/nx*ny*nz);
      speedavg/=(double)(nx*ny*nz);
#endif
/*
//    std.err of avg. = \sigma_{N-1} / \sqrt{N}
//      
//    where \sigma_{N-1} is the std.dev. of the sample
//    distribution, for N points in the sample. The std.err. of avg.
//    is the std.dev. of the distribution of the avg., the latter
//    being the best estimator for the mean.
//    DEBUG
//      printf("\nspeedavg2 >? speedavg^2 :"
//	     " %12.6g >? %12.6g\n",speedavg2,speedavg*speedavg);
//    END DEBUG
*/

   }
   
   if(type==2) {
#ifdef SGL
      phisum_b/=(float)(nx*ny*nz);
      phisum2_b/=(float)(nx*ny*nz);
      phisum_all/=(float)(nx*ny*nz);
      phisum2_all/=(float)(nx*ny*nz);
#else
      phisum_b/=(double)(nx*ny*nz);
      phisum2_b/=(double)(nx*ny*nz);
      phisum_all/=(double)(nx*ny*nz);
      phisum2_all/=(double)(nx*ny*nz);
#endif

      sd_b=sqrt(phisum2_b-phisum_b*phisum_b);
      sd_all=sqrt(phisum2_all-phisum_all*phisum_all);
   }
   
   if(type==3) {
#ifdef SGL
      phisum_b/=(float)(nx*ny*nz);
      phisum2_b/=(float)(nx*ny*nz);
      phisum_c/=(float)(nx*ny*nz);
      phisum2_c/=(float)(nx*ny*nz);
      phisum_all/=(float)(nx*ny*nz);
      phisum2_all/=(float)(nx*ny*nz);
#else
      phisum_b/=(double)(nx*ny*nz);
      phisum2_b/=(double)(nx*ny*nz);
      phisum_c/=(double)(nx*ny*nz);
      phisum2_c/=(double)(nx*ny*nz);
      phisum_all/=(double)(nx*ny*nz);
      phisum2_all/=(double)(nx*ny*nz);
#endif

      sd_b=sqrt(phisum2_b-phisum_b*phisum_b);
      sd_c=sqrt(phisum2_c-phisum_c*phisum_c);
      sd_all=sqrt(phisum2_all-phisum_all*phisum_all);
   }
   
   if(action==1) {
      printf("FIELD 1: MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n",
	     phimax,phimin,phisum,sd);
      
      if(type==2) {
	 printf("FIELD 2: MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n"
		,phimax_b,phimin_b,phisum_b,sd_b);
	 if(action==1)
	    printf("FIELD 1+FIELD 2: "
		   "MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n"
		   ,phimax_all,phimin_all,phisum_all,sd_all);
      }
      if(type==3) {
	 printf("FIELD 2: MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n"
		,phimax_b,phimin_b,phisum_b,sd_b);
	 printf("FIELD 3: MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n"
		,phimax_c,phimin_c,phisum_c,sd_c);
	 if(action==1)
	    printf("FIELD 1+FIELD 2+FIELD 3: "
		   "MAX %12.6g MIN %12.6g AVG %12.6g STDERR %12.6g\n"
		   ,phimax_all,phimin_all,phisum_all,sd_all);	      
      }
   }
   else if(action==2) {
      if(fabs(phimax)>=fabs(phimin))
	 max=fabs(phimax);
      else
	 max=fabs(phimin);
      printf("%12.6g",max);
   }
   else if(action==3) 
      printf("%12.6g %12.6g",phisum_abs,sd_abs);
   else if(action==4) 
      printf("%12.6g %12.6g",speedavg,stderr_speed);
   else if(action==5)
      printf("%12.6g",speedmax);

   fclose(f);
   if(action==4 || action==5)
      fclose(f2);

   return(1);
}










