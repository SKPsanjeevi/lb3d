# include <stdio.h>

main(int argc, char *argv[])
{
  int cnt;
 char cfr[80];
 FILE *stigfmt,*rawfmt;

 stigfmt = fopen(argv[1],"r");  
 sprintf(cfr,"%s.raw",argv[1]); 
 rawfmt = fopen(cfr,"w");      

 while (fscanf(stigfmt,"%ld",&cnt) != EOF){
   fputc(cnt,rawfmt);
 }

 fclose(stigfmt);
 fclose(rawfmt);
}
