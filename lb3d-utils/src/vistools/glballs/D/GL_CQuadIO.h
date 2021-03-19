#ifndef GL_CQUADIO_H
#define GL_CQUADIO_H

#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>

#include "GL_CQuad.h"


class GL_CQuadIO{
 public:
  GL_CQuadIO(const char *path, const char *mode = "r"){
   m_glquadfile=fopen(path,mode);
   if(!m_glquadfile)perror("QuadIO :");
  };
  ~GL_CQuadIO(void){
   if(m_glquadfile)fclose(m_glquadfile); 
  }
  GL_CQuad readquad(){
   GL_CQuad quad;
   fread(&quad,sizeof(GL_CQuad),1,m_glquadfile);
   if(errno)perror("QuadIO :"); 
   return quad;
  }
  int readnumberofquads(){
   int number;
   fread(&number,sizeof(int),1,m_glquadfile);
   if(errno)perror("QuadIO :"); 
   return number;
  }
  void writequad(const GL_CQuad quad){
   fwrite(&quad,sizeof(GL_CQuad),1,m_glquadfile);
   if(errno)perror("QuadIO :"); 
  }
  void writenumberofquads(const int number){
   m_fseek=ftell(m_glquadfile);
   fwrite(&number, sizeof(int),1,m_glquadfile);
  }

  void rewritenumberofquads(const int number){
   long pos=ftell(m_glquadfile);
   fseek(m_glquadfile,m_fseek,SEEK_SET);
   fwrite(&number, sizeof(int),1,m_glquadfile);
   fseek(m_glquadfile,pos,SEEK_SET);
  }

 private:
  FILE *m_glquadfile;
  long m_fseek;
};
#endif
