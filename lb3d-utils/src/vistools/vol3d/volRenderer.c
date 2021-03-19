#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <rpc/rpc.h>

#include "volRender.h"

static void parse_2floats_string(float *a,float *b,char *str)
{
        char *s=NULL,*p;
        size_t len;
        int i=0;

        len=strlen(str);
        if (NULL==(s=malloc(len+1))) {
                perror("parse_2floats_string: malloc");
                exit(-1);
        }

        strncpy(s,str,len+1);

        /* Read first digit */
        while ((s[i]=='.') || (s[i]=='-') || isdigit(s[i])) { i++; }
        if (',' != s[i]) {
                fprintf(stderr,"parse_2floats_string: Bad input\n");
                exit(-1);
        }
        s[i++]=0;
        if (i>=len) {
                fprintf(stderr,"parse_2floats_string: Bad input\n");
                exit(-1);
        }
        *a = (float) atof(s); p=s+i;

        /* Read second digit */
        while ((s[i]=='.') || (s[i]=='-') || isdigit(s[i])) { i++; }
        s[i++]=0;
        *b = (float) atof(p); 

        free(s);
}

int main(int argc, char *argv[])
{


    struct volRender *self=NULL;
    char *stateFilename=NULL;
    GLfloat bgColor[4];
    int setBGColor=0;
    int stereo=0;
    int nFiles=0;
    int i,c;

    int loadNewest=0;
    char *loadNewestPrefix=NULL;
    float rangeMax,rangeMin;
    int useRange=0;

    const char *optstring="sb:r:l:N:";

    while (-1!=(c=getopt(argc,argv,optstring))) {
        switch(c) { 
            case 's':
                stereo=1;
                break;
            case 'b':
                {
                    float red,grn,blu;
                    vol3d_parse_floats_string( &red,&grn,&blu,optarg);
                    bgColor[0] = (GLclampf) red;
                    bgColor[1] = (GLclampf) grn;
                    bgColor[2] = (GLclampf) blu;
                    bgColor[3] = 1.0;
                    setBGColor=1;
                }
                break;
            case 'r':
                {
                    parse_2floats_string(&rangeMax,&rangeMin,optarg);
                    if (rangeMax<rangeMin) {
                        float tmp;
                        tmp=rangeMax;
                        rangeMax = rangeMin;
                        rangeMin = tmp;
                    }
                    useRange=1;
                    fprintf(stderr,"Using range %f -- %f\n",
                            rangeMin,rangeMax);
                }
                break;
            case 'l':
                stateFilename = optarg;
                break;
            case 'N':
                loadNewest=1;
                loadNewestPrefix=optarg;
                break;

            default: return -1;
        } /* switch(c) */
    } /* while (getopt()) */

    if (loadNewest) {
        nFiles=1;
    } else {
        nFiles = argc-optind;

        if (0==nFiles) {
                fprintf(stderr,"%s <options> <infiles>\n",argv[0]);
                return -1;
        }
    }

    initCubeVectors();

    glutInit(&argc,argv);

    if (NULL==(self=volRender_new(nFiles)))
        { fprintf(stderr,"volRender_new() failed.\n"); return -1; }

    self->rangeMax=rangeMax;
    self->rangeMin=rangeMin;
    self->useRange=useRange;
    self->prefix = loadNewestPrefix;


    if (loadNewest) { /* Wait for a file. */
        char *newestFilename=NULL;
        struct stat statbuf;
        time_t newestTime=0;

        /* Load the newest file */


        /* Wait for a file.. */

        fprintf(stderr,"Waiting for a file to appear...\n");fflush(stderr);
        while (NULL==(
                newestFilename=findRipeFile(self->prefix,self->newestTime,self->minAge))) {
            sleep(1);
        }

        if (0!=stat(newestFilename,&statbuf)) { perror("stat"); return -1; }

        newestTime = statbuf.st_mtime;
        fprintf(stderr,"Newest file is %s\n",newestFilename);

        if (0!=volRender_loadVolume(self,newestFilename,0)) {
            fprintf(stderr,"volRender_loadVolume(%s) failed.\n",
                    newestFilename);
            return -1;
        }

        if (NULL==(self->name=malloc(sizeof(char *))))
                { perror("malloc"); return -1; }
        self->name[0] = newestFilename;
        self->newestTime = newestTime;

    } else {
        /* Loop over each filename on the command line. */

        for (i=0;i<nFiles;i++) {

            if (0!=volRender_loadVolume(self,argv[optind+i],i)) {
                fprintf(stderr,"volRender_loadVolume(%s) failed.\n",
                        argv[optind+i]);
                return -1;
            }
        }

        self->name = argv+optind;
    }

    if (NULL==(self->lut=u8_rgba_lut_default())) {
        fprintf(stderr,"u8_rgba_lut_default failed.\n");
        return -1;
    }

    if(NULL==(self->vrgba=vol3d_u8_lut_to_rgba(self->quantized[self->currVol],self->lut))){
        fprintf(stderr,"vol3d_u8_lut_to_rgba failed.\n");
        return -1;
    }


    if (NULL==(self->lw=lutWidget_new(self->ft[0],self->lut))) {
        fprintf(stderr,"lutWidget_new failed.\n");
        return -1;
    }

    lutWidget_setKeyCallback(self->lw,myKeyCallback,self);

    if (NULL==(self->gv=(stereo ? glutViewer_new_stereo() : glutViewer_new())))
    { fprintf(stderr,"glutViewer_new() failed.\n"); return -1; }

    /* Point glutViewer back at parent so we can find the corresponding
     * volRender object during a callback.
     */

    self->gv->data = self; 
    self->fogOn=0;
    self->wrapFrames=0;
    self->gv->fovy=20.0;

    /* Enable alpha channel */

    self->gv->alphaBufferFlag=1;

    adjFloatList_addFloat(self->afl,&self->gv->fovy,
            0.0,360.0,1.0,"FOV");
    if (self->gv->stereoFlag) {
        adjFloatList_addFloat(self->afl,&self->gv->eyesize,
                0.0,0.2,0.002,"eyesize");
    }


    /* Deserialize state from file if necessary */

    if (NULL!=stateFilename) {
        if (0!=volRender_serializeFile(self,stateFilename,VOLRENDER_READ)) {
            fprintf(stderr,"Failed to load state from %s\n",stateFilename);
        }
        else {
            int winOld;
            volRender_reload(self);
            winOld=glutGetWindow();
            glutSetWindow(self->lw->win);
            lutWidget_makeTFGraph(self->lw);
            glutPostRedisplay();
            glutSetWindow(winOld);
        }
        /* Continue nonetheless. */
    }

    if (1==setBGColor) {
        int j;
        for (j=0;j<4;j++) { self->gv->bgColor[j] = bgColor[j]; }
        fprintf(stderr,"Set BG color to %f %f %f\n",
                self->gv->bgColor[0],
                self->gv->bgColor[1],
                self->gv->bgColor[2]);
    }


    glutMotionFunc(volRender_motionFunc);

    glutViewer_setDisplayCallback(self->gv,myDisplay);
    glutViewer_setDisplayCallbackData(self->gv,self);

    glutViewer_wheelUpFunc(self->gv,wheelUpFunc);
    glutViewer_wheelDownFunc(self->gv,wheelDownFunc);

    /* Initialize GL context for viewer */

    glutSetWindow(self->gv->win);

    if (loadNewest) { /* Set up scanner callback */
        glutTimerFunc(1000*self->minAge,volRender_scannerTimerFunc,
                self->gv->win);
    }


    glutKeyboardFunc(volRender_keyboardFunc);

    glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
    glDisable(GL_LIGHTING);

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);

    volRender_setTexture(self);



    glutMainLoop();
    printf("# Done!\n");

    return 0;

}       
