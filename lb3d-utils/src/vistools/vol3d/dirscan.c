#include "dirscan.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <dirent.h>

#include <sys/types.h>
#include <sys/stat.h>

/* Scan the current directory for files which begin with "prefix",
 * which are older than minAge seconds, but which have been modified
 * since lastTime. If any such files exist, return the name of the
 * newest; otherwise return NULL.
 */

char *findRipeFile(char *prefix, time_t lastTime, time_t minAge)
{
    DIR *dir=NULL;
    time_t now;
    struct dirent *dirent=NULL;
    size_t prefixLen=0;
    struct stat statbuf;

    time_t newestFileMtime=0,minAgeTime=0;
    size_t newestFilenameLen=0;
    char  *newestFilename=NULL;

    prefixLen=strlen(prefix);

    if (NULL==(dir=opendir("."))) { perror("opendir"); return NULL; }

    now=time(NULL);
    minAgeTime=now-minAge; /* mtime must be < minAgeTime */

    while (NULL!=(dirent=readdir(dir))) {
        if (0!=strncmp(prefix,dirent->d_name,prefixLen)) {continue;}

        if (0!=stat(dirent->d_name,&statbuf))
            { perror("stat"); continue; }


        if (statbuf.st_mtime < minAgeTime) { /* older than minAge? */
            if (statbuf.st_mtime > lastTime) { /* Newer than lastTime? */
                if (statbuf.st_mtime > newestFileMtime) { /* Newest of these? */

                    /* We have a new winner. Take a copy of its filename,
                     * and update newestFileMtime.
                     */

                    newestFileMtime = statbuf.st_mtime;
                    newestFilenameLen=strlen(dirent->d_name)+1;
                    if (NULL==(newestFilename=realloc(newestFilename,
                                    newestFilenameLen)))
                        { perror("realloc"); return NULL; }

                    strncpy(newestFilename,dirent->d_name,newestFilenameLen);
         }}}
    }

    closedir(dir);

    return newestFilename;
}


#ifdef _TEST_DIRSCAN_C

int main(int argc, char *argv[])
{
    time_t minAge=3; /* Default 3 seconds */
    char *prefix=NULL;
    char *newestFilename=NULL;
    time_t newestTime=0;
    struct stat statbuf;

    if ((1==argc)||(3<argc)) {
        fprintf(stderr,"Syntax: %s <prefix> [<minAge>]\n",argv[0]);
        return -1;
    }

    if (3==argc) { minAge = (time_t) atoi(argv[2]); }

    prefix=argv[1];

    printf("Waiting for a file to appear...\n");fflush(stdout);
    while(1) {
        if (NULL!=(newestFilename=findRipeFile(prefix,newestTime,minAge))) {
            if (0!=stat(newestFilename,&statbuf)) { perror("stat"); return -1; }
            newestTime = statbuf.st_mtime;
            printf("Ripest is <%s>\n",newestFilename);
            free(newestFilename);
        } else {
            sleep(1);
        }

    }

    return 0;
}

#endif /* _TEST_DIRSCAN_C */
