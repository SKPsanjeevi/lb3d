#ifndef INCLUDED_DIRSCAN_H
#define INCLUDED_DIRSCAN_H
#include <time.h>
char *findRipeFile(char *prefix, time_t lastTime, time_t minAge);
#endif /* INCLUDED_DIRSCAN_H */
