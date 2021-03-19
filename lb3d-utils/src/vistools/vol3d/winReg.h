
#ifndef INCLUDED_WINREG_H
#define INCLUDED_WINREG_H

struct winRegEnt {
    int win;
    void *ptr;
    struct winRegEnt *next;
};

int winReg_setPtr(int win, void *ptr);
int winReg_remove(int win);
void *winReg_getPtr(int win);

#endif /* INCLUDED_WINREG_H */
