#ifndef INCLUDED_SCREENSHOT_H
#define INCLUDED_SCREENSHOT_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <png.h>

int gl_screenshot_write_png(unsigned char *pixdata, int nx, int ny, char *filename, int alphaFlag);
int gl_screenshot_png(int width, int height, char *filename, int alphaFlag);

#endif /* INCLUDED_SCREENSHOT_H */
