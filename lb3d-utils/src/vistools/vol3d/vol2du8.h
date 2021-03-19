#ifndef INCLUDED_VOL2DU8_H
#define INCLUDED_VOL2DU8_H

#include <math.h>

#include "vol2d.h"

/* vol2du8.c */
struct vol2d *vol2du8_new(int nx, int ny);
#ifdef HAVE_PNG
struct vol2d *vol2du8_new_from_gspng(char *fname);
int vol2du8_write_png(struct vol2d *v, char *filename);
#endif /* HAVE_PNG */
struct vol2d *vol2du8_new_from_2df_normalize(struct vol2d *v);
struct vol2d *vol2du8_new_from_2df(struct vol2d *v,float pmax,float pmin);
#endif /* INCLUDED_VOL2DU8_H */
