#ifndef INCLUDED_LUT_H
#define INCLUDED_LUT_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

/* Structure holding a frequency table for data which has 
 * unsigned 8-bit data. Each of the 256 possible byte values has an
 * 32-bit frequency indicating the number of voxels with that value.
 * Also contains max and min values, corresponding to byte values
 * of 255 and 0.
 */
struct u8_freqtable {
    uint32_t freq[256];
    float max,min;
};

/* RGBA Look Up Table entry and table */

struct rgba_lutent {
    GLubyte c[4];
};

struct u8_rgba_lut {
    struct rgba_lutent ent[256];
};

#endif /* INCLUDED_LUT_H */
