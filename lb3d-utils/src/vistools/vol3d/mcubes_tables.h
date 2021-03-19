#ifndef INCLUDED_MCUBES_TABLES_H
#define INCLUDED_MCUBES_TABLES_H
#include "mcubes.h"

extern int edgeTable[256];
extern int triTable[256][16];
extern const struct vertex cube_vert[8];
extern const int cube_edge[12][2];
extern const int cube_edge_delta[12][4];

#endif /* INCLUDED_MCUBES_TABLES_H */
