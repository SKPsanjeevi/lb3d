#include<stdio.h>
#include<float.h>
#include "main.h"

int
extrema(Stats *stats, float *field) {

   static float  min = FLT_MAX, max = -FLT_MAX;

   if(*field > max)  {
      max = *field;
      stats->max = *field;
   }
   if(*field < min)  {
      min = *field;
      stats->min = *field;
   }

   return 1;
}
