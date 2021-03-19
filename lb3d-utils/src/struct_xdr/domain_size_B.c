#include <stdio.h>
#include <math.h>
#include "main.h"

int
domain_size(CASF *casf, Params *params, Stats *stats)
{
   int m;
   double fac;
   double sumS, sum_errS, sum_kerrS, sum_k2errS;

   casf->k1   = 0.;
   casf->k2   = 0.;
   sumS       = 0.;
   sum_errS   = 0.;
   sum_kerrS  = 0.;
   sum_k2errS = 0.;

   /* NORMALISE and ACCUM MOMENTA */
   for(m = 0 ; m < params->kradius; m++)
      if(casf->sum_shell[m] != 0) {

/*       Standard error of avg of {x_j} is the standard *       */
/*       deviation for the distribution of the average		*/
/*       (the latter is the best estimator for the mean, and	*/
/*       becomes a dirac-delta centred on the mean for N\to\infty), */
/*       and is \sigma_{N-1} / \sqrt{N}, where \sigma_{N-1} is the standard */
/*       deviation of the sample, and N is the number of elements   */
/*       of the sample. In other words:				*/

/*       \sqrt(							*/
/*       \frac{1}{N-1}\sum_j x_j^2 -				*/
/*       \frac{1}{N(N-1)}(\sum_j x_j)^2				*/
/*       / */
/*       N) */

/*       where instead of \sigma we've used \sigma_{N-1}, which the */
/*       sample							*/

	 fac         = DPI / params->nx * m;
	 casf->k1   += fac * (casf->S[m]);
	 casf->k2   += fac * fac * (casf->S[m]);
	 sumS       += casf->S[m];
	 sum_errS   += casf->errS[m];
	 sum_kerrS  += fac * casf->errS[m];
	 sum_k2errS += fac * fac * casf->errS[m];
     
      }
   
   casf->k1 /= sumS;
   casf->k2 /= sumS;
   casf->R1 = DPI / casf->k1;
   casf->R2 = DPI / sqrt(casf->k2);

   casf->errk1 = fabs(sum_kerrS - casf->k1 * sum_errS)/sumS;
   casf->errk2 = fabs(sum_k2errS - casf->k2 * sum_errS)/sumS;
   casf->errR1 = DPI / (casf->k1 * casf->k1) * casf->errk1;
   casf->errR2 = PI / pow(casf->k2, 1.5) * casf->errk2;

   return(1);

}



