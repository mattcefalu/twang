/*
 * Copied from R 3.6 source
 
   *  R : A Computer Language for Statistical Data Analysis
*  Copyright (C) 1999-2016   The R Core Team.
*
   *  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
   *  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
   *  You should have received a copy of the GNU General Public License
*  along with this program; if not, a copy is available at
*  https://www.R-project.org/Licenses/
   */
   
   /* ks.c
Compute the asymptotic distribution of the one- and two-sample
two-sided Kolmogorov-Smirnov statistics, and the exact distributions
in the two-sided one-sample and two-sample cases.
*/
   
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>		/* constants */

/* Two-sample two-sided asymptotic distribution */
static void 
pkstwo(int n, double *x, double tol)
{
   /* x[1:n] is input and output
   *
      * Compute
   *   \sum_{k=-\infty}^\infty (-1)^k e^{-2 k^2 x^2}
   *   = 1 + 2 \sum_{k=1}^\infty (-1)^k e^{-2 k^2 x^2}
   *   = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^\infty \exp(-(2k-1)^2\pi^2/(8x^2))
   *
      * See e.g. J. Durbin (1973), Distribution Theory for Tests Based on the
   * Sample Distribution Function.  SIAM.
   *
      * The 'standard' series expansion obviously cannot be used close to 0;
   * we use the alternative series for x < 1, and a rather crude estimate
   * of the series remainder term in this case, in particular using that
   * ue^(-lu^2) \le e^(-lu^2 + u) \le e^(-(l-1)u^2 - u^2+u) \le e^(-(l-1))
   * provided that u and l are >= 1.
   *
      * (But note that for reasonable tolerances, one could simply take 0 as
         * the value for x < 0.2, and use the standard expansion otherwise.)
   *
      */
      double new, old, s, w, z;
   int i, k, k_max;
   
   k_max = (int) sqrt(2 - log(tol));
   
   for(i = 0; i < n; i++) {
      if(x[i] < 1) {
         z = - (M_PI_2 * M_PI_4) / (x[i] * x[i]);
         w = log(x[i]);
         s = 0;
         for(k = 1; k < k_max; k += 2) {
            s += exp(k * k * z - w);
         }
         x[i] = s / M_1_SQRT_2PI;
      }
      else {
         z = -2 * x[i] * x[i];
         s = -1;
         k = 1;
         old = 0;
         new = 1;
         while(fabs(old - new) > tol) {
            old = new;
            new += 2 * s * exp(z * k * k);
            s *= -1;
            k++;
         }
         x[i] = new;
      }
   }
}

/* Two-sided two-sample */
static double psmirnov2x(double *x, int m, int n)
{
   double md, nd, q, *u, w;
   int i, j;
   
   if(m > n) {
      i = n; n = m; m = i;
   }
   md = (double) m;
   nd = (double) n;
   /*
      q has 0.5/mn added to ensure that rounding error doesn't
   turn an equality into an inequality, eg abs(1/2-4/5)>3/10 
   */
   q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
   u = (double *) R_alloc(n + 1, sizeof(double));
   
   for(j = 0; j <= n; j++) {
   u[j] = ((j / nd) > q) ? 0 : 1;
   }
   for(i = 1; i <= m; i++) {
   w = (double)(i) / ((double)(i + n));
   if((i / md) > q)
   u[0] = 0;
   else
   u[0] = w * u[0];
   for(j = 1; j <= n; j++) {
   if(fabs(i / md - j / nd) > q) 
   u[j] = 0;
   else
   u[j] = w * u[j] + u[j - 1];
   }
   }
   return u[n];
}
   

/* Two-sided two-sample */
SEXP pSmirnov2x(SEXP statistic, SEXP snx, SEXP sny)
{
   int nx = asInteger(snx), ny = asInteger(sny);
   double st = asReal(statistic);
   return ScalarReal(psmirnov2x(&st, nx, ny));
   }

/* Two-sample two-sided asymptotic distribution */
SEXP pKS2(SEXP statistic, SEXP stol)
{
   int n = LENGTH(statistic);
   double tol = asReal(stol);
   SEXP ans = duplicate(statistic);
   pkstwo(n, REAL(ans), tol);
   return ans;
}

