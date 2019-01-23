/*
 *  psmirnov2x.h
 *  
 *  Copied from R 2.15 source
 *
 */

#include <R.h>
#include <Rmath.h>		/* constants */

#include "kspval.h"

/* Two-sided two-sample */
void
psmirnov2x(double *x, Sint *m, Sint *n)
{
    double md, nd, q, *u, w;
    Sint i, j;
	
    if(*m > *n) {
		i = *n; *n = *m; *m = i;
    }
    md = (double) (*m);
    nd = (double) (*n);
    /*
	 q has 0.5/mn added to ensure that rounding error doesn't
	 turn an equality into an inequality, eg abs(1/2-4/5)>3/10 
	 
	 */
    q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
    u = (double *) R_alloc(*n + 1, sizeof(double));
	
    for(j = 0; j <= *n; j++) {
		u[j] = ((j / nd) > q) ? 0 : 1;
    }
    for(i = 1; i <= *m; i++) {
		w = (double)(i) / ((double)(i + *n));
		if((i / md) > q)
			u[0] = 0;
		else
			u[0] = w * u[0];
		for(j = 1; j <= *n; j++) {
			if(fabs(i / md - j / nd) > q) 
				u[j] = 0;
			else
				u[j] = w * u[j] + u[j - 1];
		}
    }
    *x = u[*n];
}
