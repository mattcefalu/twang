/*
   *  ks.h
*  
   *
   *  Copied from R 3.6 source
*
   */
#include <R.h>
#include <Rinternals.h>

SEXP pSmirnov2x(SEXP statistic, SEXP snx, SEXP sny);
SEXP pKS2(SEXP sn, SEXP stol);