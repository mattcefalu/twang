#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "ks.h"

/* .C calls */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef CEntries[]  = {
   {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
   CALLDEF(pSmirnov2x, 3),
   CALLDEF(pKS2, 2),
   {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
   {NULL, NULL, 0}
};

static const R_ExternalMethodDef ExtEntries[] = {
   {NULL, NULL, 0}
};

void attribute_visible R_init_stats(DllInfo *dll)
{
   R_registerRoutines(dll, CEntries, CallEntries, FortEntries, ExtEntries);
   R_useDynamicSymbols(dll, FALSE);
   R_forceSymbols(dll, TRUE);
}