#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void psmirnov2x(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"psmirnov2x", (DL_FUNC) &psmirnov2x, 3},
    {NULL, NULL, 0}
};

void R_init_twang(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

