#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "000.api.h"

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static R_CallMethodDef callMethods[]  = {
  CALLDEF(sumsq, 1),
  CALLDEF(sumprod, 2),
  CALLDEF(negCDbeta, 3),
  {NULL, NULL, 0}
};

void R_init_hsa(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
