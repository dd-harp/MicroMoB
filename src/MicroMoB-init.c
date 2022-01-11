#include "utils_sample.h"
#include <R_ext/Rdynload.h>

static const
  R_CallMethodDef callMethods[] = {
    {"C_draw_multinom", (DL_FUNC) &C_draw_multinom, 2},
    {NULL, NULL, 0}
  };

void R_init_MicroMoB(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
