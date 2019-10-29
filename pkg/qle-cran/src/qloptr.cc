/**
 * @file qlopt.cc
 *
 * @date 24.06.2016
 * @author: M. Baaske
 */

#include "qsoptim.h"
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

//extern void (*expm)(double *x, int n, double *z, precond_type type);

/* R Interface functions  */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, (n)}

R_NativePrimitiveArgType argType_posdef[] = { REALSXP, INTSXP, INTSXP };
R_NativePrimitiveArgType argType_eigen[] = { INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP };

static R_CMethodDef CEntries[] = {
	{"isPositiveDefinite", (DL_FUNC) &isPositiveDefinite, 3, argType_posdef},
	{"geneigen", (DL_FUNC) &geneigen, 8, argType_eigen},
	{ NULL, NULL, 0, NULL}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(kriging,5),
    CALLDEF(estimateJacobian,5),
    CALLDEF(getDualKrigingWeights,3),
    CALLDEF(Fmat,2),
    CALLDEF(Pmat,1),
    CALLDEF(QSopt,5),
	CALLDEF(initQL,3),
    CALLDEF(finalizeQL,0),
    CALLDEF(qDValue,1),
	CALLDEF(mahalValue,1),
	CALLDEF(internQD,1),
	CALLDEF(internMD,1),
	CALLDEF(mahalanobis,6),
	CALLDEF(quasiDeviance,6),
    CALLDEF(covMatrix,2),
	CALLDEF(covValue,2),
//	CALLDEF(invertMatrix,2),
//	CALLDEF(RSolve,3),
    {NULL, NULL, 0}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_qle(DllInfo *info)
{
  R_registerRoutines(info, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


/*
 void R_unload_qle(DllInfo *info){
  /// Release resources.
}
*/


#ifdef  __cplusplus
}
#endif



