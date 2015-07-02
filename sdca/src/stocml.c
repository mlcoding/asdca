//
//  stocml.c
//  
//
//  Created by Xingguo Li on 7/1/15.
//
//

#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#include "mymath.h"

SEXP ridge_sdca_call(SEXP Y_, SEXP X_, SEXP X_maxrn_, SEXP n_, SEXP d_, SEXP lambda_, SEXP nnlambda_, SEXP mmax_ite_, SEXP pprec_, SEXP idx_);


static R_CallMethodDef callMethods[] = {
    {"ridge_sdca_call", (DL_FUNC) &ridge_sdca_call, 10},
    {NULL, NULL, 0}
};

void R_init_ncvreg(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

