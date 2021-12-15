#ifndef MATRIX_PAMATRIX_H
#define MATRIX_PAMATRIX_H

#include "Lapack-etc.h"
#include "Mutils.h"

SEXP packedMatrix_t(SEXP obj);
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms);
SEXP packedMatrix_diag_set(SEXP obj, SEXP val);

#endif
