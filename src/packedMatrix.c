#include "packedMatrix.h"

SEXP packedMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) < 2) {
        return mkString(_("'Dim' slot has length less than two"));
    }
    int n = INTEGER(val)[0];
    if (INTEGER(val)[1] != n) {
        return mkString(_("Matrix is not square"));
    }
    val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym), "LU", "uplo");
    if (isString(val)) {
        return val;
    }
    val = GET_SLOT(obj, Matrix_xSym);
    if (LENGTH(val) != PACKED_LENGTH(n)) {
        return mkString(_("'x' slot has invalid length"));
    }
    return ScalarLogical(1);
}
