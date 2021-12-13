#include "packedMatrix.h"

SEXP packedMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) < 2)
    {
        return mkString(_("'Dim' slot has length less than two"));
    }
    int n = INTEGER(val)[0];
    if (INTEGER(val)[1] != n)
    {
        return mkString(_("Matrix is not square"));
    }
    val = check_scalar_string(GET_SLOT(obj, Matrix_uploSym), "LU", "uplo");
    if (isString(val))
    {
        return val;
    }
    val = GET_SLOT(obj, Matrix_xSym);
    if (LENGTH(val) != PACKED_LENGTH(n))
    {
        return mkString(_("'x' slot has invalid length"));
    }
    return ScalarLogical(1);
}

SEXP packedMatrix_t(SEXP obj)
{
    SEXP val = PROTECT(duplicate(obj));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    const char *uplo = CHAR(STRING_ELT(GET_SLOT(obj, Matrix_uploSym), 0));
    int up = uplo[0] == 'U';

#define PM_T(_x0_, _x1_, _n_)						\
    int k = 0;								\
    /* Loop over [i,j] 0-indices of triangle of 'val' */		\
    if (up)								\
    {									\
	for (int j = 0; j < (_n_); ++j)					\
	{								\
	    for (int i = j; i < (_n_); ++i)				\
	    {								\
		(_x1_)[k++] = (_x0_)[j + (i * (i + 1)) / 2];		\
	    }								\
	}								\
    }									\
    else								\
    {									\
	int n2 = 2 * (_n_);						\
	for (int j = 0; j < (_n_); ++j)					\
	{								\
	    for (int i = 0; i <= j; ++i)				\
	    {								\
		(_x1_)[k++] = (_x0_)[j + (i * (n2 - i - 1)) / 2];	\
	    }								\
	}								\
    }									\

    if (n > 1)
    {
	SEXP x0 = GET_SLOT(obj, Matrix_xSym);
	SEXP x1 = GET_SLOT(val, Matrix_xSym);
	if (isReal(x0)) {
	    double *rx0 = REAL(x0);
	    double *rx1 = REAL(x1);
	    PM_T(rx0, rx1, n)
	} else {
	    int *lx0 = LOGICAL(x0);
	    int *lx1 = LOGICAL(x1);
	    PM_T(lx0, lx1, n)
	}
    }
#undef PM_T
    
    /* Toggle 'uplo' */
    SET_SLOT(val, Matrix_uploSym, mkString(up ? "L" : "U"));
    /* Reverse 'Dimnames' */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = GET_SLOT(val, Matrix_DimNamesSym);
    SET_VECTOR_ELT(dn1, 0, duplicate(VECTOR_ELT(dn0, 1)));
    SET_VECTOR_ELT(dn1, 1, duplicate(VECTOR_ELT(dn0, 0)));
    /* Reverse 'names(Dimnames)' if non-NULL */
    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (ndn0 != R_NilValue)
    {
	SEXP ndn1 = getAttrib(dn1, R_NamesSymbol);
	SET_STRING_ELT(ndn1, 0, STRING_ELT(ndn0, 1));
	SET_STRING_ELT(ndn1, 1, STRING_ELT(ndn0, 0));
    }
    UNPROTECT(1);
    return val;
}
