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
    Rboolean up = uplo_P(obj)[0] == 'U';

#define PM_T(_x0_, _x1_, _n_)						\
    int k = 0;								\
    /* Loop over [i,j] 0-indices of stored triangle of 'val' */		\
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
    
    /* Permute 'x' if not a no-op */
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
    if (!isNull(ndn0))
    {
	SEXP ndn1 = PROTECT(getAttrib(dn1, R_NamesSymbol));
	SET_STRING_ELT(ndn1, 0, STRING_ELT(ndn0, 1));
	SET_STRING_ELT(ndn1, 1, STRING_ELT(ndn0, 0));
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return val;
}

/* This utility improves existing *_getDiag machinery insofar 
   as it handles "names" precisely as documented in ?diag ... 
   
   On the other hand, it conditions on class (double/logical, 
   triangular/symmetric) which is very slightly wasteful ...
   [1 call to isReal, 1 call to R_has_slot]
   
   The existing *_getDiag machinery avoids this by special-casing 
   everything but seems much harder to maintain  ... ??
*/
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    SEXP val;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';
    
#define PM_D_G_UDIAG(_dest_, _one_)		\
    for (int j = 0; j < n; ++j)			\
    {						\
	(_dest_)[j] = _one_;			\
    }						\
    
#define PM_D_G_NDIAG(_dest_, _x_)					\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_dest_)[j] = (_x_)[pos];					\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
        {								\
	    (_dest_)[j] = (_x_)[pos];					\
        }								\
    }									\

    if (isReal(x)) { /* d[st]pMatrix */
	val = PROTECT(allocVector(REALSXP, n));
	double *rval = REAL(val);
	if (utr)
	{
	    PM_D_G_UDIAG(rval, 1.0);
	}
	else
	{
	    double *rx = REAL(x);
	    PM_D_G_NDIAG(rval, rx);
	}
    }
    else /* [ln][st]pMatrix */
    {
	val = PROTECT(allocVector(LGLSXP, n));
	int *lval = LOGICAL(val);
	if (utr)
	{
	    PM_D_G_UDIAG(lval, 1);
	}
	else
	{
	    int *lx = LOGICAL(x);
	    PM_D_G_NDIAG(lval, lx);
	}
    }

#undef PM_D_G_UDIAG
#undef PM_D_G_NDIAG

    /* Get 'names(val)' from 'dimnames(obj)' as documented in ?diag
       with symmetric dimnames as a special case
     */
    if (LOGICAL(nms)[0])
    {
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym);
	SEXP rn = VECTOR_ELT(dn, 0);
	SEXP cn = VECTOR_ELT(dn, 1);
	if (isNull(rn))
	{
	    if (!tr && !isNull(cn))
	    {
		setAttrib(val, R_NamesSymbol, cn);
	    }
	}
	else
	{
	    if (!tr || R_compute_identical(rn, cn, 16))
	    {
		setAttrib(val, R_NamesSymbol, rn);
	    }
	}
    }
    UNPROTECT(1);
    return val;
}
