#include "packedMatrix.h"

SEXP packedMatrix_validate(SEXP obj)
{
    SEXP val = GET_SLOT(obj, Matrix_DimSym);
    if (LENGTH(val) != 2)
    {
        return mkString(_("'Dim' slot does not have length 2"));
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
    if (R_has_slot(obj, Matrix_diagSym))
    {
	val = check_scalar_string(GET_SLOT(obj, Matrix_diagSym), "NU", "diag");
	if (isString(val))
	{
	    return val;
	}
    }
    val = GET_SLOT(obj, Matrix_xSym);
    if (LENGTH(val) != PACKED_LENGTH(n))
    {
        return mkString(_("'x' slot does not have length 'Dim[1]*(Dim[1]+1)/2'"));
    }
    return ScalarLogical(1);
}

SEXP packedMatrix_t(SEXP obj)
{
    /* Initialize result of same class */
    const char *cl = class_P(obj);
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? upper->lower : lower->upper */
    Rboolean up = uplo_P(obj)[0] == 'U';
    
#define PM_T(_x0_, _x1_)						\
    int k = 0;								\
    /* Loop over [i,j] 0-indices of stored triangle of result */	\
    if (up)								\
    {									\
	for (int j = 0; j < n; ++j)					\
	{								\
	    for (int i = j; i < n; ++i)					\
	    {								\
		(_x1_)[k++] = (_x0_)[j + (i * (i + 1)) / 2];		\
	    }								\
	}								\
    }									\
    else								\
    {									\
	int n2 = 2 * n;							\
	for (int j = 0; j < n; ++j)					\
	{								\
	    for (int i = 0; i <= j; ++i)				\
	    {								\
		(_x1_)[k++] = (_x0_)[j + (i * (n2 - i - 1)) / 2];	\
	    }								\
	}								\
    }									\
    
    if (n > 1)
    {
	/* Permute 'x' */
	SEXP x0 = GET_SLOT(obj, Matrix_xSym);
	SEXP x1;
	if (isReal(x0)) {
	    x1 = PROTECT(allocVector(REALSXP, LENGTH(x0)));
	    double *rx0 = REAL(x0);
	    double *rx1 = REAL(x1);
	    PM_T(rx0, rx1)
	} else {
	    x1 = PROTECT(allocVector(LGLSXP,  LENGTH(x0)));
	    int *lx0 = LOGICAL(x0);
	    int *lx1 = LOGICAL(x1);
	    PM_T(lx0, lx1)
	}
	SET_SLOT(res, Matrix_xSym, x1);
	UNPROTECT(1);
    }
    else
    {
	/* Copy 'x' */
	slot_dup(res, obj, Matrix_xSym);
    }
    
#undef PM_T

    /* Toggle 'uplo' */
    SET_SLOT(res, Matrix_uploSym, mkString(up ? "L" : "U"));
    /* Copy 'Dim' */
    slot_dup(res, obj, Matrix_DimSym);
    /* Reverse 'Dimnames' and (if not absent) 'names(Dimnames)' */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn1, 0, duplicate(VECTOR_ELT(dn0, 1)));
    SET_VECTOR_ELT(dn1, 1, duplicate(VECTOR_ELT(dn0, 0)));
    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (!isNull(ndn0))
    {
	SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(ndn1, 0, duplicate(STRING_ELT(ndn0, 1)));
	SET_STRING_ELT(ndn1, 1, duplicate(STRING_ELT(ndn0, 0)));
	setAttrib(dn1, R_NamesSymbol, ndn1);
	UNPROTECT(1);
    }
    SET_SLOT(res, Matrix_DimNamesSym, dn1);
    UNPROTECT(2);
    return res;
} /* packedMatrix_t */

/* This utility improves existing *_getDiag machinery insofar 
   as it handles 'names' precisely as documented in ?diag ... 
   
   On the other hand, it conditions on class (double/logical, 
   triangular/symmetric) which is very slightly wasteful ...
   [1 call to isReal, 1 call to R_has_slot]
   
   Existing *_getDiag machinery avoids this by special-casing 
   everything but seems much harder to maintain ... ??
*/
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int nms1 = asLogical(nms);
    if (nms1 == NA_LOGICAL)
    {
	error(_("'names' must be TRUE or FALSE"));
    }
	
    SEXP res;
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

    SEXP x = GET_SLOT(obj, Matrix_xSym);
    if (isReal(x)) { /* d[st]pMatrix */
	res = PROTECT(allocVector(REALSXP, n));
	double *rres = REAL(res);
	if (utr)
	{
	    PM_D_G_UDIAG(rres, 1.0);
	}
	else
	{
	    double *rx = REAL(x);
	    PM_D_G_NDIAG(rres, rx);
	}
    }
    else /* [ln][st]pMatrix */
    {
	res = PROTECT(allocVector(LGLSXP, n));
	int *lres = LOGICAL(res);
	if (utr)
	{
	    PM_D_G_UDIAG(lres, 1);
	}
	else
	{
	    int *lx = LOGICAL(x);
	    PM_D_G_NDIAG(lres, lx);
	}
    }

#undef PM_D_G_UDIAG
#undef PM_D_G_NDIAG

    if (nms1)
    {
	SEXP dn = GET_SLOT(obj, Matrix_DimNamesSym);
	SEXP rn = VECTOR_ELT(dn, 0);
	SEXP cn = VECTOR_ELT(dn, 1);
	if (isNull(rn))
	{
	    if (!tr && !isNull(cn))
	    {
		setAttrib(res, R_NamesSymbol, cn);
	    }
	}
	else
	{
	    if (!tr || R_compute_identical(rn, cn, 16))
	    {
		setAttrib(res, R_NamesSymbol, rn);
	    }
	}
    }
    UNPROTECT(1);
    return res;
} /* packedMatrix_diag_get */

SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int nv = LENGTH(val);
    Rboolean nv1 = (nv == 1);
    if (!(nv1 || nv == n))
    {
	error(_("replacement diagonal has wrong length"));
    }

    /* PROTECT counter */
    int ptc = 0;
    /* ? double result : logical result */
    Rboolean dbl = TRUE;
    
    SEXP res = PROTECT(duplicate(obj)); ++ptc;
    SEXP x = GET_SLOT(res, Matrix_xSym);
    switch (TYPEOF(x))
    {
    case LGLSXP:
	switch (TYPEOF(val))
	{
	case LGLSXP:
	    dbl = FALSE;
	    break;
	case INTSXP:
	    val = PROTECT(coerceVector(val, REALSXP)); ++ptc;
	case REALSXP:
	{
	    /* [ln][st]pMatrix->d[st]pMatrix */
	    char *cl = class_P(res);
	    cl[0] = 'd';
	    SET_SLOT(res, Matrix_xSym, coerceVector(x, REALSXP));
	    x = GET_SLOT(res, Matrix_xSym);
	    break;
	}
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    case REALSXP:
	switch (TYPEOF(val))
	{
	case LGLSXP:
	case INTSXP:
	    val = PROTECT(coerceVector(val, REALSXP)); ++ptc;
	case REALSXP:
	    break;
	default:
	    error(_("replacement diagonal has incompatible type '%s'"),
		  type2char(TYPEOF(val)));
	}
	break;
    default:
	error(_("'x' slot is not of type '%s' or '%s', which should never happen; please report"),
	      type2char(LGLSXP), type2char(REALSXP));
    }

    /* ? upper : lower */
    Rboolean up = uplo_P(res)[0] == 'U';    

#define PM_D_S_ONE(_x_, _val_)						\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_x_)[pos] = _val_;						\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
        {								\
	    (_x_)[pos] = _val_;						\
        }								\
    }									\
    
#define PM_D_S_FULL(_x_, _val_)						\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_x_)[pos] = (_val_)[j];					\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
	{								\
	    (_x_)[pos] = (_val_)[j];					\
	}								\
    }									\
    
    if (dbl) /* d[st]pMatrix */
    {
	double *rx = REAL(x);
	double *rval = REAL(val);
	if (nv1)
	{
	    PM_D_S_ONE(rx, rval[0])
	}
	else
	{
	    PM_D_S_FULL(rx, rval)
	}
    }
    else /* [ln][st]pMatrix */
    {
	int *lx = LOGICAL(x);
	int *lval = LOGICAL(val);
	if (nv1)
	{
	    PM_D_S_ONE(lx, lval[0])
	}
	else
	{
	    PM_D_S_FULL(lx, lval)
	}
    }

#undef PM_D_S_ONE
#undef PM_D_S_FULL
    
    /* Assigning to unit triangular Matrix toggles 'diag' */
    if (Diag_P(res)[0] == 'U')
    {
	SET_SLOT(res, Matrix_diagSym, mkString("N"));
    }
    UNPROTECT(ptc);
    return res;
} /* packedMatrix_diag_set */
