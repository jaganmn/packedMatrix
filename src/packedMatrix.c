#ifdef __GLIBC__
#define _POSIX_C_SOURCE 200809L
#endif
#include <string.h>
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

#define PM_AR21_UP(i, j) i + (j * (j + 1)) / 2
#define PM_AR21_LO(i, j, n2) i + (j * (n2 - j - 1)) / 2

SEXP packedMatrix_t(SEXP obj)
{
    /* Initialize result of same class */
    const char *cl = class_P(obj);
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl));
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? upper -> lower : lower -> upper */
    Rboolean up = uplo_P(obj)[0] == 'U';

#define PM_T_LOOP(_px0_, _px1_)						\
    /* Loop over [i,j] 0-indices of stored triangle of _result_ */	\
    if (up)								\
    {									\
	for (int j = 0, k = 0; j < n; ++j)				\
	{								\
	    for (int i = j; i < n; ++i)					\
	    {								\
		(_px1_)[k++] = (_px0_)[PM_AR21_UP(j, i)]; \
	    }								\
	}								\
    }									\
    else								\
    {									\
	int n2 = 2 * n;							\
	for (int j = 0, k = 0; j < n; ++j)				\
	{								\
	    for (int i = 0; i <= j; ++i)				\
	    {								\
		(_px1_)[k++] = (_px0_)[PM_AR21_LO(j, i, n2)]; \
	    }								\
	}								\
    }

#define PM_T(_datatype_, _sexptype_, _accessor_)			\
    SEXP x1 = PROTECT(allocVector(_sexptype_, LENGTH(x0)));		\
    _datatype_ *px0 = _accessor_(x0);					\
    _datatype_ *px1 = _accessor_(x1);					\
    PM_T_LOOP(px0, px1);						\
    SET_SLOT(res, Matrix_xSym, x1);					\
    UNPROTECT(1)

    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    if (n > 1)
    {
	/* Permute 'x' slot */
	if (isReal(x0))
	{
	    PM_T(double, REALSXP, REAL);
	}
	else
	{
	    PM_T(int, LGLSXP, LOGICAL);
	}
    }
    else
    {
	/* Preserve 'x' slot */
	SET_SLOT(res, Matrix_xSym, x0);
    }

#undef PM_T
#undef PM_T_LOOP

    /* Toggle 'uplo' slot */
    SET_SLOT(res, Matrix_uploSym, mkString(up ? "L" : "U"));
    /* Preserve 'Dim' slot */
    SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
    /* Reverse 'Dimnames' slot and (if not absent) 'names(Dimnames)' */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn1, 0, VECTOR_ELT(dn0, 1));
    SET_VECTOR_ELT(dn1, 1, VECTOR_ELT(dn0, 0));
    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (!isNull(ndn0))
    {
	SEXP ndn1 = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(ndn1, 0, STRING_ELT(ndn0, 1));
	SET_STRING_ELT(ndn1, 1, STRING_ELT(ndn0, 0));
	setAttrib(dn1, R_NamesSymbol, ndn1);
	UNPROTECT(1);
    }
    SET_SLOT(res, Matrix_DimNamesSym, dn1);
    UNPROTECT(2);
    return res;
} /* packedMatrix_t */

/* 'packedMatrix_diag_get' handles 'names' precisely as documented
   in '?diag', whereas current '*_getDiag' machinery neglects 'names'.

   It combines the work of 6 existing C utilities:
   '[dl][st]pMatrix_detDiag' and '[dl]_packed_getDiag'.

   It does so at the cost of conditioning on class
   (double/logical, triangular/symmetric).
   However, this conditioning happens outside of any loop,
   so there is minimal cost to performance
   (1 call to isReal, 1 call to R_has_slot) ...
*/
SEXP packedMatrix_diag_get(SEXP obj, SEXP nms)
{
    int do_nms = asLogical(nms);
    if (do_nms == NA_LOGICAL)
    {
	error(_("'names' must be TRUE or FALSE"));
    }

    SEXP res;
    SEXP x = GET_SLOT(obj, Matrix_xSym);
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';

#define PM_D_G_UDIAG(_pres_, _one_)		\
    for (int j = 0; j < n; ++j)			\
    {						\
	(_pres_)[j] = _one_;			\
    }

#define PM_D_G_NDIAG(_pres_, _px_)					\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_pres_)[j] = (_px_)[pos];					\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
        {								\
	    (_pres_)[j] = (_px_)[pos];					\
        }								\
    }

#define PM_D_G(_datatype_, _sexptype_, _accessor_, _one_)	\
    res = PROTECT(allocVector(_sexptype_, n));			\
    _datatype_ *pres = _accessor_(res);				\
    if (utr)							\
    {								\
	PM_D_G_UDIAG(pres, _one_);				\
    }								\
    else							\
    {								\
	_datatype_ *px = _accessor_(x);				\
	PM_D_G_NDIAG(pres, px);					\
    }

    if (isReal(x)) { /* d[st]pMatrix */
	PM_D_G(double, REALSXP, REAL, 1.0);
    }
    else /* [ln][st]pMatrix */
    {
	PM_D_G(int, LGLSXP, LOGICAL, 1);
    }

#undef PM_D_G
#undef PM_D_G_NDIAG
#undef PM_D_G_UDIAG

    if (do_nms)
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


/* 'packedMatrix_diag_set' handles coercions similarly to base R's
   'diag<-' and performs explicit tests for type compatibility,
   whereas current '*_setDiag' machinery does neither. It also does
   more to avoid memory allocation ...

   It combines the work of 8 existing C utilities:
   '[dl][st]pMatrix_setDiag' and '(tr_)?[dl]_packed_setDiag'

   Like 'packedMatrix_diag_get', it conditions on class _outside_
   of any loop ...
*/
SEXP packedMatrix_diag_set(SEXP obj, SEXP val)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int nv = LENGTH(val);
    Rboolean nv1 = (nv == 1);
    if (!(nv1 || nv == n))
    {
	error(_("replacement diagonal has wrong length"));
    }

    SEXP res;
    int nprotect = 0;
    if (MAYBE_REFERENCED(obj)) /* MAYBE_SHARED seems less safe ... */
    {
	/* Avoids allocation where possible */
	const char *cl = class_P(obj);
	res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
	SET_SLOT(res, Matrix_DimSym, GET_SLOT(obj, Matrix_DimSym));
	SET_SLOT(res, Matrix_DimNamesSym, GET_SLOT(obj, Matrix_DimNamesSym));
	SET_SLOT(res, Matrix_uploSym, GET_SLOT(obj, Matrix_uploSym));
	slot_dup(res, obj, Matrix_xSym);
    }
    else
    {
	/* Avoids allocation _entirely_ */
	res = obj;
    }
    if (Diag_P(res)[0] == 'U')
    {
	/* Assigning to unit diagonal "[dln]tpMatrix" toggles 'diag' slot */
	SET_SLOT(res, Matrix_diagSym, mkString("N"));
    }

    SEXP x = GET_SLOT(res, Matrix_xSym);
    /* ? double result : logical result */
    Rboolean dbl = TRUE;
    /* ? upper : lower */
    Rboolean up = uplo_P(res)[0] == 'U';

    /* Test that LHS and RHS of assignment have compatible types,
       coercing one or the other from logical to double if necessary
     */
    switch (TYPEOF(x))
    {
    case LGLSXP:
	switch (TYPEOF(val))
	{
	case LGLSXP:
	    dbl = FALSE;
	    break;
	case INTSXP:
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
	case REALSXP:
	{
	    /* [ln][st]pMatrix -> d[st]pMatrix */
	    SEXP strcl = getAttrib(res, R_ClassSymbol);
	    char *cl = strdup(CHAR(STRING_ELT(strcl, 0)));
	    cl[0] = 'd';
	    SET_STRING_ELT(strcl, 0, mkChar(cl));
	    free(cl);
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
	    /* logical, integer -> double */
	    val = PROTECT(coerceVector(val, REALSXP)); ++nprotect;
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

#define PM_D_S_ONE(_px_, _d_)						\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_px_)[pos] = _d_;						\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
        {								\
	    (_px_)[pos] = _d_;						\
        }								\
    }

#define PM_D_S_FULL(_px_, _pval_)					\
    if (up)								\
    {									\
        for (int j = 0, pos = 0; j < n; pos += 1 + (++j))		\
	{								\
	    (_px_)[pos] = (_pval_)[j];					\
	}								\
    }									\
    else								\
    {									\
	for (int j = 0, pos = 0; j < n; pos += n - (j++))		\
	{								\
	    (_px_)[pos] = (_pval_)[j];					\
	}								\
    }

#define PM_D_S(_datatype_, _accessor_)		\
    _datatype_ *px = _accessor_(x);		\
    _datatype_ *pval = _accessor_(val);		\
    if (nv1)					\
    {						\
	_datatype_ d = pval[0];			\
	PM_D_S_ONE(px, d);			\
    }						\
    else					\
    {						\
	PM_D_S_FULL(px, pval);			\
    }

    if (dbl) /* d[st]pMatrix */
    {
	PM_D_S(double, REAL);
    }
    else /* [ln][st]pMatrix */
    {
	PM_D_S(int, LOGICAL);
    }

#undef PM_D_S
#undef PM_D_S_FULL
#undef PM_D_S_ONE

    UNPROTECT(nprotect);
    return res;
} /* packedMatrix_diag_set */


/* Low level implementation of 'x[i]' in the "nice" situation 
   in which 'i' supplies only in-bounds integers and NA ...
*/
SEXP packedMatrix_sub0(SEXP obj, SEXP index)
{
    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    SEXP x = GET_SLOT(obj, Matrix_xSym);

    int nindex = LENGTH(index);
    SEXP res = PROTECT(allocVector(TYPEOF(x), nindex));

    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';

#define PM_SUB0_LOOP(_na_, _zero_, _one_)				\
    if (up)								\
    {									\
	for (int pos = 0, i, j, k; pos < nindex; ++pos)			\
	{								\
	    k = pindex[pos];						\
	    if (k == NA_INTEGER)					\
	    {								\
		pres[pos] = _na_;					\
									\
	    }								\
            else							\
	    {								\
		k -= 1;	/* 1-index to 0-index */			\
		i = k % n;						\
		j = k / n;						\
		if (i > j)						\
		{							\
		    if (tr)						\
		    {							\
			pres[pos] = _zero_;				\
		    }							\
		    else						\
		    {							\
			pres[pos] = px[PM_AR21_UP(j, i)];		\
		    }							\
		}							\
		else if (utr && i == j)					\
		{							\
		    pres[pos] = _one_;					\
		}							\
		else							\
		{							\
		    pres[pos] = px[PM_AR21_UP(i, j)];			\
		}							\
	    }								\
	}								\
    }									\
    else								\
    {									\
	int n2 = 2 * n;							\
	for (int k = 0, i, j, pos; k < nindex; ++k)			\
	{								\
	    pos = pindex[k];						\
	    if (pos == NA_INTEGER)					\
	    {								\
		pres[k] = _na_;						\
	    }								\
	    else							\
	    {								\
		pos -= 1;						\
		i = pos % n;						\
		j = pos / n;						\
		if (tr && i < j)					\
		{							\
		    if (tr)						\
		    {							\
			pres[k] = _zero_;				\
		    }							\
		    else						\
		    {							\
			pres[k] = px[PM_AR21_LO(j, i, n2)];		\
		    }							\
		}							\
		else if (utr && i == j)					\
		{							\
		    pres[k] = _one_;					\
		}							\
		else							\
		{							\
		    pres[k] = px[PM_AR21_LO(i, j, n2)];			\
		}							\
	    }								\
	}								\
    }

    int *pindex = INTEGER(index);
    if (isReal(x))
    {
	double *px = REAL(x);
	double *pres = REAL(res);
	PM_SUB0_LOOP(NA_REAL, 0.0, 1.0);
    }
    else
    {
	int *px = LOGICAL(x);
	int *pres = LOGICAL(res);
	PM_SUB0_LOOP(NA_LOGICAL, 0, 1);
    }

#undef PM_SUB0

    UNPROTECT(1);
    return res;
} /* packedMatrix_sub0 */

/* Low level implementation of 'x[i, ]' and 'x[, i]' in the "nice" situation 
   in which 'i' supplies only in-bounds integers and NA ...
*/
SEXP packedMatrix_sub1(SEXP obj, SEXP index, SEXP drop, SEXP col)
{
    int do_drop = asLogical(drop);
    if (do_drop == NA_LOGICAL)
    {
	error(_("'drop' must be TRUE or FALSE"));
    }

    int n = INTEGER(GET_SLOT(obj, Matrix_DimSym))[0];
    int n2 = 2 * n;
    int nindex = LENGTH(index);
    int nprotect = 0;

    /* ? triangular : symmetric */
    Rboolean tr = R_has_slot(obj, Matrix_diagSym);
    /* ? unit triangular : other */
    Rboolean utr = tr && diag_P(obj)[0] == 'U';
    /* ? upper : lower */
    Rboolean up = uplo_P(obj)[0] == 'U';
    /* ? index columns : index rows */
    Rboolean do_col = LOGICAL(col)[0];
    /* redundant but slightly easier to conceptualize */
    int mar = !do_col;

    /* Initialize result of same type but "general" class */
    char *cl = strdup(class_P(obj));
    cl[1] = 'g';
    cl[2] = 'e';
    SEXP res = PROTECT(NEW_OBJECT_OF_CLASS(cl)); ++nprotect;
    free(cl);

    /* Set 'x' slot */
    SEXP x0 = GET_SLOT(obj, Matrix_xSym);
    SEXP x1 = PROTECT(allocVector(TYPEOF(x0), n * nindex)); ++nprotect;
    SET_SLOT(res, Matrix_xSym, x1);

    /* Set 'Dim' slot */
    SEXP d1 = PROTECT(GET_SLOT(res, Matrix_DimSym));
    INTEGER(d1)[mar] = n;
    INTEGER(d1)[!mar] = nindex;
    UNPROTECT(1);

    /* Set 'Dimnames' slot and (if not absent) 'names(Dimnames)' ...
       _could_ use 'symmetric_DimNames' but it copies unnecessarily ...
     */
    SEXP dn0 = GET_SLOT(obj, Matrix_DimNamesSym);
    SEXP dn1 = GET_SLOT(res, Matrix_DimNamesSym);

    SEXP ndn0 = getAttrib(dn0, R_NamesSymbol);
    if (!isNull(ndn0))
    {
	SEXP ndn1 = PROTECT(duplicate(ndn0));
	if (!tr)
	{
	    /* Enforce symmetry of 'names(Dimnames)' */
	    if (strcmp(CHAR(STRING_ELT(ndn1, 0)), CHAR(STRING_ELT(ndn1, 1))))
	    {
		int nonemp = LENGTH(STRING_ELT(ndn1, 0)) == 0;
		SET_STRING_ELT(ndn1, !nonemp, STRING_ELT(ndn1, nonemp));
	    }
	}
	setAttrib(dn1, R_NamesSymbol, ndn1);
	UNPROTECT(1);
    }

    SEXP srcnms, destnms;
    srcnms = VECTOR_ELT(dn0, !mar);
    Rboolean do_nms = !isNull(srcnms);
    if (tr)
    {
	SET_VECTOR_ELT(dn1, mar, VECTOR_ELT(dn0, mar));
    }
    else
    {
	/* Enforce symmetry of 'Dimnames' */
	if (!do_nms)
	{
	    srcnms = VECTOR_ELT(dn0, mar);
	    do_nms = !isNull(srcnms);
	}
	SET_VECTOR_ELT(dn1, mar, srcnms);
    }
    if (nindex == 0)
    {
	do_nms = FALSE;
    }
    if (do_nms)
    {
	destnms = PROTECT(allocVector(STRSXP, nindex)); nprotect++;
	SET_VECTOR_ELT(dn1, !mar, destnms);
    }
    else
    {
	SET_VECTOR_ELT(dn1, !mar, R_NilValue);
    }

    /* FIXME: condense or separate different cases to reduce repetition
       and avoid conditioning inside of loops
     */
#define PM_SUB1_LOOP(_na_, _zero_, _one_)				\
    int incr = (do_col ? 1 : nindex);					\
    for (int i, k = 0, pos = 0; k < nindex; ++k)			\
    {									\
	i = pindex[k];							\
	if (i == NA_INTEGER)						\
	{								\
	    if (do_nms)							\
	    {								\
		SET_STRING_ELT(destnms, k, NA_STRING);			\
	    }								\
	    for (int j = 0; j < n; ++j, pos += incr)			\
	    {								\
		px1[pos] = _na_;					\
	    }								\
	}								\
	else								\
	{								\
	    i -= 1; /* 1-index -> 0-index */				\
	    if (do_nms)							\
	    {								\
		SET_STRING_ELT(destnms, k, STRING_ELT(srcnms, i));	\
	    }								\
	    if (tr)							\
	    {								\
		if (up)							\
		{							\
		    if (do_col)						\
		    {							\
			for (int j = 0; j < i; ++j, pos += incr)	\
			{						\
			    px1[pos] = px0[PM_AR21_UP(i, j)];		\
			}						\
			for (int j = i + 1; j < n; ++j, pos += incr)	\
			{						\
			    px1[pos] = _zero_;				\
			}						\
		    }							\
		    else						\
		    {							\
			for (int j = 0; j < i; ++j, pos += incr)	\
			{						\
			    px1[pos] = _zero_;				\
			}						\
			for (int j = i; j < n; ++j, pos += incr)	\
			{						\
			    px1[pos] = px0[PM_AR21_UP(i, j)];		\
			}						\
		    }							\
		}							\
		else							\
		{							\
		    if (do_col)						\
		    {							\
			for (int j = 0; j < i; ++j, pos += incr)	\
			{						\
			    px1[pos] = _zero_;				\
			}						\
			for (int j = i; j < n; ++j, pos += incr)	\
			{						\
			    px1[pos] = px0[PM_AR21_LO(i, j, n2)];	\
			}						\
		    }							\
		    else						\
		    {							\
			for (int j = 0; j < i + 1; ++j, pos += incr)	\
			{						\
			    px1[pos] = px0[PM_AR21_LO(i, j, n2)];	\
			}						\
			for (int j = i + 1; j < n; ++j, pos += incr)	\
			{						\
			    px1[pos] = _zero_;				\
			}						\
		    }							\
		}							\
		if (utr)						\
		{							\
		    px1[pos - (n - i) * incr] = _one_;			\
		}							\
	    }								\
	    else							\
	    {								\
		if (up)							\
		{							\
		    for (int j = 0; j < i; ++j, pos += incr)		\
		    {							\
		        px1[pos] = px0[PM_AR21_UP(j, i)];		\
		    }							\
		    for (int j = i; j < n; ++j, pos += incr)		\
		    {							\
			px1[pos] = px0[PM_AR21_UP(i, j)];		\
		    }							\
		}							\
		else							\
		{							\
		    for (int j = 0; j < i + 1; ++j, pos += incr)	\
		    {							\
			px1[pos] = px0[PM_AR21_LO(i, j, n2)];		\
		    }							\
		    for (int j = i + 1; j < n; ++j, pos += incr)	\
		    {							\
			px1[pos] = px0[PM_AR21_LO(j, i, n2)];		\
		    }							\
		}							\
	    }								\
	}								\
	if (!do_col)							\
	{								\
	    pos = (pos + 1) % nindex;					\
	}								\
    }

    int *pindex = INTEGER(index);
    if (isReal(x0))
    {
	double *px0 = REAL(x0);
	double *px1 = REAL(x1);
	PM_SUB1_LOOP(NA_REAL, 0.0, 1.0);
    }
    else
    {
	int *px0 = LOGICAL(x0);
	int *px1 = LOGICAL(x1);
	PM_SUB1_LOOP(NA_LOGICAL, 0, 1);
    }

#undef PM_SUB1_LOOP

    /* Drop dimensions in this special case */
    if (nindex == 1 && do_drop)
    {
	SEXP nms = VECTOR_ELT(dn1, mar);
	if (!isNull(nms))
	{
	    setAttrib(x1, R_NamesSymbol, nms);
	}
	res = x1;
    }
    UNPROTECT(nprotect);
    return res;
} /* packedMatrix_sub1 */

#undef PM_AR21_UP
#undef PM_AR21_LO
