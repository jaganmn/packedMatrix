## TODO/IDEAS:
## * we would really benefit from reusing some of the machinery
##   in R_MAIN_DIR/subscript.c to convert supplied subscripts to
##   "nice" in-bounds integer subscripts, which are necessary to
##   compute "triangular" subscripts
##
## * methods for x[i, j] and x[i, j, drop=] with neither i nor j missing ??
## * haven't given much thought to long vector support ...
## * ditto subassignment ...
## * an efficient 'apply' analogue for "packedMatrix" should be relatively
##   easy to implement now that row/column extraction are fast
##
## NTS:
## * need to use dimnames(x) instead of x@Dimnames to get symmetric
##   dimnames when dealing with "symmetricMatrix"


## UTILITIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.pM.error.oob <- function() {
    stop("subscript out of bounds")
}
.pM.error.ist <- function(i) {
    if (isS4(i)) {
        stop(sprintf("invalid subscript (S4) class '%s%'", class(i)))
    } else {
        stop(sprintf("invalid subscript type '%s%'", typeof(i)))
    }
}
.pM.error.dim <- function() {
    stop("incorrect number of dimensions")
}

## convert [i, j] integer index of triangular 'n'-by-'n' matrix 'x' to
## [k] integer index of x[upper.tri(x, TRUE)] or x[lower.tri(x, TRUE)]
.pM.arity21 <- function(i, j, n, up) {
    if (up) {
        i + ((j - 1L) * j) %/% 2L            # needs 1 <= i <= j <= n
    } else {
        i + ((j - 1L) * (2L * n - j)) %/% 2L # needs 1 <= j <= i <= n
    }
}

## do x[ij] efficiently where 'ij' is a two-column integer matrix
## containing only integers in '1:dim(x)[1L]' or NA
.pM.sub0.ij <- function(x, ij) {
    ## deal with [i, j] outside of ("opposite") stored triangle
    up <- x@uplo == "U"
    sy <- is(x, "symmetricMatrix") # otherwise "triangularMatrix"
    op <-
        if (up) {
            ij[, 1L] > ij[, 2L]
        } else {
            ij[, 1L] < ij[, 2L]
        }
    op <- !is.na(op) & op
    aop <- any(op)
    if (aop) {
        if (sy) {
            ## transpose
            ij[op, ] <- ij[op, 2:1, drop = FALSE]
        } else {
            ## discard ... we know the indexed element is zero
            ij <- ij[!op, , drop = FALSE]
        }
    }
    ## now subset
    k <- .pM.arity21(ij[, 1L], ij[, 2L], n = x@Dim[1L], up = up)
    xk <- x@x[k]
    if (aop && !sy) {
        replace(vector(typeof(xk), length(op)), !op, xk)
    } else {
        xk
    }
}

## NOTATION:
## utility .pM.sub<code>.<class> subsets "packedMatrix"
## code=
## 0: 1-ary indexing, as in x[i]
## 1: 2-ary indexing with one of 'i', 'j' missing, as in x[i, , drop=]
## 2: 2-ary indexing with neither 'i' nor 'j' missing, as in x[i, j, drop=]
## class=
##        index: vector index [includes factor, dispatches to chr,logi,num]
## chr,logi,num: vector index
##          mat: matrix, Matrix, or array index
##         null: NULL index

## Mainly to avoid defining methods for each "index" subclass
## separately, at the (reasonable?) cost of having to dispatch
## on index type ourselves ...
.pM.sub0.index <- function(x, i) {
    switch(mode(i),
           numeric =
               {
                   if (!is.numeric(i)) {
                       class(i) <- NULL
                   }
                   .pM.sub0.num(x, i)
               },
           logical = .pM.sub0.logi(x, i),
           character = .pM.sub0.chr(x, i),
           .pM.error.ist(i))
}
.pM.sub1.index <- function(x, i, drop, col) {
    switch(mode(i),
           numeric =
               {
                   if (!is.numeric(i)) {
                       class(i) <- NULL
                   }
                   .pM.sub1.num(x, i, drop = drop, col = col)
               },
           logical = .pM.sub1.logi(x, i, drop = drop, col = col),
           character = .pM.sub1.chr(x, i, drop = drop, col = col),
           .pM.error.ist(i))
}

## Emulating 'stringSubscript' in R_MAIN_DIR/subscript.c,
## though the first case, 'x[<character but not array>]',
## is quite pathological and really need not be supported ...
.pM.sub0.chr <- function(x, i) {
    rep.int(x@x[1L][NA], length(i))
}
.pM.sub1.chr <- function(x, i, drop, col) {
    if (length(i) > 0L) {
        nms <- dimnames(x)[[1L + col]]
        if (is.null(nms) || anyNA(i <- match(i, nms))) {
            .pM.error.oob()
        }
        .pM.sub1.num(x, i, drop = drop, col = col)
    } else {
        .pM.sub1.num(x, integer(0L), drop = drop, col = col)
    }
}

.pM.sub0.logi <- function(x, i) {
    ni <- length(i)
    if (ni == 0L) {
        return(x@x[0L])
    }
    n <- x@Dim[1L]
    ## optimize when identical(x@x, as(x, "[dln]geMatrix")@x)
    if (n <= 1L) {
        return(x@x[i])
    }
    if (anyNA(i)) {
        ## optimize x[NA], etc.
        if (all(is.na(i))) {
            return(rep.int(x@x[1L][NA], max(n * n, length(i))))
        }
    } else {
        ## optimize x[FALSE], etc.
        if (!any(i)) {
            return(x@x[0L])
        }
        ## optimize x[TRUE], etc.
        if (all(i)) {
            xi <- as(x, geClass(x))@x
            if (ni > n * n) {
                length(xi) <- ni
            }
            return(xi)
        }
    }
    ## dispatch
    ## FIXME: inefficient ... see TODO about adapting 'makeSubscript'
    ## though notably still an improvement over 'as(x, "matrix")[i]'
    ## which allocates a double vector of length n*n when 'x' is a "dMatrix"
    .pM.sub0.num(x, seq_len(n * n)[i])
}
.pM.sub1.logi <- function(x, i, drop, col) {
    p <- 1L + col # subset on this dimension
    pp <- 1L + !col # but not this one
    d <- x@Dim
    dn <- dimnames(x)
    empty <- function() {
        d[p] <- 0L
        dn[p] <- list(NULL)
        new(geClass(x), x = x@x[0L], Dim = d, Dimnames = dn)
    }
    ni <- length(i)
    if (ni == 0L) {
        return(empty())
    }
    n <- d[1L]
    if (ni > n) {
        stop("logical subscript too long")
    }
    if (anyNA(i)) {
        ## optimize x[NA, ], etc.
        if (all(is.na(i))) {
            x0 <- rep.int(x@x[1L][NA], n * n)
            if (!is.null(dn[[p]])) {
                dn[[p]][] <- NA
            }
            return(new(geClass(x), x = x0, Dim = d, Dimnames = dn))
        }
    } else {
        ## optimize x[FALSE, ], etc.
        if (!any(i)) {
            return(empty())
        }
        ## optimize x[TRUE, ], etc.
        if (all(i)) {
            if (n == 1L && !isFALSE(drop[1L])) {
                x0 <- x@x
                if (!is.null(dn[[pp]])) {
                    names(x0) <- dn[[pp]]
                }
                return(x0)
            } else {
                return(x)
            }
        }
    }
    ## dispatch
    .pM.sub1.num(x, seq_len(n)[i], drop = drop, col = col)
}

.pM.sub0.num <- function(x, i) {
    if (length(i) == 0L) {
        return(x@x[0L])
    }
    n <- x@Dim[1L]
    ## optimize when identical(x@x, as(x, "[dln]geMatrix")@x)
    if (n <= 1L) {
        return(x@x[i])
    }
    ## take care of nonpositive and noninteger values:
    ## after this, 'i' contains only elements in 'c(NA, 1:(n*n))'
    if (any(i < 0, na.rm = TRUE)) {
        ## FIXME: inefficient ... see TODO about adapting 'makeSubscript'
        i <- seq_len(n * n)[i]
    } else {
        if (is.double(i)) {
            i <- as.integer(i)
        }
        i <- i[i > 0L]
        i[i > n * n] <- NA
    }
    .Call(packedMatrix_sub0, x, i)
}

## A lot of acrobatics/complexity here ... should simplify but
## perhaps not worth the effort if implementing in C anyway ...
.pM.sub1.num <- function(x, i, drop, col) {
    d <- x@Dim
    n <- d[1L]
    if (any(i >= n + 1, na.rm = TRUE)) {
        .pM.error.oob()
    }
    ## take care of nonpositive and noninteger values:
    ## after this, 'i' contains only elements in 'c(NA, 1:n)'
    i <- seq_len(n)[i]
    if (length(i) == 0L) { # x[integer(0), , drop=]
        dn <- dimnames(x)
        if (n > 0L) {
            p <- 1L + col # subset on this dimension
            d[p] <- 0L
            dn[p] <- list(NULL)
        }
        return(new(geClass(x), x = x@x[0L], Dim = d, Dimnames = dn))
    }
    .Call(packedMatrix_sub1, x, i, drop, col)
}

## Could support "[dn]Matrix" and "array" ... leaving out for now
.pM.sub0.mat <- function(x, i) {
    ## if (is(i, "lMatrix") || is(i, "nMatrix") || is.logical(i)) {
    if (is(i, "lMatrix") || is.logical(i)) {
        return(.pM.sub0.logi(x, as.vector(i)))
    }
    d <- dim(i)
    if (length(d) == 2L && d[2L] == 2L) {
        ## if (is(i, "dMatrix")) {
        ##     i <- as(i, "matrix")
        ## }
        if (is.double(i)) {
            ## coerce to integer while preserving 'dim'
            storage.mode(i) <- "integer"
        }
        if (is.integer(i)) {
            ## rows containing 0 are deleted, rows containing NA result in NA,
            ## rows containing both are handled according to the first column
            i <- i[i[, 1L] != 0L, , drop = FALSE] # NA,j -> NA,NA
            i <- i[i[, 2L] != 0L, , drop = FALSE]
            if (dim(i)[1L] == 0L) {
                return(x@x[0L])
            }
            if (any(i < 1L, na.rm = TRUE)) {
                stop("negative values are not allowed in a matrix subscript")
            }
            if (any(i > x@Dim[1L], na.rm = TRUE)) {
                .pM.error.oob()
            }
            .pM.sub0.ij(x, i)
        } else if (is.character(i)) {
            if (d[1L] == 0L) {
                return(x@x[0L])
            }
            dn <- dimnames(x)
            m <- c(match(i[, 1L], dn[[1L]]), match(i[, 2L], dn[[2L]]))
            dim(m) <- d
            ## error if character row contains zero NA but integer row
            ## contains at least one NA, indicating nonmatch that should
            ## not be ignored
            if (any(rowSums(is.na(i)) == 0L & rowSums(is.na(m)) > 0L)) {
                .pM.error.oob()
            }
            .pM.sub0.ij(x, m)
        } else {
            .pM.error.ist(i)
        }
    } else {
        ## if (is(i, "dMatrix") || is.numeric(i)) {
        if (is.numeric(i)) {
            .pM.sub0.num(x, as.vector(i))
        } else if (is.character(i)) {
            .pM.sub0.chr(x, as.vector(i))
        } else {
            .pM.error.ist(i)
        }
    }
}

## Could support "[dln]Matrix" and "array" ... leaving out for now
.pM.sub1.mat <- function(x, i, drop, col) {
    ## if (is.numeric(i) || is(i, "dMatrix")) {
    if (is.numeric(i)) {
        .pM.sub1.num(x, as.vector(i), drop = drop, col = col)
    ## } else if (is.logical(i) || is(i, "lMatrix") || is(i, "nMatrix")) {
    } else if (is.logical(i)) {
        .pM.sub1.logi(x, as.vector(i), drop = drop, col = col)
    } else if (is.character(i)) {
        .pM.sub1.chr(x, as.vector(i), drop = drop, col = col)
    } else {
        .pM.error.ist(i)
    }
}

## unused argument 'i' included for consistency with other utilities
.pM.sub0.null <- function(x, i) {
    .pM.sub0.num(x, integer(0L))
}
.pM.sub1.null <- function(x, i, drop, col) {
    .pM.sub1.num(x, integer(0L), drop = drop, col = col)
}


## METHOD DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "missing"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[m, m, m] : nargs() = %d", na), .M.level = 2)
	      if (na == 2L) { # x[]
                  x
              } else if (na == 3L) { # x[, ]
                  ## drop=TRUE implicit
                  if (x@Dim[1L] == 1L) x@x else x
              } else { # x[, , ], and so on
                  .pM.error.dim()
              }
          })

setMethod("[", signature(x = "packedMatrix", i = "missing", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop) {
              na <- nargs()
              Matrix.msg(sprintf("pM[m, m, %s] : nargs() = %d", drop, na), .M.level = 2)
              if (na == 4L) { # x[, , drop=], x[, drop=, ], x[drop=, , ]
                  if (!isFALSE(drop[1L]) && x@Dim[1L] == 1L) x@x else x
              } else if (na < 4L) { # x[drop=], x[, drop=], x[drop=, ]
                  x
              } else { # x[, , , drop=], and so on
                  .pM.error.dim()
              }
          })

## Generalizing here to avoid repetition ...
.j.map <- c(index = "index", matrix = "mat", `NULL` = "null")
.i.map <- c(.j.map, lMatrix = "mat")
.op <- options(keep.source = FALSE)
for (.i.cl in names(.i.map)) {
    .nm0 <- as.name(paste0(".pM.sub0.", .i.map[[.i.cl]]))
    .nm1 <- as.name(paste0(".pM.sub1.", .i.map[[.i.cl]]))
    eval(bquote({
        setMethod("[", signature(x = "packedMatrix", i = .i.cl, j = "missing", drop = "missing"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[%s, m, m] : nargs() = %d", .(.i.cl), na), .M.level = 2)
                      if (na == 2L) { # x[i]
                          .(.nm0)(x, i)
                      } else if (na == 3L) { # x[i, ]
                          .(.nm1)(x, i, drop = TRUE, col = FALSE)
                      } else { # x[i, , ], etc.
                          .pM.error.dim()
                      }
                  })
        setMethod("[", signature(x = "packedMatrix", i = .i.cl, j = "missing", drop = "logical"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[%s, m, %s] : nargs() = %d", .(.i.cl), drop, na), .M.level = 2)
                      if (na == 3L) { # x[i, drop=]
                          .(.nm0)(x, i)
                      } else if (na == 4L) { # x[i, , drop=]
                          .(.nm1)(x, i, drop = drop, col = FALSE)
                      } else { # x[i, , , drop=], etc.
                          .pM.error.dim()
                      }
                  })
    })) # eval(bquote({
}
for (.j.cl in names(.j.map)) {
    .nm0 <- as.name(paste0(".pM.sub0.", .j.map[[.j.cl]]))
    .nm1 <- as.name(paste0(".pM.sub1.", .j.map[[.j.cl]]))
    eval(bquote({
        setMethod("[", signature(x = "packedMatrix", i = "missing", j = .j.cl, drop = "missing"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[m, %s, m] : nargs() = %d", .(.j.cl), na), .M.level = 2)
                      if (na == 2L) { # x[j=]
                          .(.nm0)(x, j)
                      } else if (na == 3L) { # x[, j]
                          .(.nm1)(x, j, drop = TRUE, col = TRUE)
                      } else { # x[, j, ], etc.
                          .pM.error.dim()
                      }
                  })
        setMethod("[", signature(x = "packedMatrix", i = "missing", j = .j.cl, drop = "logical"),
                  function(x, i, j, ..., drop) {
                      na <- nargs()
                      Matrix.msg(sprintf("pM[m, %s, %s] : nargs() = %d", .(.j.cl), drop, na), .M.level = 2)
                      if (na == 3L) { # x[j=, drop=]
                          .(.nm0)(x, j)
                      } else if (na == 4L) { # x[, j, drop=]
                          .(.nm1)(x, j, drop = drop, col = TRUE)
                      } else { # x[, j, , drop=], etc.
                          .pM.error.dim()
                      }
                  })
    })) # eval(bquote({
}
options(.op)
rm(.i.map, .j.map, .i.cl, .j.cl, .op, .nm0, .nm1)

setMethod("t", signature(x = "packedMatrix"),
          function(x) .Call(packedMatrix_t, x))
setMethod("diag", signature(x = "packedMatrix"),
          function(x, nrow, ncol, names) .Call(packedMatrix_diag_get, x, names))
setMethod("diag<-", signature(x = "packedMatrix"),
          function(x, value) .Call(packedMatrix_diag_set, x, value))
