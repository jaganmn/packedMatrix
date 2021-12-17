packedMatrix
================
2021-12-17

Based on discussion on R-devel beginning
[here](https://stat.ethz.ch/pipermail/r-devel/2021-November/081261.html)
and continued privately with Martin Mächler.

## Viewing changes relative to r3419

``` bash
$ git diff --compact-summary 2f64f5f HEAD
```

## Class definition

``` r
library("Matrix")
showClass("packedMatrix")
```

    ## Virtual Class "packedMatrix" [package "Matrix"]
    ## 
    ## Slots:
    ##                                     
    ## Name:       uplo       Dim  Dimnames
    ## Class: character   integer      list
    ## 
    ## Extends: 
    ## Class "Matrix", directly
    ## Class "mMatrix", by class "Matrix", distance 2
    ## Class "replValueSp", by class "Matrix", distance 2
    ## 
    ## Known Subclasses: 
    ## Class "dtpMatrix", directly
    ## Class "dspMatrix", directly
    ## Class "ltpMatrix", directly
    ## Class "lspMatrix", directly
    ## Class "ntpMatrix", directly
    ## Class "nspMatrix", directly
    ## Class "dppMatrix", by class "dspMatrix", distance 2
    ## Class "pCholesky", by class "dtpMatrix", distance 2
    ## Class "pBunchKaufman", by class "dtpMatrix", distance 2

## Utilities for tests

`checkEquiv` takes as an argument a list of calls `FUN(x, ...)` and
returns a logical vector of the same length, indicating whether the
values of `FUN(x, ...)` are consistent with the values of
`FUN(as(x, "matrix"), ...)`. (Here, *consistent* means identical modulo
a coercion from `"Matrix"` to `"matrix"`.)

``` r
checkEquiv <- function(exprs) {
  f <- function(x1, envir) {
    stopifnot(is.call(x1), length(x1) > 1L)
    x2 <- x1
    x2[[2L]] <- call("as", x2[[2L]], "matrix") 
    m1 <- eval(x1, envir)
    m2 <- eval(x2, envir)
    if (is(m1, "Matrix")) {
      m1 <- as(m1, "matrix")
    }
    identical(m1, m2)
  }
  res <- vapply(exprs, f, NA, parent.frame())
  names(res) <- vapply(exprs, deparse, "")
  res
}
```

`checkTimes` takes as an argument a list of calls `FUN(x, ...)` and
returns a benchmark comparing their evaluation times to those of the
calls `FUN(as(x, "matrix"), ...)`.

``` r
checkTimes <- function(exprs, times = 100L) {
  stopifnot(requireNamespace("microbenchmark"))
  f <- function(x) {
    stopifnot(is.call(x), length(x) > 1L)
    x[[2L]] <- call("as", x[[2L]], "matrix")
    x
  }
  l <- vector("list", 2L * length(exprs))
  l[c(TRUE, FALSE)] <- exprs
  l[c(FALSE, TRUE)] <- lapply(exprs, f)
  l <- c(list(quote(microbenchmark::microbenchmark)), l, list(times = times))
  eval(as.call(l), parent.frame())
}
```

## `[` method tests

For `"packedMatrix"`, the subset operator `[` avoids coercion to
`"matrix"` in many special cases. Here is a list of operations to be
tested.

``` r
exprs <- alist(S[],
               S[NULL],
               S[integer(0L)],
               S[c(0:6, NA)],
               S[-(0:6)],
               S[logical(0L)],
               S[TRUE],
               S[FALSE],
               S[NA],
               S[c(TRUE, FALSE, NA)],
               S[, ],
               S[NULL, ],
               S[integer(0L), ],
               S[1L, ],
               S[1L, , drop = FALSE],
               S[c(0:6, NA), ],
               S[-seq_len(nrow(S))[-(1:6)], ],
               S[logical(0L), ],
               S[TRUE, ],
               S[FALSE, ],
               S[NA, ],
               S[c(TRUE, FALSE, NA), ],
               S[character(0L), ],
               S[rownames(S)[1:6], ],
               S[matrix(0L, 0L, 2L)],
               S[cbind(c(0:6, NA), c(NA, 6:0))],
               S[cbind(c(rownames(S), NA), c(NA, colnames(S)))],
               S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))])
```

We start by checking that the `"packedMatrix"` subset is consistent with
the subset via `"matrix"`.

``` r
n <- 10L
S <- new("dspMatrix", uplo = "U", x = as.double(seq_len(choose(n + 1L, 2L))),
         Dim = c(n, n), Dimnames = rep(list(A = sprintf("i%d", seq_len(n))), 2L))
all(checkEquiv(exprs))
```

    ## [1] TRUE

Now we test its performance against the subset via `"matrix"`.

``` r
n <- 1000L
S <- new("dspMatrix", uplo = "U", x = as.double(seq_len(choose(n + 1L, 2L))),
         Dim = c(n, n), Dimnames = rep(list(A = sprintf("i%d", seq_len(n))), 2L))
checkTimes(exprs)
```

    ## Loading required namespace: microbenchmark

    ## Warning in microbenchmark::microbenchmark(S[], as(S, "matrix")[], S[NULL], : less accurate nanosecond times to avoid potential integer overflows

    ## Unit: microseconds
    ##                                                            expr      min        lq        mean    median         uq       max neval
    ##                                                             S[]    1.476    3.8950    10.95561    9.6760    16.3385    33.292   100
    ##                                               as(S, "matrix")[] 1784.607 2143.0290  4176.08206 3336.3340  3967.7545 51560.657   100
    ##                                                         S[NULL]    2.829    7.2775    15.05602    9.7580    23.2470    51.701   100
    ##                                           as(S, "matrix")[NULL] 1371.409 1637.2735  3642.36374 2575.5175  3603.1210 49641.242   100
    ##                                                  S[integer(0L)]    3.444    9.4300    18.87148   20.6025    25.6250    43.583   100
    ##                                    as(S, "matrix")[integer(0L)] 1317.863 1685.3255  3145.08007 2232.3680  3429.3220 51081.162   100
    ##                                                   S[c(0:6, NA)]    4.059   14.8625    24.82796   24.1490    34.7885    76.711   100
    ##                                     as(S, "matrix")[c(0:6, NA)] 1382.930 1648.7945  3085.26886 2554.9765  3314.9115 50178.465   100
    ##                                                       S[-(0:6)] 4360.760 4835.4170  6481.58586 6079.7260  6762.9910 54005.692   100
    ##                                         as(S, "matrix")[-(0:6)] 3481.966 5098.8830  7679.67925 5605.5610  6766.1275 55111.421   100
    ##                                                  S[logical(0L)]    3.116    8.3435    18.69067   20.6640    26.4450    50.348   100
    ##                                    as(S, "matrix")[logical(0L)] 1430.080 1631.9845  3534.11718 1775.5665  3476.4515 53543.622   100
    ##                                                         S[TRUE]  929.880 1087.0330  1566.16433 1147.2005  1397.0340  4753.335   100
    ##                                           as(S, "matrix")[TRUE] 3778.970 4055.1460  5878.53039 5338.9790  6086.3885 52231.540   100
    ##                                                        S[FALSE]    4.510    9.1430    19.27328   19.3110    27.5110    51.373   100
    ##                                          as(S, "matrix")[FALSE] 2463.690 2630.2320  5376.01389 2805.1585  4160.5570 56049.829   100
    ##                                                           S[NA]  886.420  982.2370  1508.58803 1022.2735  2245.5905  4100.984   100
    ##                                             as(S, "matrix")[NA] 3450.642 4025.5235  7319.87596 5407.3670  6477.0775 53630.706   100
    ##                                           S[c(TRUE, FALSE, NA)] 7771.468 8176.0970 10067.52089 9795.2895 10401.4745 57341.575   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)] 3219.566 3456.0745  5191.47084 4886.5850  5476.6980 53665.761   100
    ##                                                           S[, ]    1.517    4.1820    10.98308    9.7785    17.5275    40.180   100
    ##                                             as(S, "matrix")[, ] 2892.632 3462.3475  4992.47406 4635.1115  5110.8960 51348.441   100
    ##                                                       S[NULL, ]    5.002   18.9420    30.82995   32.4310    39.2985   103.935   100
    ##                                         as(S, "matrix")[NULL, ] 1374.812 1647.0520  2428.71085 1859.7190  3300.6435  4954.973   100
    ##                                                S[integer(0L), ]    6.929   16.8100    31.79017   32.4925    44.5875    80.606   100
    ##                                  as(S, "matrix")[integer(0L), ] 1228.237 1661.2790  3424.43849 1767.3665  3068.3375 63492.149   100
    ##                                                         S[1L, ]    7.831   19.3930    33.54046   36.1620    45.8175    66.502   100
    ##                                           as(S, "matrix")[1L, ] 1379.896 1634.9570  3267.67335 1775.6485  3170.2430 50390.599   100
    ##                                           S[1L, , drop = FALSE]    6.888   19.1880    32.85617   33.9480    43.2755    94.956   100
    ##                             as(S, "matrix")[1L, , drop = FALSE] 1438.444 1681.4100  2996.75191 2044.4445  3262.4930 48790.492   100
    ##                                                 S[c(0:6, NA), ]   12.464   26.0965    41.43993   41.1230    54.1815   103.648   100
    ##                                   as(S, "matrix")[c(0:6, NA), ] 1370.384 1671.9390  2959.25044 2015.6010  3262.0215 49955.220   100
    ##                                  S[-seq_len(nrow(S))[-(1:6)], ]   21.402   45.0590    60.90673   58.8760    75.5630   138.867   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:6)], ] 1322.947 1680.2005  2398.16093 1814.2500  3184.3060  5167.476   100
    ##                                                S[logical(0L), ]  144.279  291.0590   309.27817  311.7230   330.9520   640.584   100
    ##                                  as(S, "matrix")[logical(0L), ] 1464.151 1666.1785  3414.20120 1850.3505  3265.9780 50459.438   100
    ##                                                       S[TRUE, ]    7.298   16.2155    27.23999   27.3675    34.8295    66.420   100
    ##                                         as(S, "matrix")[TRUE, ] 3170.325 3671.0580  5158.44821 4766.8650  5203.6585 52747.607   100
    ##                                                      S[FALSE, ]  144.484  279.3945   295.20492  300.3250   322.1370   504.546   100
    ##                                        as(S, "matrix")[FALSE, ] 1444.594 1670.7910  4020.62359 2046.8430  3487.4600 51898.620   100
    ##                                                         S[NA, ] 1073.626 1269.9750  1656.75834 1322.7625  1431.4125  4595.239   100
    ##                                           as(S, "matrix")[NA, ] 2585.255 2879.8605  4036.02360 4047.9505  4519.8605  7744.162   100
    ##                                         S[c(TRUE, FALSE, NA), ]  654.647  789.8650  1574.83460  869.0360   980.0845 51118.021   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ] 2293.253 2680.9285  5222.09989 3949.6325  4509.8975 52631.782   100
    ##                                              S[character(0L), ]    6.355   16.7895    34.23049   31.4265    46.3505   133.865   100
    ##                                as(S, "matrix")[character(0L), ] 1397.075 1678.4170  3542.04043 1912.4245  3244.7195 57480.524   100
    ##                                           S[rownames(S)[1:6], ]   29.274   46.5555    67.33225   68.7570    81.1800   136.858   100
    ##                             as(S, "matrix")[rownames(S)[1:6], ] 1379.158 1690.1020  4000.68447 1934.5440  3194.6995 51914.282   100
    ##                                           S[matrix(0L, 0L, 2L)]   18.122   45.8790    55.47628   55.4525    65.1490   154.365   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)] 1383.340 1617.8190  2377.58016 1788.2765  3100.6865  6192.804   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]   25.051   60.2905    74.33259   75.6450    90.5690   179.990   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))] 1444.348 1639.1595  3569.60514 2296.5945  3350.6430 51971.518   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  128.576  197.4560   213.28897  212.9335   227.7550   377.282   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))] 1501.338 1733.0495  2869.78967 1937.1885  3128.0950 49401.433   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 7653.183 8929.0620 11525.65801 9664.8890 10265.4775 58984.035   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 3958.345 5366.7975  8686.72822 5853.9595  6785.4180 54478.914   100

The `"packedMatrix"` subset is slower in some cases because:

-   It takes extra steps to avoid a coercion to `"matrix"`.
-   Some operations are implemented in R rather than C (for now).

However, for large enough `n`, there is no comparison…

``` r
n <- 30000L
i <- replace(logical(1000L), 1L, TRUE)
S <- new("dspMatrix", uplo = "U", x = as.double(seq_len(choose(n + 1L, 2L))), Dim = c(n, n))
{S[i]; cat("done!\n")}
```

    ## done!

``` r
{as(S, "matrix")[i]; cat("done!\n")}
```

    ## Error: vector memory exhausted (limit reached?)

## `t` method tests

``` r
selectMethod("t", signature(x = "packedMatrix"))
```

    ## Method Definition:
    ## 
    ## function (x) 
    ## .Call(packedMatrix_t, x)
    ## <bytecode: 0x139075d98>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
n <- 1000L
L <- new("dtpMatrix", uplo = "L", x = as.double(seq_len(choose(n + 1L, 2L))), 
         Dim = c(n, n), Dimnames = list(A = NULL, B = sprintf("j%d", seq_len(n))))
identical(as(t(L), "matrix"), t(as(L, "matrix")))
```

    ## [1] TRUE

``` r
microbenchmark::microbenchmark(t(L), as(t(as(L, "dtrMatrix")), "dtpMatrix"))
```

    ## Unit: microseconds
    ##                                    expr      min        lq      mean    median       uq      max neval
    ##                                    t(L)  615.246  644.4585  754.8707  755.4865  830.947 1222.620   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4407.090 5175.1840 5376.5166 5322.1075 5507.305 6787.427   100

``` r
n <- 30000L
L <- new("dtpMatrix", uplo = "L", x = as.double(seq_len(choose(n + 1L, 2L))), Dim = c(n, n))
{t(L); cat("done!\n")}
```

    ## done!

``` r
{as(t(as(L, "dtrMatrix")), "dtpMatrix"); cat("done!\n")}
```

    ## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': vector memory exhausted (limit reached?)

## `diag` method tests

``` r
selectMethod("diag", signature(x = "packedMatrix"))
```

    ## Method Definition:
    ## 
    ## function (x = 1, nrow, ncol, names = TRUE) 
    ## .Call(packedMatrix_diag_get, x, names)
    ## <bytecode: 0x1123369b8>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
n <- 4L
S <- new("lspMatrix", uplo = "U", x = rep_len(c(TRUE, FALSE, NA), choose(n + 1L, 2L)),
         Dim = c(n, n), Dimnames = list(letters[seq_len(n)], NULL))
S
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE    NA FALSE FALSE
    ## c  TRUE FALSE    NA    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(S)
```

    ##    a    b    c    d 
    ## TRUE   NA   NA TRUE

``` r
diag(S, names = FALSE)
```

    ## [1] TRUE   NA   NA TRUE

## `diag<-` method tests

``` r
selectMethod("diag<-", signature(x = "packedMatrix"))
```

    ## Method Definition:
    ## 
    ## function (x, value) 
    ## .Call(packedMatrix_diag_set, x, value)
    ## <bytecode: 0x1123e57b8>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
S
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE    NA FALSE FALSE
    ## c  TRUE FALSE    NA    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(S) <- logical(n)
S
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a FALSE FALSE  TRUE  TRUE
    ## b FALSE FALSE FALSE FALSE
    ## c  TRUE FALSE FALSE    NA
    ## d  TRUE FALSE    NA FALSE

``` r
diag(S) <- TRUE
S
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE  TRUE FALSE FALSE
    ## c  TRUE FALSE  TRUE    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(S) <- 1 # 'x' coerced to "dspMatrix"
S
```

    ## 4 x 4 Matrix of class "dspMatrix"
    ##   a b  c  d
    ## a 1 0  1  1
    ## b 0 1  0  0
    ## c 1 0  1 NA
    ## d 1 0 NA  1

``` r
diag(S) <- TRUE # 'value' coerced to double
S
```

    ## 4 x 4 Matrix of class "dspMatrix"
    ##   a b  c  d
    ## a 1 0  1  1
    ## b 0 1  0  0
    ## c 1 0  1 NA
    ## d 1 0 NA  1

``` r
diag(S) <- double(2L) # right type, wrong length
```

    ## Error in `diag<-`(`*tmp*`, value = c(0, 0)): replacement diagonal has wrong length

``` r
diag(S) <- vector("list", n) # right length, wrong type
```

    ## Error in `diag<-`(`*tmp*`, value = list(NULL, NULL, NULL, NULL)): replacement diagonal has incompatible type 'list'
