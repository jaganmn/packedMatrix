packedMatrix
================
2021-12-15

Based on discussion on R-devel beginning
[here](https://stat.ethz.ch/pipermail/r-devel/2021-November/081261.html)
and continued privately with Martin Mächler.

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
exprs <- alist(S[integer(0L)],
               S[c(0:6, NA)],
               S[-(0:6)],
               S[logical(0L)],
               S[TRUE],
               S[FALSE],
               S[NA],
               S[c(TRUE, FALSE, NA)],
               S[integer(0L), ],
               S[1L, ],
               S[1L, , drop = FALSE],
               S[c(0:1, NA), ],
               S[c(0:2, NA), ],
               S[-seq_len(nrow(S)), ],
               S[-seq_len(nrow(S))[-1L], ],
               S[-seq_len(nrow(S))[-(1:2)], ],
               S[logical(0L), ],
               S[TRUE, ],
               S[FALSE, ],
               S[NA, ],
               S[c(TRUE, FALSE, NA), ],
               S[character(0L), ],
               S[rownames(S)[1L], ],
               S[rownames(S)[1L], , drop = FALSE],
               S[rownames(S)[1:2], ],
               S[matrix(0L, 0L, 2L)],
               S[cbind(c(0:6, NA), c(NA, 6:0))],
               S[cbind(c(rownames(S), NA), c(NA, colnames(S)))],
               S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))],
               S[NULL],
               S[NULL, ],
               S[, NULL])
```

We start by checking that the `"packedMatrix"` subset is consistent with
the subset via `"matrix"`.

``` r
n <- 10L
S <- new("dspMatrix", uplo = "U", x = as.double(seq_len(choose(n + 1L, 2L))),
         Dim = c(n, n), Dimnames = list(A = sprintf("i%d", seq_len(n)), B = NULL))
all(checkEquiv(exprs))
```

    ## [1] TRUE

Now we test its performance against the subset via `"matrix"`.

``` r
n <- 1000L
S <- new("dspMatrix", uplo = "U", x = as.double(seq_len(choose(n + 1L, 2L))),
         Dim = c(n, n), Dimnames = list(A = sprintf("i%d", seq_len(n)), B = NULL))
checkTimes(exprs)
```

    ## Loading required namespace: microbenchmark

    ## Warning in microbenchmark::microbenchmark(S[integer(0L)], as(S, "matrix")[integer(0L)], : less accurate nanosecond times to avoid potential integer overflows

    ## Unit: microseconds
    ##                                                            expr       min         lq        mean     median         uq        max neval
    ##                                                  S[integer(0L)]     2.829     7.5440    16.55867    16.9740    24.5590     37.228   100
    ##                                    as(S, "matrix")[integer(0L)]  1114.216  1595.0435  2360.84970  1699.6550  2726.6640   8375.972   100
    ##                                                   S[c(0:6, NA)]    17.343    36.4080    53.20406    53.1770    65.7845    103.443   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1347.834  1611.8125  2703.37969  1697.4410  2367.4015  52401.895   100
    ##                                                       S[-(0:6)] 47505.757 49914.3430 61245.36745 51870.5760 55647.8855 120729.625   100
    ##                                         as(S, "matrix")[-(0:6)]  3193.408  3779.1955  6275.89132  5194.1260  6604.9360  59263.901   100
    ##                                                  S[logical(0L)]     2.952     8.2820    16.38565    16.2975    21.9760     45.715   100
    ##                                    as(S, "matrix")[logical(0L)]  1454.885  1628.2945  2755.24100  1799.8590  3505.6845  10504.774   100
    ##                                                         S[TRUE]   704.093  1088.6935  1514.53918  1150.5830  1228.8520   7681.842   100
    ##                                           as(S, "matrix")[TRUE]  3435.677  4008.7750  6262.10179  5209.3985  6190.0775  54349.272   100
    ##                                                        S[FALSE]     2.665    10.2295    18.54676    18.3065    25.0305     46.453   100
    ##                                          as(S, "matrix")[FALSE]  2190.138  2590.2980  3883.32156  2694.8275  3208.9880  54253.578   100
    ##                                                           S[NA]   820.943   964.9555  1338.79883   990.9495  1111.4075   6299.158   100
    ##                                             as(S, "matrix")[NA]  2784.187  3760.6430  4932.75264  4032.2885  5765.1945  11685.000   100
    ##                                           S[c(TRUE, FALSE, NA)] 32872.816 35881.3140 46317.40685 38201.0735 41629.4115 100798.172   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  2791.895  3427.0055  7062.79981  3687.6835  5242.4240  63528.803   100
    ##                                                S[integer(0L), ]   140.343   226.0125   284.29195   283.9865   331.3005    679.493   100
    ##                                  as(S, "matrix")[integer(0L), ]  1392.032  1639.3030  4002.41795  1763.2870  3310.6475  56055.528   100
    ##                                                         S[1L, ]    27.265    51.4345   119.11156    71.2785    82.0615   5205.278   100
    ##                                           as(S, "matrix")[1L, ]  1168.377  1626.7365  2360.23593  1728.1090  3066.0825   7144.824   100
    ##                                           S[1L, , drop = FALSE]   162.606   287.5945   316.71311   330.4395   359.9390    753.334   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1376.657  1647.9335  2370.98777  1714.8250  2761.7805   8879.862   100
    ##                                                 S[c(0:1, NA), ]   168.305   252.9290   314.88205   331.9360   364.4900    559.404   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]  1016.390  1627.8230  2705.74867  1709.3925  2137.8630  52995.083   100
    ##                                                 S[c(0:2, NA), ]  1507.037  1801.5400  2545.52231  1948.6070  2773.7525   7338.549   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1431.351  1623.8870  2486.58809  1742.2540  3277.9500   7789.180   100
    ##                                          S[-seq_len(nrow(S)), ]   146.657   242.1050   280.16366   290.3620   319.5745    731.768   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1142.342  1658.9625  2356.25155  1748.2400  2996.1775   8547.926   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    38.581    63.0785    82.62361    80.5650   100.1015    189.379   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1182.153  1674.6245  2896.94479  1777.3500  3157.3485  54513.067   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1311.016  1862.6710  3434.43101  1968.5740  3585.9010  54799.165   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1245.990  1625.6090  2478.85713  1774.5210  3215.0355   8108.201   100
    ##                                                S[logical(0L), ]   141.081   218.9605   266.61767   282.5515   319.3900    457.642   100
    ##                                  as(S, "matrix")[logical(0L), ]  1214.625  1632.3535  2318.44299  1708.1010  2726.1310   7534.570   100
    ##                                                       S[TRUE, ]     6.724    14.6165    25.17359    23.8415    31.7750    109.183   100
    ##                                         as(S, "matrix")[TRUE, ]  2611.085  3276.3920  4853.76942  3613.1865  5061.7985  61020.710   100
    ##                                                      S[FALSE, ]   141.532   220.1495   307.78003   289.3370   323.3055   4647.391   100
    ##                                        as(S, "matrix")[FALSE, ]  1464.274  1612.7760  2240.55898  1712.6315  2319.2470   6759.588   100
    ##                                                         S[NA, ]   780.148  1214.7480  1638.06316  1300.3150  1415.8120   6770.084   100
    ##                                           as(S, "matrix")[NA, ]  2234.008  2760.5915  5640.71522  3141.7685  4647.4115  65879.333   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2596.038  3167.3935  5572.35305  3358.0640  4985.0260  55714.572   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  2312.564  2639.6825  4144.34273  2945.8705  4290.7730  55293.953   100
    ##                                              S[character(0L), ]   141.696   220.0470   275.70819   286.2415   313.1990    623.036   100
    ##                                as(S, "matrix")[character(0L), ]  1249.885  1636.0845  2804.11054  1747.0920  2928.9990  54148.741   100
    ##                                            S[rownames(S)[1L], ]    36.736    60.8235    81.10292    81.2005    97.7850    226.115   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1474.032  1645.8015  2485.29372  1762.8155  3054.9305   8629.598   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   173.512   276.9345   330.27714   349.5250   369.7380    622.790   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1181.415  1635.8590  2261.48579  1704.9030  2149.1585   7529.855   100
    ##                                           S[rownames(S)[1:2], ]  1066.082  1824.7050  2782.15545  1968.1435  3453.5120   8791.466   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1213.067  1642.6445  2437.95102  1746.1080  3085.9675   8140.960   100
    ##                                           S[matrix(0L, 0L, 2L)]    17.343    39.4010    52.49271    50.9220    61.4795    154.324   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1403.758  1619.6640  3513.13174  1680.1800  2833.1205  59882.837   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    25.912    55.7395    73.46339    76.0140    88.2320    107.953   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1156.077  1628.4175  2253.10170  1762.5695  2800.4640   5716.794   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   133.906   186.7140   212.13810   209.6535   227.8985    370.558   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1396.214  1708.3880  2922.81579  1848.7720  3036.6240  49076.303   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 28076.513 30782.0005 37923.42314 31994.7600 33649.3560  84775.167   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3538.997  4068.6965  5313.49340  4653.4590  5811.1555  14663.117   100
    ##                                                         S[NULL]     2.870     7.7080    15.84363    15.5390    21.7300     45.633   100
    ##                                           as(S, "matrix")[NULL]  1239.881  1620.6890  3805.48470  1749.6545  3115.0570  54400.153   100
    ##                                                       S[NULL, ]   146.411   251.8220   278.41624   286.2210   313.8755    709.833   100
    ##                                         as(S, "matrix")[NULL, ]  1389.654  1613.0630  2203.03332  1733.8285  2715.6555   5962.835   100
    ##                                                       S[, NULL]   138.826   173.0610   258.05728   273.1420   311.6410    483.431   100
    ##                                         as(S, "matrix")[, NULL]  1368.785  1588.7705  2907.08245  1714.7020  3017.4975  52898.610   100

The `"packedMatrix"` subset is slower in some cases because:

-   It is implemented in R rather than C (for now).
-   It takes extra steps to avoid a coercion to `"matrix"`.

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
    ## <bytecode: 0x1287f5450>
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
    ##                                    expr      min       lq      mean   median        uq      max neval
    ##                                    t(L)  540.626  643.249  748.0696  765.757  838.7985 1039.760   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4313.733 5188.468 5349.4922 5349.824 5502.0155 7211.859   100

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
    ## <bytecode: 0x10e5d4ef8>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
n <- 4L
U <- new("lspMatrix", uplo = "U", x = rep_len(c(TRUE, FALSE, NA), choose(n + 1L, 2L)),
         Dim = c(n, n), Dimnames = list(letters[seq_len(n)], NULL))
U
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE    NA FALSE FALSE
    ## c  TRUE FALSE    NA    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(U)
```

    ##    a    b    c    d 
    ## TRUE   NA   NA TRUE

``` r
diag(U, names = FALSE)
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
    ## <bytecode: 0x10e683cf8>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
U
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE    NA FALSE FALSE
    ## c  TRUE FALSE    NA    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(U) <- logical(n)
U
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a FALSE FALSE  TRUE  TRUE
    ## b FALSE FALSE FALSE FALSE
    ## c  TRUE FALSE FALSE    NA
    ## d  TRUE FALSE    NA FALSE

``` r
diag(U) <- TRUE
U
```

    ## 4 x 4 Matrix of class "lspMatrix"
    ##       a     b     c     d
    ## a  TRUE FALSE  TRUE  TRUE
    ## b FALSE  TRUE FALSE FALSE
    ## c  TRUE FALSE  TRUE    NA
    ## d  TRUE FALSE    NA  TRUE

``` r
diag(U) <- 1 # 'x@x' coerced to 'typeof(value)'
U
```

    ## 4 x 4 Matrix of class "dspMatrix"
    ##   a b  c  d
    ## a 1 0  1  1
    ## b 0 1  0  0
    ## c 1 0  1 NA
    ## d 1 0 NA  1

``` r
diag(U) <- TRUE # 'value' coerced to 'typeof(x@x)'
U
```

    ## 4 x 4 Matrix of class "dspMatrix"
    ##   a b  c  d
    ## a 1 0  1  1
    ## b 0 1  0  0
    ## c 1 0  1 NA
    ## d 1 0 NA  1

``` r
diag(U) <- double(2L) # right type, wrong length
```

    ## Error in `diag<-`(`*tmp*`, value = c(0, 0)): replacement diagonal has wrong length

``` r
diag(U) <- vector("list", n) # right length, wrong type
```

    ## Error in `diag<-`(`*tmp*`, value = list(NULL, NULL, NULL, NULL)): replacement diagonal has incompatible type 'list'
