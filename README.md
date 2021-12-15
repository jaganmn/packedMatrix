packedMatrix
================
2021-12-15

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
    ##                                                  S[integer(0L)]     3.526     8.1180    19.29788    18.8395    25.9940     66.912   100
    ##                                    as(S, "matrix")[integer(0L)]  1410.072  1652.1770  3002.90191  1893.7285  2974.9190  57876.174   100
    ##                                                   S[c(0:6, NA)]    16.482    39.1345    60.32043    58.6710    71.7500    183.065   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1409.580  1662.8165  3672.95302  1918.3490  3363.6400  61961.045   100
    ##                                                       S[-(0:6)] 47213.263 51260.5985 68110.74652 54567.8020 98315.0685 147071.961   100
    ##                                         as(S, "matrix")[-(0:6)]  3514.766  3946.3320  7171.20094  4962.8450  6086.4705  62930.367   100
    ##                                                  S[logical(0L)]     4.100     8.9585    19.38357    18.1015    25.0100     68.593   100
    ##                                    as(S, "matrix")[logical(0L)]  1161.530  1653.1405  2461.32266  1827.9645  3009.7280   7273.564   100
    ##                                                         S[TRUE]   791.464  1117.1475  1548.32564  1213.7230  1486.7215   4876.089   100
    ##                                           as(S, "matrix")[TRUE]  3556.340  4086.2035  5972.47246  5038.7565  6259.0395  56531.251   100
    ##                                                        S[FALSE]     4.551     8.9790    19.03999    18.2245    27.4905     53.915   100
    ##                                          as(S, "matrix")[FALSE]  2468.036  2661.9660  3926.60936  2828.9795  4050.4105  56946.130   100
    ##                                                           S[NA]   654.442   975.2875  2029.58651  1052.2650  1474.7905  52737.439   100
    ##                                             as(S, "matrix")[NA]  3393.611  3887.7840  5305.01706  5295.1705  6198.7285  12540.301   100
    ##                                           S[c(TRUE, FALSE, NA)] 34276.410 36861.8905 43016.92366 38396.4590 40522.9240  97695.251   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  3255.646  3503.3680  4817.07319  4395.2205  5677.1060  10617.319   100
    ##                                                S[integer(0L), ]   137.514   213.5690   284.10663   292.8425   334.0475    642.429   100
    ##                                  as(S, "matrix")[integer(0L), ]  1245.457  1658.6345  2636.83751  1783.9510  2202.3970  51515.393   100
    ##                                                         S[1L, ]    29.643    59.0810    79.73639    81.2005    95.8170    206.230   100
    ##                                           as(S, "matrix")[1L, ]  1386.169  1703.4270  4137.21816  2604.1355  4076.7120  52137.445   100
    ##                                           S[1L, , drop = FALSE]   172.446   306.6185   353.06658   350.9600   384.6210    954.890   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1110.157  1654.1245  2295.35466  1809.0635  2399.7915  10113.880   100
    ##                                                 S[c(0:1, NA), ]   173.389   290.4235   332.15207   337.7375   389.4385    620.207   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]  1213.641  1660.4795  3200.18653  1769.6010  3277.3965  62483.918   100
    ##                                                 S[c(0:2, NA), ]  1075.184  1859.3910  3067.79507  2043.8910  2959.7900  53691.837   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1333.443  1643.7925  2315.60169  1759.5150  2226.7100   8634.559   100
    ##                                          S[-seq_len(nrow(S)), ]   149.896   231.0760   288.78596   302.6210   336.6305    619.428   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1277.683  1668.2900  3519.31003  1839.0550  3120.6125  52463.600   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    42.025    65.9690    85.06270    83.9680   104.5705    193.110   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1235.658  1712.0780  4045.64589  1888.9110  3228.1760  55467.629   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1616.958  1805.2915  2773.90133  1943.2155  3413.7625  10362.750   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1454.352  1671.3855  2462.16972  1834.9755  3007.6165   6897.717   100
    ##                                                S[logical(0L), ]   141.286   227.9395   283.15543   296.9015   335.7900    550.179   100
    ##                                  as(S, "matrix")[logical(0L), ]  1371.737  1666.2195  3484.64084  1903.7325  3300.2950  49231.652   100
    ##                                                       S[TRUE, ]     5.043    13.0175    24.95301    22.8165    35.2395     59.327   100
    ##                                         as(S, "matrix")[TRUE, ]  2746.631  3316.6130  5874.11223  3649.8200  5212.2480  56871.756   100
    ##                                                      S[FALSE, ]   143.582   245.9385   291.23407   282.3055   327.1390    696.262   100
    ##                                        as(S, "matrix")[FALSE, ]  1004.705  1647.2570  2619.56052  1791.1465  3222.2310   8353.545   100
    ##                                                         S[NA, ]   990.478  1247.1585  1730.11595  1328.3385  1706.4815   5145.131   100
    ##                                           as(S, "matrix")[NA, ]  2460.738  2869.3030  4383.19807  3375.2020  4728.7965  53369.495   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2610.429  3285.1660  4657.55736  3960.7435  5829.5235  11154.009   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  2493.046  2705.0160  4483.48489  3118.9520  4845.5235  56061.678   100
    ##                                              S[character(0L), ]   149.322   230.4200   295.48249   301.2885   334.3345    647.800   100
    ##                                as(S, "matrix")[character(0L), ]  1178.463  1657.8965  2619.45433  1829.0305  3404.5990   7042.980   100
    ##                                            S[rownames(S)[1L], ]    32.267    59.9625   133.07698    80.6880   100.1835   5118.809   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1169.238  1698.3840  3223.17810  1912.6910  3387.9325  63109.209   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   185.033   304.2610   338.98390   342.7600   371.8085    715.614   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1319.298  1665.8300  3121.12910  1877.4720  3276.6380  55832.119   100
    ##                                           S[rownames(S)[1:2], ]  1478.460  1886.1845  2609.41630  2015.1705  3094.7415   9496.953   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1245.826  1681.4510  2960.03149  1873.8435  3036.8905  52215.181   100
    ##                                           S[matrix(0L, 0L, 2L)]    17.343    38.7655    51.75307    51.9265    62.4020     86.674   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1296.338  1657.1175  3161.62193  1813.7170  3465.9965  54507.122   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    27.306    54.9400    74.92053    70.0280    91.9220    200.613   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1257.839  1672.3285  2911.86059  1834.8115  2970.5320  52819.644   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   130.872   192.0440   216.99824   212.6055   234.2125    387.573   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1307.162  1754.3285  2412.77784  1870.4200  3014.3610   7603.245   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 28634.646 31515.9005 41562.47080 33054.4460 35805.3205  93733.216   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3178.853  4214.1645  6780.31719  5571.0185  7095.5215  55507.481   100
    ##                                                         S[NULL]     2.542     7.7695    16.94448    14.8830    22.2015     71.996   100
    ##                                           as(S, "matrix")[NULL]  1125.081  1649.9015  2797.23484  1785.4680  2998.6580  51380.175   100
    ##                                                       S[NULL, ]   140.548   252.9905   291.96059   296.3070   332.5305    649.071   100
    ##                                         as(S, "matrix")[NULL, ]  1271.123  1664.2720  2408.47899  1875.6270  3201.0545   7891.598   100
    ##                                                       S[, NULL]   141.286   229.5180   283.75977   288.4555   324.9865    695.073   100
    ##                                         as(S, "matrix")[, NULL]  1165.015  1648.9585  3046.94944  1816.0745  3214.0515  58389.207   100

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
    ## <bytecode: 0x13140ac18>
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
    ##                                    expr      min        lq      mean   median       uq      max neval
    ##                                    t(L)  541.405  629.0835  746.0786  731.604  828.118 1193.838   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4141.492 4976.4160 5152.0141 5078.957 5203.392 7421.697   100

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
    ## <bytecode: 0x12a6b91a0>
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
    ## <bytecode: 0x12a767fa0>
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
diag(S) <- 1 # 'x' coerced to "d?pMatrix"
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
