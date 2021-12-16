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
    ##                                                  S[integer(0L)]     3.034     8.3435    17.08142    15.7645    24.3950     48.954   100
    ##                                    as(S, "matrix")[integer(0L)]  1117.291  1616.5275  2647.23798  1677.1050  1985.4660  50059.237   100
    ##                                                   S[c(0:6, NA)]    16.318    43.9725    59.37005    58.8350    70.6430    174.947   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1166.286  1642.9315  2887.82311  1742.0900  2973.8735  54056.368   100
    ##                                                       S[-(0:6)] 47611.414 51033.6430 62014.65578 52720.1165 55852.7830 117184.847   100
    ##                                         as(S, "matrix")[-(0:6)]  3193.326  3826.3045  6228.40635  4704.4835  6238.5395  55049.347   100
    ##                                                  S[logical(0L)]     3.444     9.0200    19.16217    17.3430    24.7435     80.073   100
    ##                                    as(S, "matrix")[logical(0L)]  1350.827  1621.4065  2860.31375  1702.0945  2073.9235  52198.904   100
    ##                                                         S[TRUE]   698.640  1089.6980  1513.36084  1147.9795  1326.1655   6479.230   100
    ##                                           as(S, "matrix")[TRUE]  3637.561  3947.9515  5297.10980  4149.2615  5614.2325  55636.385   100
    ##                                                        S[FALSE]     3.567     8.4665    18.42663    19.2905    25.1125     55.309   100
    ##                                          as(S, "matrix")[FALSE]  2319.124  2570.3515  3234.50435  2678.7555  2998.9450  10079.727   100
    ##                                                           S[NA]   647.677   976.3125  1323.98061  1004.7460  1076.1065   7604.639   100
    ##                                             as(S, "matrix")[NA]  3430.019  3810.4580  6475.36657  4285.3405  6173.2470  56168.852   100
    ##                                           S[c(TRUE, FALSE, NA)] 33956.774 36561.3195 45172.92342 37922.7245 40443.8350  92331.918   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  2602.926  3396.3170  5002.94177  3659.2705  5328.4420  54234.431   100
    ##                                                S[integer(0L), ]   137.514   168.2025   271.46182   283.0435   320.3535    886.133   100
    ##                                  as(S, "matrix")[integer(0L), ]  1474.975  1622.0215  2537.93239  1732.0245  3414.8695   6779.268   100
    ##                                                         S[1L, ]    28.577    58.2200    74.23378    74.8250    84.7880    180.564   100
    ##                                           as(S, "matrix")[1L, ]  1400.601  1631.5745  2261.63995  1729.8105  2073.4110   8749.892   100
    ##                                           S[1L, , drop = FALSE]   174.619   277.9390   318.75696   328.6355   360.1030    569.162   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1139.841  1640.6150  2904.48551  1759.4535  3028.6085  55586.037   100
    ##                                                 S[c(0:1, NA), ]   171.257   305.0400   327.17713   333.9450   372.1365    601.142   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]  1312.984  1659.2085  3408.33369  1721.9385  2995.6650  54448.861   100
    ##                                                 S[c(0:2, NA), ]  1074.364  1822.2860  3367.27793  1985.1380  3926.0780  51329.376   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1059.727  1649.2045  2086.63801  1718.6175  1967.2415   5044.189   100
    ##                                          S[-seq_len(nrow(S)), ]   150.716   172.4255   277.03864   304.6505   333.9655    398.028   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1436.722  1641.4965  2383.62971  1736.6370  3039.5145   8762.807   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    37.228    70.2125    86.37142    86.8175   103.3405    191.839   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1215.445  1637.2940  2129.32844  1702.8325  1910.5180   8113.285   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1623.108  1823.0240  3232.79793  1912.4245  3551.9325  54020.001   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1437.419  1640.1230  2401.63240  1763.3690  2799.8900   8879.616   100
    ##                                                S[logical(0L), ]   140.384   214.1635   312.97473   291.8175   317.7090   4198.646   100
    ##                                  as(S, "matrix")[logical(0L), ]  1020.736  1640.3280  2346.57268  1733.7875  2897.7980   7004.604   100
    ##                                                       S[TRUE, ]     7.216    13.8375    28.07721    27.2240    38.3965     90.405   100
    ##                                         as(S, "matrix")[TRUE, ]  2369.841  3264.1535  5495.16604  3522.1665  5815.8910  57160.273   100
    ##                                                      S[FALSE, ]   140.548   217.3205   268.45734   286.9795   313.5885    570.515   100
    ##                                        as(S, "matrix")[FALSE, ]  1120.407  1624.2150  2521.06991  1709.2900  2970.4910   9460.176   100
    ##                                                         S[NA, ]  1061.326  1228.0730  2153.90630  1279.5895  1369.7280  54956.277   100
    ##                                           as(S, "matrix")[NA, ]  2536.096  2776.0690  4933.60626  2953.8450  4564.4275  58038.903   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2536.383  3180.0420  4243.75420  3416.2225  5119.6290   9872.513   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  2254.508  2626.8495  3768.32107  2938.2855  4710.8590   9292.199   100
    ##                                              S[character(0L), ]   147.641   253.2365   278.64256   292.2685   324.1255    399.463   100
    ##                                as(S, "matrix")[character(0L), ]  1224.014  1609.2295  2944.99392  1707.5475  2845.4000  54724.299   100
    ##                                            S[rownames(S)[1L], ]    32.308    58.1995    78.26900    78.9455    94.8945    204.877   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1356.444  1661.8940  2922.42383  1756.2555  3062.4130  54740.535   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   174.619   267.5455   320.32972   341.0380   379.2910    438.700   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1209.951  1633.6040  2446.37570  1738.6050  2538.1460  10148.730   100
    ##                                           S[rownames(S)[1:2], ]  1447.013  1850.0225  2694.86768  1993.7070  3454.9675   8008.366   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1452.220  1655.6825  2889.50493  1759.0230  2782.3215  50141.114   100
    ##                                           S[matrix(0L, 0L, 2L)]    20.172    41.1845    53.57060    54.0380    62.1355    104.345   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1419.871  1646.1705  2931.24047  1757.8750  3043.5325  51248.073   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    24.518    54.9400    75.24689    75.7475    90.1385    178.350   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1166.450  1607.2205  2362.19409  1712.5700  2845.2155   9827.003   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   136.079   186.9395   216.01875   214.4710   234.4995    433.001   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1462.470  1709.5155  2852.13466  1799.2030  2026.7940  55249.304   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 28326.080 31085.0110 35030.94284 32265.9750 33731.9095  84624.984   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3595.495  4114.1450  6794.70122  4943.6365  6560.4100  58461.900   100
    ##                                                         S[NULL]     3.895     9.1840    17.25813    15.9900    23.3905     44.772   100
    ##                                           as(S, "matrix")[NULL]  1216.839  1607.4255  3209.01547  1713.9640  2453.8295  52361.961   100
    ##                                                       S[NULL, ]   141.327   219.9445   269.57418   278.0825   306.5160    582.487   100
    ##                                         as(S, "matrix")[NULL, ]  1344.308  1613.0425  2861.38959  1711.0325  3076.8040  55114.660   100
    ##                                                       S[, NULL]   140.876   236.7545   269.25684   284.4990   316.6635    483.964   100
    ##                                         as(S, "matrix")[, NULL]  1278.913  1627.5770  2447.32608  1743.2585  2941.7295   8649.278   100

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
    ## <bytecode: 0x124bef450>
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
    ##                                    expr      min       lq      mean   median       uq      max neval
    ##                                    t(L)  539.068  637.837  735.8873  720.124  817.089 1240.865   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4379.046 5015.017 5185.0297 5115.652 5320.550 6504.240   100

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
    ## <bytecode: 0x123a717a0>
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
    ## <bytecode: 0x123b205a0>
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
