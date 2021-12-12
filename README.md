packedMatrix
================
2021-12-12

Based on discussion on R-devel beginning
[here](https://stat.ethz.ch/pipermail/r-devel/2021-November/081261.html)
and continued privately with Martin Mächler.

## Utilities

A basic helper for constructing `"packedMatrix"`, a virtual superclass
of `"[dln][st]pMatrix"`:

``` r
## Construct an 'n'-by-'n' "dspMatrix" storing a sequence of integers
newDsp <- function(n) {
  x <- as.double(seq_len(choose(n + 1L, 2L)))
  new("dspMatrix", uplo = "U", x = x, Dim = c(n, n),
      Dimnames = list(A = sprintf("i%d", seq_len(n)), B = NULL))
}
```

And two nonstandard evaluation tricks to avoid code duplication in
tests:

``` r
## For each call FUN(<packedMatrix>, ...) passed as argument,
## test equivalence to FUN(as(<packedMatrix>, "matrix"), ...)
checkEquiv <- function(...) {
  e <- parent.frame()
  equiv1 <- function(x1) {
    stopifnot(is.call(x1) && length(x1) > 1L)
    x2 <- x1
    x2[[2L]] <- call("as", x2[[2L]], "matrix") 
    m1 <- eval(x1, e)
    m2 <- eval(x2, e)
    if (is(m1, "Matrix")) {
      m1 <- as(m1, "matrix")
    }
    identical(m1, m2)
  }
  l <- as.list(match.call())[-1L]
  res <- vapply(l, equiv1, NA)
  names(res) <- vapply(l, deparse1, "")
  res
}

## Benchmark the calls FUN(<packedMatrix>, ...) passed as arguments,
## together with the calls FUN(as(<packedMatrix>, "matrix"), ...)
checkTimes <- function(..., times = 100L) {
  e <- parent.frame()
  f <- function(x) {
    stopifnot(is.call(x) && length(x) > 1L)
    x[[2L]] <- call("as", x[[2L]], "matrix")
    x
  }
  l1 <- as.list(match.call())[-1L]
  l1[["times"]] <- NULL
  l2 <- vector("list", 2L * length(l1))
  l2[c(TRUE, FALSE)] <- l1
  l2[c(FALSE, TRUE)] <- lapply(l1, f)
  cl <- as.call(c(list(quote(microbenchmark::microbenchmark)), l2, list(times = times)))
  eval(cl, e)
}
```

## Subset tests

For `"packedMatrix"`, the subset operator `[` avoids coercion to
`"matrix"` in many special cases. Here is a list of subset operations to
be tested.

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
library("Matrix")
S <- newDsp(10L)
all(do.call(checkEquiv, exprs))
```

    ## [1] TRUE

Now we check that, in many cases, the `"packedMatrix"` subset is faster.

``` r
S <- newDsp(1000L)
do.call(checkTimes, exprs)
```

    ## Warning in microbenchmark::microbenchmark(S[integer(0L)], as(S, "matrix")[integer(0L)], : less accurate nanosecond times to avoid potential integer overflows

    ## Unit: microseconds
    ##                                                            expr       min         lq        mean     median         uq        max neval
    ##                                                  S[integer(0L)]     2.583     7.8720    18.24910    18.5115    24.3950     49.815   100
    ##                                    as(S, "matrix")[integer(0L)]  1361.487  1715.1940  4168.09813  1851.9085  3193.4490  56884.138   100
    ##                                                   S[c(0:6, NA)]    16.359    37.3920    59.17981    57.9330    73.7180    150.757   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1538.935  1677.8430  2296.29028  1773.4550  2427.4255   8852.556   100
    ##                                                       S[-(0:6)] 48555.767 51241.4720 63699.05632 53559.7965 57627.2425 149982.264   100
    ##                                         as(S, "matrix")[-(0:6)]  3372.496  4025.7490  7111.66525  5646.5815  6506.0235  61167.736   100
    ##                                                  S[logical(0L)]     2.788     7.7695    17.36186    15.7440    23.3085     92.004   100
    ##                                    as(S, "matrix")[logical(0L)]  1515.442  1690.4300  2574.82911  1820.4205  3189.6565  11346.791   100
    ##                                                         S[TRUE]  1031.888  1164.8715  1485.48084  1224.2805  1336.8050   7308.127   100
    ##                                           as(S, "matrix")[TRUE]  3833.254  4169.3310  6908.09902  4798.1070  6519.6355  60050.240   100
    ##                                                        S[FALSE]     2.706     8.5280    19.60989    18.7165    26.6295     60.434   100
    ##                                          as(S, "matrix")[FALSE]  2539.991  2713.4005  3448.42390  2835.5395  3691.6810   8577.569   100
    ##                                                           S[NA]   628.202   995.5005  1764.93274  1024.7130  1117.2295  56477.664   100
    ##                                             as(S, "matrix")[NA]  3468.682  3885.3650  5413.65230  4374.1055  5883.6025  54463.047   100
    ##                                           S[c(TRUE, FALSE, NA)] 34425.445 37314.3870 43297.21237 38535.2235 40945.0395  91239.309   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  3212.309  3525.3850  5126.78145  3869.2520  5419.3185  58249.643   100
    ##                                                S[integer(0L), ]   140.179   251.3505   298.10649   305.6755   341.2020    668.997   100
    ##                                  as(S, "matrix")[integer(0L), ]  1442.995  1736.2270  3010.78908  1898.4435  3121.8220  56686.805   100
    ##                                                         S[1L, ]    29.930    60.0855    77.96929    78.4330    93.5825    233.372   100
    ##                                           as(S, "matrix")[1L, ]  1414.172  1738.6050  3036.73757  1912.2195  3018.0715  54079.000   100
    ##                                           S[1L, , drop = FALSE]   176.505   288.5170   359.53064   354.2400   390.9145   2773.199   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1620.115  1741.0650  2688.33761  1902.7280  3584.4455   8673.386   100
    ##                                                 S[c(0:1, NA), ]   164.861   216.1110   331.21358   349.3405   393.1695    663.052   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]  1616.384  1726.2435  3146.39658  1937.1270  3480.4490  52775.405   100
    ##                                                 S[c(0:2, NA), ]  1764.394  1956.4380  2844.76204  2186.8170  3659.2500   8645.916   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1545.782  1707.0760  2704.66545  1873.1875  3371.9015  10598.131   100
    ##                                          S[-seq_len(nrow(S)), ]   152.028   260.6165   334.67234   315.0235   341.2635   3668.762   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1502.363  1729.9335  3997.27409  1894.5280  3157.1845  56251.877   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    38.130    72.9390    91.29839    94.1360   106.7435    256.988   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1511.137  1710.9710  2989.99019  1824.3770  3026.5585  58203.641   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1375.714  1915.0280  3935.68553  2141.4915  3415.0950  57510.987   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1494.163  1712.3855  2469.40868  1798.3215  2999.2730   8981.624   100
    ##                                                S[logical(0L), ]   145.222   247.4760   292.93188   306.0650   332.0590    664.897   100
    ##                                  as(S, "matrix")[logical(0L), ]  1585.306  1686.0635  2441.85135  1831.6750  3123.4210   7174.303   100
    ##                                                       S[TRUE, ]     5.248    12.4845    27.07107    24.0055    37.1665    108.404   100
    ##                                         as(S, "matrix")[TRUE, ]  3068.481  3454.7215  5467.28686  4521.4390  6012.6500  60430.556   100
    ##                                                      S[FALSE, ]   140.261   233.2695   283.98650   304.7325   337.5735    530.786   100
    ##                                        as(S, "matrix")[FALSE, ]  1578.213  1721.8975  3118.13200  1899.7555  3588.4225  50959.761   100
    ##                                                         S[NA, ]  1131.682  1292.3200  2258.07500  1344.2055  1488.4435  50110.897   100
    ##                                           as(S, "matrix")[NA, ]  2691.732  2922.8695  4254.56508  3554.4745  4992.3650  10430.277   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2881.931  3296.3795  5020.05435  3642.7065  5431.4955  57551.044   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  2524.821  2744.0685  4569.45820  2980.0850  4855.8760  55523.594   100
    ##                                              S[character(0L), ]   143.623   226.2790   290.47762   304.9990   340.9150    644.889   100
    ##                                as(S, "matrix")[character(0L), ]  1461.117  1692.8695  2637.15608  1889.2390  3201.1160   9743.076   100
    ##                                            S[rownames(S)[1L], ]    35.096    59.1835    82.13858    83.3530    99.9785    208.403   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1542.953  1716.0755  3077.34151  1882.9455  3308.0235  54528.934   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   182.860   305.4090   336.46855   358.6475   380.6645    605.488   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1462.470  1723.6400  3073.64208  1942.5390  3452.2820  57726.688   100
    ##                                           S[rownames(S)[1:2], ]  1645.822  1896.9060  3299.91985  2055.9245  3552.6295  54773.458   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1579.894  1739.8145  2701.14232  1864.4750  3300.1105   9929.216   100
    ##                                           S[matrix(0L, 0L, 2L)]    16.892    38.2940    56.06422    54.8170    68.9415    142.393   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1409.170  1700.1880  2497.55682  1804.8610  3293.8170   8029.112   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    26.486    63.2220    84.84663    86.2025   101.2905    158.465   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1398.346  1692.3980  2507.09588  1767.0795  2650.5270  12068.350   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   127.633   189.5840   218.59027   212.9540   233.1260    417.585   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1226.433  1777.0630  3357.80693  1993.3995  3433.7090  56919.972   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 29360.633 31630.9465 38316.45324 32588.9935 35163.5270  83531.063   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3928.702  4212.4425  6044.87149  5451.6675  6436.4260  55959.301   100
    ##                                                         S[NULL]     2.870     8.3435    18.44713    17.5480    26.0760     46.576   100
    ##                                           as(S, "matrix")[NULL]  1526.348  1691.8650  2805.53488  1781.8190  2489.8480  55191.699   100
    ##                                                       S[NULL, ]   139.728   169.9245   272.54996   296.1225   338.5985    480.438   100
    ##                                         as(S, "matrix")[NULL, ]  1499.206  1710.8275  2596.25571  1824.4795  3045.2750   9890.553   100
    ##                                                       S[, NULL]   142.516   263.9990   292.39970   299.3205   330.3165    597.288   100
    ##                                         as(S, "matrix")[, NULL]  1111.633  1706.5430  3191.19605  1892.1910  3215.3635  57164.537   100

The `"packedMatrix"` subset is slower in some cases because priority is
given to memory efficiency:

``` r
S <- newDsp(30000L)
i <- replace(logical(1000L), 1L, TRUE)
{S[i]; cat("done!\n")}
```

    ## done!

``` r
{as(S, "matrix")[i]; cat("done!\n")}
```

    ## Error: vector memory exhausted (limit reached?)

## Transpose tests

``` r
S <- newDsp(1000L)
identical(as(t(S), "matrix"), t(as(S, "matrix")))
```

    ## [1] TRUE

``` r
microbenchmark::microbenchmark(t(S), as(t(as(S, "dsyMatrix")), "dspMatrix"))
```

    ## Unit: milliseconds
    ##                                    expr      min       lq     mean   median       uq      max neval
    ##                                    t(S) 7.981962 8.175646 8.299750 8.243132 8.367362 8.865881   100
    ##  as(t(as(S, "dsyMatrix")), "dspMatrix") 5.007166 5.200050 5.450021 5.307614 5.562491 6.850977   100

``` r
S <- newDsp(30000L)
{t(S); cat("done!\n")}
```

    ## done!

``` r
{as(t(as(S, "dsyMatrix")), "dspMatrix"); cat("done!\n")}
```

    ## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': vector memory exhausted (limit reached?)
