packedMatrix
================
2021-12-13

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
    ##                                                  S[integer(0L)]     3.895     8.3435    19.06869    18.7165    25.0715     50.389   100
    ##                                    as(S, "matrix")[integer(0L)]  1223.645  1660.9510  2930.44425  1764.7425  2135.9565  59765.700   100
    ##                                                   S[c(0:6, NA)]    16.359    37.5560    56.35942    58.5480    68.3470    147.641   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1182.522  1633.4400  2813.14571  1761.9545  2477.2405  52300.543   100
    ##                                                       S[-(0:6)] 47633.185 50722.2685 63076.11175 52390.1485 56674.9765 120609.413   100
    ##                                         as(S, "matrix")[-(0:6)]  3438.096  4055.3100  5730.53843  5053.8035  5710.9310  59975.169   100
    ##                                                  S[logical(0L)]     2.747     8.9995    17.72922    18.0400    22.8575     69.085   100
    ##                                    as(S, "matrix")[logical(0L)]  1354.845  1664.4975  2444.97637  1794.0985  3057.3700   9616.550   100
    ##                                                         S[TRUE]   662.027  1108.4965  1616.21180  1184.0800  1389.1825   6544.953   100
    ##                                           as(S, "matrix")[TRUE]  3671.427  4085.7320  5581.28695  4841.5670  5816.8750  51219.947   100
    ##                                                        S[FALSE]     4.756     7.9335    16.78294    15.5595    23.7595     63.837   100
    ##                                          as(S, "matrix")[FALSE]  1992.559  2631.0725  3348.96651  2753.4780  4164.1240   7788.442   100
    ##                                                           S[NA]   637.099   974.8570  1947.24826  1010.6090  1144.0845  49208.077   100
    ##                                             as(S, "matrix")[NA]  2951.672  3883.4790  5130.73959  5009.7285  6039.8740  10916.127   100
    ##                                           S[c(TRUE, FALSE, NA)] 33392.081 36681.8185 45613.07851 38114.3995 41225.4385 129927.565   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  2925.473  3540.0015  5457.66416  4083.1900  5896.7020  56696.686   100
    ##                                                S[integer(0L), ]   142.106   255.8605   290.54322   302.9490   332.9610    751.612   100
    ##                                  as(S, "matrix")[integer(0L), ]  1104.417  1682.5990  3414.32420  1847.9930  3345.4975  51646.675   100
    ##                                                         S[1L, ]    26.363    54.4890    76.28296    73.5130    95.7555    185.074   100
    ##                                           as(S, "matrix")[1L, ]  1302.201  1636.8840  2279.62296  1765.2345  2070.3770   7921.938   100
    ##                                           S[1L, , drop = FALSE]   166.419   222.1995   319.35966   332.5920   373.2025    568.752   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1211.222  1636.2895  2506.35173  1769.1910  3193.6745   7082.586   100
    ##                                                 S[c(0:1, NA), ]   171.872   273.9415   328.93316   338.6805   369.9430    827.585   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]  1314.583  1634.4855  2333.49614  1769.2115  2582.4260   6432.203   100
    ##                                                 S[c(0:2, NA), ]  1471.859  1839.2190  3719.98781  2065.6620  3486.6810  54068.504   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1412.409  1696.1290  3639.97098  1829.9735  2764.8760  58275.965   100
    ##                                          S[-seq_len(nrow(S)), ]   145.058   232.6545   293.74737   293.3960   343.5390    646.242   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1269.032  1651.6235  2259.11066  1753.6520  2076.8755   8567.893   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    39.032    71.6885    88.17255    85.7925   103.8530    178.596   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1208.680  1659.9875  2721.60296  1783.2540  2386.8560  52372.539   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1616.056  1898.1155  3918.63076  2070.7665  3662.3045  54886.495   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1141.153  1702.2380  3441.97706  1864.3110  3802.9960  56019.325   100
    ##                                                S[logical(0L), ]   137.801   249.0955   285.29563   296.3275   324.2895    608.563   100
    ##                                  as(S, "matrix")[logical(0L), ]  1377.477  1677.1460  3001.41484  1787.2720  2987.5265  53954.565   100
    ##                                                       S[TRUE, ]     4.961    13.6530    27.35397    22.3450    34.8090    162.483   100
    ##                                         as(S, "matrix")[TRUE, ]  2968.441  3371.1840  4321.91537  3678.3560  4865.6750  10503.708   100
    ##                                                      S[FALSE, ]   141.901   221.7075   276.98206   294.2570   323.1620    529.310   100
    ##                                        as(S, "matrix")[FALSE, ]  1356.690  1660.5000  2795.50259  1801.4580  2847.9420  51784.722   100
    ##                                                         S[NA, ]  1070.633  1245.3135  2286.88775  1333.7300  1453.0605  61002.260   100
    ##                                           as(S, "matrix")[NA, ]  2621.540  2806.7370  3908.66120  3123.2570  4835.6630   8876.705   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2892.427  3257.9215  5377.02741  3645.9045  5356.5475  57134.853   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  1988.213  2646.2630  5000.96311  2824.7770  4216.9935  58780.675   100
    ##                                              S[character(0L), ]   146.575   259.3455   283.62119   292.1455   327.3850    627.587   100
    ##                                as(S, "matrix")[character(0L), ]  1315.157  1667.1830  2445.11741  1819.3545  3158.6195   7993.114   100
    ##                                            S[rownames(S)[1L], ]    33.702    62.0125    86.76707    82.7585   103.9145    266.828   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1523.847  1738.5435  3971.31494  1909.7185  3872.3475  52095.953   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   183.270   304.9375   340.96338   346.4705   383.3705    642.511   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1146.196  1678.0275  2801.10606  1794.8365  2482.1400  55653.605   100
    ##                                           S[rownames(S)[1:2], ]  1529.915  1874.5815  3202.62849  2002.9320  3169.9765  57451.742   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1242.300  1664.8050  2415.51008  1780.6710  2243.8890   9755.909   100
    ##                                           S[matrix(0L, 0L, 2L)]    16.933    40.3645    55.19953    56.6825    69.4130    109.675   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1337.420  1675.5265  2921.69526  1776.4070  2945.1735  54668.498   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    28.413    65.7845    84.79866    83.0865    96.2680    171.134   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1290.311  1684.5260  2996.72157  1798.5060  2789.3325  51695.219   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   125.337   186.8165   211.87242   212.0930   230.8710    336.446   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1296.461  1762.4055  2460.27511  1900.3910  2829.6970   9252.019   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 28900.408 31521.2510 39450.95112 33294.7470 34512.1395 101331.213   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3522.392  4271.0725  6611.87320  5419.1340  6478.2460  53338.171   100
    ##                                                         S[NULL]     4.141     8.4460    18.09207    17.5685    24.4360     45.346   100
    ##                                           as(S, "matrix")[NULL]  1330.368  1632.3740  2249.98775  1798.8955  2322.6500   6677.629   100
    ##                                                       S[NULL, ]   137.883   162.7085   249.95896   263.5070   307.5205    548.252   100
    ##                                         as(S, "matrix")[NULL, ]  1446.316  1653.7145  2272.41065  1785.3655  2528.3265   6693.824   100
    ##                                                       S[, NULL]   138.785   203.3190   276.85988   289.8495   316.9505    671.211   100
    ##                                         as(S, "matrix")[, NULL]  1121.678  1655.6620  2263.32833  1781.9830  2223.6760   7944.119   100

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

    ## Unit: microseconds
    ##                                    expr      min       lq      mean   median       uq      max neval
    ##                                    t(S)  623.159  648.374  685.5077  666.045  698.763  910.856   100
    ##  as(t(as(S, "dsyMatrix")), "dspMatrix") 4940.582 5128.731 5367.2514 5239.636 5430.040 6989.762   100

``` r
S <- newDsp(30000L)
{t(S); cat("done!\n")}
```

    ## done!

``` r
{as(t(as(S, "dsyMatrix")), "dspMatrix"); cat("done!\n")}
```

    ## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': vector memory exhausted (limit reached?)
