packedMatrix
================
2021-12-14

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
    ##                                                  S[integer(0L)]     5.002     7.9745    15.52998    14.9035    21.6070     41.123   100
    ##                                    as(S, "matrix")[integer(0L)]  1022.581  1592.9115  2435.70299  1742.1515  3067.4970   8663.341   100
    ##                                                   S[c(0:6, NA)]    14.186    38.4785    58.38031    54.9810    73.5130    134.972   100
    ##                                     as(S, "matrix")[c(0:6, NA)]  1074.979  1554.7405  2194.03095  1669.7045  2201.8435   6626.502   100
    ##                                                       S[-(0:6)] 45093.809 49920.1650 63339.73888 51659.6515 58721.5325 153498.506   100
    ##                                         as(S, "matrix")[-(0:6)]  2828.344  3907.7510  5531.53673  4714.5900  5897.9320  55137.784   100
    ##                                                  S[logical(0L)]     2.788     6.8880    15.27824    13.8785    21.5455     57.892   100
    ##                                    as(S, "matrix")[logical(0L)]  1228.647  1594.1825  2399.94320  1739.1790  2334.1095   9366.901   100
    ##                                                         S[TRUE]   744.109  1068.1320  3028.11035  1160.7715  1436.4965  50220.490   100
    ##                                           as(S, "matrix")[TRUE]  3299.229  3923.1670  5058.47094  4395.1385  5555.9715  14397.232   100
    ##                                                        S[FALSE]     4.592     8.1590    17.24050    15.9695    25.6455     39.319   100
    ##                                          as(S, "matrix")[FALSE]  2425.765  2568.5065  3696.70801  2688.6160  3083.4460  56082.055   100
    ##                                                           S[NA]   634.926   974.9595  1845.10742  1011.9825  1113.1500  52143.923   100
    ##                                             as(S, "matrix")[NA]  3471.634  3796.6820  4909.49990  4148.2160  5562.9210  12378.884   100
    ##                                           S[c(TRUE, FALSE, NA)] 33120.415 36299.4320 44590.19698 37296.6135 39538.5960  93124.448   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)]  2706.123  3424.9145  5360.67046  3636.5565  5282.5835  58244.313   100
    ##                                                S[integer(0L), ]   143.254   256.0860   276.20224   286.4875   307.3155    488.187   100
    ##                                  as(S, "matrix")[integer(0L), ]  1041.646  1567.5325  2246.25101  1707.5270  2505.0180   8510.452   100
    ##                                                         S[1L, ]    26.035    50.6965    68.03130    69.0850    82.5945    147.887   100
    ##                                           as(S, "matrix")[1L, ]  1412.409  1582.0875  2378.51455  1764.9270  2898.4540  14194.733   100
    ##                                           S[1L, , drop = FALSE]   168.387   301.6985   320.41746   330.8905   360.3900    615.082   100
    ##                             as(S, "matrix")[1L, , drop = FALSE]  1313.353  1598.7130  3580.22537  1757.7725  3512.7570  53489.174   100
    ##                                                 S[c(0:1, NA), ]   170.191   206.9475   303.61074   325.1300   368.6720    458.175   100
    ##                                   as(S, "matrix")[c(0:1, NA), ]   964.730  1580.6320  2641.54759  1695.0015  2196.9235  49501.022   100
    ##                                                 S[c(0:2, NA), ]  1558.984  1769.4165  3518.20016  1929.0910  3216.6140  52227.604   100
    ##                                   as(S, "matrix")[c(0:2, NA), ]  1218.643  1609.6190  2825.01316  1711.4630  2935.2515  56171.722   100
    ##                                          S[-seq_len(nrow(S)), ]   144.525   183.5365   267.62176   280.8295   314.6340    555.714   100
    ##                            as(S, "matrix")[-seq_len(nrow(S)), ]  1062.105  1590.2055  2787.36942  1741.5980  2214.0000  51596.368   100
    ##                                     S[-seq_len(nrow(S))[-1L], ]    39.565    64.2265    83.47108    85.5055   102.7050    189.051   100
    ##                       as(S, "matrix")[-seq_len(nrow(S))[-1L], ]  1240.783  1578.2950  2266.18398  1706.1125  2240.4450   7358.844   100
    ##                                  S[-seq_len(nrow(S))[-(1:2)], ]  1158.209  1804.2870  2602.29706  1979.4185  3198.1025   8087.291   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:2)], ]  1040.703  1627.5155  3226.22686  1764.6400  2551.2455  49769.941   100
    ##                                                S[logical(0L), ]   137.801   162.6880   248.86016   262.9535   314.5315    492.533   100
    ##                                  as(S, "matrix")[logical(0L), ]  1026.107  1570.9355  2579.50885  1676.4695  1937.4345  50621.552   100
    ##                                                       S[TRUE, ]     6.929    12.7305    23.13097    24.6615    31.4880     51.373   100
    ##                                         as(S, "matrix")[TRUE, ]  2769.058  3279.0160  4662.06203  3515.2990  4749.3990  55423.144   100
    ##                                                      S[FALSE, ]   140.835   170.0270   253.13113   269.6775   302.6825    517.748   100
    ##                                        as(S, "matrix")[FALSE, ]  1086.582  1589.5700  2415.52033  1711.1350  3065.7955   7270.407   100
    ##                                                         S[NA, ]   953.209  1240.7625  2258.91837  1318.6010  1425.5085  49795.197   100
    ##                                           as(S, "matrix")[NA, ]  2512.111  2790.2755  4498.52205  3211.6120  4727.6485  56413.130   100
    ##                                         S[c(TRUE, FALSE, NA), ]  2699.112  3172.5800  4343.77001  3694.2230  5292.1980  10927.156   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ]  1980.013  2602.1060  4002.62746  2826.8270  4201.6595  56642.853   100
    ##                                              S[character(0L), ]   145.714   255.7375   281.42400   291.6740   317.9960    589.744   100
    ##                                as(S, "matrix")[character(0L), ]  1309.868  1595.9455  2355.33602  1741.0035  2710.2435   9250.420   100
    ##                                            S[rownames(S)[1L], ]    32.759    58.3840    75.92749    75.3785    92.6805    134.070   100
    ##                              as(S, "matrix")[rownames(S)[1L], ]  1006.263  1606.7900  2883.21881  1736.2680  2957.3095  49076.344   100
    ##                              S[rownames(S)[1L], , drop = FALSE]   177.284   284.7245   330.40998   335.6670   368.7130    641.691   100
    ##                as(S, "matrix")[rownames(S)[1L], , drop = FALSE]  1136.889  1585.0600  2602.50657  1706.9120  2020.5005  50695.639   100
    ##                                           S[rownames(S)[1:2], ]  1421.306  1790.5725  2645.63570  1969.5580  3163.6010   9471.861   100
    ##                             as(S, "matrix")[rownames(S)[1:2], ]  1361.979  1594.8385  2651.43310  1693.0540  2051.1275  56875.815   100
    ##                                           S[matrix(0L, 0L, 2L)]    18.532    41.0205    54.57264    55.6575    66.9120    100.737   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)]  1218.151  1592.5835  2703.58510  1710.8275  2194.4430  55735.318   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]    26.240    62.1970    76.21613    74.8455    90.2410    162.975   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))]  1304.374  1592.7475  3448.96264  1741.2085  2987.9775  53358.220   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   129.109   189.7685   209.41898   207.8905   229.4975    289.337   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  1308.310  1658.3065  2603.86285  1800.8020  3486.2300  13669.687   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 28415.132 30725.4205 36614.08035 31914.1130 33942.0550  99550.788   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))]  3834.607  4091.5130  5625.39188  4446.2860  5872.5940  52062.907   100
    ##                                                         S[NULL]     4.674     8.5280    16.80221    14.6985    21.8325     47.519   100
    ##                                           as(S, "matrix")[NULL]  1452.548  1596.6015  2931.39791  1760.9910  2960.4460  50727.414   100
    ##                                                       S[NULL, ]   142.557   225.0490   264.92560   272.6295   307.9510    475.354   100
    ##                                         as(S, "matrix")[NULL, ]  1060.178  1564.3550  2658.70978  1712.8775  2106.0060  50371.247   100
    ##                                                       S[, NULL]   142.598   232.3880   271.81442   287.7175   312.3995    394.133   100
    ##                                         as(S, "matrix")[, NULL]   994.947  1586.3105  2812.39828  1700.2085  2491.0165  55356.396   100

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
    ## <bytecode: 0x11c9ab798>
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
    ##                                    t(L)  549.195  699.870  727.1034  718.156  751.694  855.916   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4235.874 5128.342 5197.0468 5202.306 5351.648 5885.878   100

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
    ## {
    ##     .Call(packedMatrix_diag_get, x, isTRUE(names[1L]))
    ## }
    ## <bytecode: 0x11890af68>
    ## <environment: namespace:Matrix>
    ## 
    ## Signatures:
    ##         x             
    ## target  "packedMatrix"
    ## defined "packedMatrix"

``` r
n <- 6L
U <- new("lspMatrix", uplo = "U", x = sample(c(TRUE, FALSE), choose(n + 1L, 2L), TRUE),
         Dim = c(n, n), Dimnames = list(letters[1:6], NULL))
U
```

    ## 6 x 6 Matrix of class "lspMatrix"
    ##       a     b     c     d     e     f
    ## a  TRUE  TRUE FALSE FALSE  TRUE FALSE
    ## b  TRUE  TRUE FALSE FALSE FALSE FALSE
    ## c FALSE FALSE FALSE FALSE FALSE  TRUE
    ## d FALSE FALSE FALSE  TRUE  TRUE  TRUE
    ## e  TRUE FALSE FALSE  TRUE FALSE FALSE
    ## f FALSE FALSE  TRUE  TRUE FALSE FALSE

``` r
diag(U)
```

    ##     a     b     c     d     e     f 
    ##  TRUE  TRUE FALSE  TRUE FALSE FALSE

``` r
diag(U, names = FALSE)
```

    ## [1]  TRUE  TRUE FALSE  TRUE FALSE FALSE

## `diag<-` method tests
