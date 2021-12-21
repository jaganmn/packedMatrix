packedMatrix
================
2021-12-21

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
    ##                                                            expr      min         lq        mean     median         uq       max neval
    ##                                                             S[]    1.927     4.5920    17.81327    15.2930    23.3495    54.530   100
    ##                                               as(S, "matrix")[] 2219.740  2641.4455  5177.32297  4167.5475  5882.5980 64189.928   100
    ##                                                         S[NULL]    2.419     9.6965    20.85834    17.5275    27.0805    70.069   100
    ##                                           as(S, "matrix")[NULL] 1615.359  1926.8975  4381.48017  2924.8170  4090.3035 60598.492   100
    ##                                                  S[integer(0L)]    3.116     9.1635    23.84560    21.3815    30.0530   111.397   100
    ##                                    as(S, "matrix")[integer(0L)] 1752.504  1917.6110  3347.42450  3003.7625  3892.4375 12882.815   100
    ##                                                   S[c(0:6, NA)]    4.264    15.5185    32.11694    31.6110    40.0365    99.343   100
    ##                                     as(S, "matrix")[c(0:6, NA)] 1610.849  1897.3775  3721.25184  2941.4015  3783.3980 63111.382   100
    ##                                                       S[-(0:6)] 5415.526  5515.5045  8627.81737  7043.0825  9596.4805 70577.974   100
    ##                                         as(S, "matrix")[-(0:6)] 4046.208  5624.0520  8862.60141  6580.8075  8500.1815 69278.479   100
    ##                                                  S[logical(0L)]    4.182     9.3070    26.10675    23.7800    31.5495   110.700   100
    ##                                    as(S, "matrix")[logical(0L)] 1619.008  1919.4355  5209.99423  3065.4265  4257.4195 74022.179   100
    ##                                                         S[TRUE] 1018.604  1243.6735  2218.67441  1376.8620  2321.0305  9553.287   100
    ##                                           as(S, "matrix")[TRUE] 4298.071  5191.7685  7906.52077  6289.2975  7967.8580 61929.229   100
    ##                                                        S[FALSE]    4.264     9.6145    22.23225    18.4910    29.4995    87.576   100
    ##                                          as(S, "matrix")[FALSE] 2718.464  2980.8435  5573.42725  3961.5020  5006.1615 68514.362   100
    ##                                                           S[NA] 1060.465  1117.5985  2518.92233  1141.8090  1979.3365 70038.742   100
    ##                                             as(S, "matrix")[NA] 4024.150  5566.6725  8663.72394  6542.7185  7920.0520 76400.384   100
    ##                                           S[c(TRUE, FALSE, NA)] 8703.521 11043.2270 14716.28539 12230.4025 14565.5575 78580.231   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)] 3650.845  4206.2310  7190.57262  5556.8530  7183.5690 70019.185   100
    ##                                                           S[, ]    1.640     4.4075    11.94699    10.4755    16.8715    54.571   100
    ##                                             as(S, "matrix")[, ] 3318.294  3971.5470  5607.07431  5401.8525  6175.5020 11414.441   100
    ##                                                       S[NULL, ]    5.207    15.0880    35.59005    28.2695    41.9225   207.214   100
    ##                                         as(S, "matrix")[NULL, ] 1649.061  1938.0905  4698.29988  3054.3975  4550.4875 62419.671   100
    ##                                                S[integer(0L), ]    5.904    19.8440    44.07746    32.7590    53.5460   206.107   100
    ##                                  as(S, "matrix")[integer(0L), ] 1739.876  1942.5800  3822.66862  3061.6135  3829.3795 56549.086   100
    ##                                                         S[1L, ]    7.380    16.9125    40.38377    32.9435    47.3140   263.302   100
    ##                                           as(S, "matrix")[1L, ] 1699.860  1987.9875  4686.99700  3140.3950  4227.5100 69776.178   100
    ##                                           S[1L, , drop = FALSE]    9.676    19.7825    44.49771    37.7200    54.6530   150.388   100
    ##                             as(S, "matrix")[1L, , drop = FALSE] 1545.003  1915.1510  3485.51414  3047.0380  4532.2630 12148.546   100
    ##                                                 S[c(0:6, NA), ]   13.366    29.4380    69.37364    46.5145    67.0965  1387.481   100
    ##                                   as(S, "matrix")[c(0:6, NA), ] 1666.609  1947.1925  3810.64414  3028.0140  4126.4655 57020.258   100
    ##                                  S[-seq_len(nrow(S))[-(1:6)], ]   27.552    43.6855    93.66368    61.6230    82.8610  1780.138   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:6)], ] 1729.995  1923.3305  4391.16683  2218.6125  3726.7770 86434.929   100
    ##                                                S[logical(0L), ]  155.554   312.3175   401.43469   337.8810   367.2165  1634.793   100
    ##                                  as(S, "matrix")[logical(0L), ] 1605.396  1911.9325  4009.27233  3038.5100  4471.5830 58499.948   100
    ##                                                       S[TRUE, ]    6.437    16.6255    34.10380    30.4835    39.9545   204.426   100
    ##                                         as(S, "matrix")[TRUE, ] 3461.384  3784.7510  7893.47170  5783.1730  7799.6760 70492.858   100
    ##                                                      S[FALSE, ]  160.802   308.5250   425.35573   331.4850   371.8290  2135.936   100
    ##                                        as(S, "matrix")[FALSE, ] 1635.121  1949.0580  5268.65252  3106.0165  4000.5135 68636.624   100
    ##                                                         S[NA, ] 1285.760  1439.1615  2378.90528  1495.2700  2903.6405  9128.445   100
    ##                                           as(S, "matrix")[NA, ] 2956.551  3336.2725  6588.06983  5175.9015  6916.5770 68375.290   100
    ##                                         S[c(TRUE, FALSE, NA), ]  820.943   894.9890  1742.50820  1020.3875  2550.2820  7183.692   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ] 2758.562  3022.0690  6537.28969  4690.5845  5576.6560 71641.104   100
    ##                                              S[character(0L), ]    7.093    21.5455    49.01468    39.4215    55.5140   452.681   100
    ##                                as(S, "matrix")[character(0L), ] 1663.903  1913.6955  4300.04638  2905.1370  4065.1910 57402.583   100
    ##                                           S[rownames(S)[1:6], ]   27.839    48.5235   138.34917    70.6635    93.6850  4950.668   100
    ##                             as(S, "matrix")[rownames(S)[1:6], ] 1662.140  1920.7065  4257.10175  2755.0155  3614.4985 66693.839   100
    ##                                           S[matrix(0L, 0L, 2L)]   18.819    45.7560    69.59627    61.0900    79.8680   240.096   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)] 1673.292  1923.1255  3185.05958  2104.3865  3946.8445 12089.711   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]   19.762    53.8125    84.74290    72.6315    87.5760   303.031   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))] 1736.555  1901.0880  3773.71011  2990.8885  3942.0065 64182.958   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  100.122   155.9230   257.02490   179.8055   327.4260  2177.920   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))] 1717.039  2003.7110  3460.37335  3126.0245  4182.9635 11913.124   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 9378.586 11754.8640 20600.98013 14128.7230 16949.7280 81360.974   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 4347.886  5945.6150  9388.88602  6932.3415  9509.7040 73731.899   100

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
    ## <bytecode: 0x133867ef8>
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
    ##                                    expr      min        lq     mean   median       uq       max neval
    ##                                    t(L)  651.121  684.6385 1085.222  855.260  913.808  4554.157   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 5101.507 5359.4790 6237.843 5522.762 6430.215 13670.671   100

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
    ## <bytecode: 0x1315db890>
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
    ## <bytecode: 0x131688690>
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
