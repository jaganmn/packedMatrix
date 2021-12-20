packedMatrix
================
2021-12-20

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
    ##                                                             S[]    2.091     5.5965    17.75997    16.0720    21.2175    67.035   100
    ##                                               as(S, "matrix")[] 2177.387  2575.4355  5442.39412  3930.6700  4890.7670 68035.687   100
    ##                                                         S[NULL]    2.665     8.9380    20.02891    17.9785    25.2150    70.397   100
    ##                                           as(S, "matrix")[NULL] 1571.612  1926.8770  3348.42572  2946.7930  4280.8715 10081.900   100
    ##                                                  S[integer(0L)]    3.854    10.7830    25.06166    22.0785    32.0620   118.408   100
    ##                                    as(S, "matrix")[integer(0L)] 1693.915  1936.5530  4105.22053  2396.5935  3562.9820 65749.855   100
    ##                                                   S[c(0:6, NA)]    5.002    16.9330    35.74626    31.7340    42.7220   107.994   100
    ##                                     as(S, "matrix")[c(0:6, NA)] 1792.192  1933.3550  3382.80873  2743.5150  3933.3350 12161.297   100
    ##                                                       S[-(0:6)] 5336.765  5947.2550  8533.69531  7988.4195  9984.5660 22997.310   100
    ##                                         as(S, "matrix")[-(0:6)] 4117.753  5435.5955  7548.90565  6361.0885  8624.5960 58698.839   100
    ##                                                  S[logical(0L)]    3.198    11.1930    23.80009    23.2470    30.3605    77.531   100
    ##                                    as(S, "matrix")[logical(0L)] 1661.853  2035.3015  3788.28069  3439.9410  4344.2575 10874.922   100
    ##                                                         S[TRUE] 1122.416  1285.4320  2498.25710  1382.2945  2193.7870 62166.578   100
    ##                                           as(S, "matrix")[TRUE] 4212.709  5794.3660  8466.66195  6843.9660  8587.6140 69642.682   100
    ##                                                        S[FALSE]    3.813    13.2430    26.10101    24.3745    32.4720   132.102   100
    ##                                          as(S, "matrix")[FALSE] 2757.373  3016.3290  4457.41709  3891.9045  5069.5065 11842.235   100
    ##                                                           S[NA] 1075.717  1133.3220  2400.69637  1176.1670  2106.2110 59041.476   100
    ##                                             as(S, "matrix")[NA] 3900.781  5548.0175  9209.91610  6758.4400  8345.5295 67783.906   100
    ##                                           S[c(TRUE, FALSE, NA)] 8720.044 10378.6580 14938.10933 12496.7795 15377.2960 86321.523   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)] 3629.935  4033.6825  6472.16283  5214.1135  6835.0690 77534.485   100
    ##                                                           S[, ]    1.722     5.0635    16.46806    13.2020    21.4840    57.359   100
    ##                                             as(S, "matrix")[, ] 3535.389  3777.5760  6320.65266  5280.0005  6536.2200 65670.192   100
    ##                                                       S[NULL, ]    5.863    17.7325    40.79049    33.6815    50.4915   258.095   100
    ##                                         as(S, "matrix")[NULL, ] 1501.502  1954.9210  4598.08645  2656.8820  3757.1990 73696.639   100
    ##                                                S[integer(0L), ]    6.191    14.3500    35.78029    34.7885    47.6420   115.866   100
    ##                                  as(S, "matrix")[integer(0L), ] 1784.935  1971.8745  5161.98815  3123.7695  4152.8490 65971.419   100
    ##                                                         S[1L, ]    7.626    24.3745    42.46001    38.8065    52.4390   127.469   100
    ##                                           as(S, "matrix")[1L, ] 1727.248  1953.4860  4043.53849  3283.9155  4217.6905 57049.860   100
    ##                                           S[1L, , drop = FALSE]    8.364    28.2900    50.38654    42.0865    55.0220   405.818   100
    ##                             as(S, "matrix")[1L, , drop = FALSE] 1752.299  1983.8055  4600.27995  2867.2120  3748.4660 68832.194   100
    ##                                                 S[c(0:6, NA), ]   15.211    33.1485   684.92058    55.9240    81.5490 59234.586   100
    ##                                   as(S, "matrix")[c(0:6, NA), ] 1787.026  1989.5660  3654.03808  3227.1715  4619.4085 11189.392   100
    ##                                  S[-seq_len(nrow(S))[-(1:6)], ]   23.493    43.1935   159.34076    66.5020   107.9735  2397.106   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:6)], ] 1754.185  1969.1685  4160.83621  2836.5645  4441.0585 72363.893   100
    ##                                                S[logical(0L), ]  155.677   307.5615   415.47432   341.4890   394.2560  1136.643   100
    ##                                  as(S, "matrix")[logical(0L), ] 1683.583  1969.1890  3773.73430  2915.0385  3861.6055 62746.769   100
    ##                                                       S[TRUE, ]    8.528    17.6505    34.54291    32.2055    41.7175   172.364   100
    ##                                         as(S, "matrix")[TRUE, ] 3600.251  4113.2840  6942.96542  5917.0175  7304.8880 66788.180   100
    ##                                                      S[FALSE, ]  152.930   317.3195   447.58060   342.2885   365.5560  5760.500   100
    ##                                        as(S, "matrix")[FALSE, ] 1559.025  1997.4585  3559.92709  3223.4610  4208.2195  9307.041   100
    ##                                                         S[NA, ] 1308.679  1452.6095  2287.30267  1595.7200  2682.3430  7581.474   100
    ##                                           as(S, "matrix")[NA, ] 2989.351  3306.0350  6509.63642  4832.4240  6357.0705 74680.270   100
    ##                                         S[c(TRUE, FALSE, NA), ]  842.017   906.1000  1640.53792  1044.3315  1819.7850  7822.308   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ] 2568.650  3060.7935  6075.14343  4698.1490  5848.0965 66682.851   100
    ##                                              S[character(0L), ]   10.824    28.9050    53.52878    45.8380    62.6070   170.396   100
    ##                                as(S, "matrix")[character(0L), ] 1744.386  1938.6850  3678.61963  2154.6730  3817.1000 66426.970   100
    ##                                           S[rownames(S)[1:6], ]   26.158    51.2090   113.99066    73.0825    97.5185  2008.016   100
    ##                             as(S, "matrix")[rownames(S)[1:6], ] 1578.623  2000.8000  3600.21861  3059.6455  4907.3720  8848.005   100
    ##                                           S[matrix(0L, 0L, 2L)]   17.671    55.1040    81.85691    70.9915    89.9130   258.710   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)] 1745.944  1962.4445  3421.74397  2908.3555  4077.1015 13629.220   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]   20.910    54.5505    82.00041    69.9870    92.0655   227.386   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))] 1681.246  1965.1300  3339.08469  2811.6160  4213.5905 10731.053   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  101.352   162.8315   354.52495   183.3110   237.3080  7545.476   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))] 1758.859  2023.6575  5118.33176  2920.4300  4484.2110 69354.247   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 9560.913 11858.4915 19788.17071 13214.3820 16812.5010 79532.251   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 4220.827  5733.6450  8366.77447  6726.6445  8513.1785 71999.485   100

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
    ## <bytecode: 0x12167ccf8>
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
    ##                                    expr      min       lq      mean   median        uq       max neval
    ##                                    t(L)  637.550  684.044  971.4811  863.091  907.3915  3785.325   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 5385.145 5621.284 6643.3546 5797.913 7005.7110 13478.545   100

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
    ## <bytecode: 0x1220dc690>
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
    ## <bytecode: 0x122189490>
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
