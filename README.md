packedMatrix
================
2021-12-23

Based on discussion on R-devel beginning
[here](https://stat.ethz.ch/pipermail/r-devel/2021-November/081261.html)
and continued privately with Martin Mächler.

This machinery was added to `Matrix` in
[`r3424`](https://r-forge.r-project.org/scm/viewvc.php?view=rev&root=matrix&revision=3424).

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
    ##                                                            expr      min        lq        mean     median         uq       max neval
    ##                                                             S[]    1.476    3.9565    11.09706    11.4595    15.2930    32.759   100
    ##                                               as(S, "matrix")[] 1792.643 2116.9530  3709.40243  2724.0810  4013.4285 52053.600   100
    ##                                                         S[NULL]    2.624    7.9540    16.35736    15.3955    22.0375    60.106   100
    ##                                           as(S, "matrix")[NULL] 1107.328 1593.2190  3005.35576  1728.8675  3404.0045 54967.593   100
    ##                                                  S[integer(0L)]    3.239    7.8925    16.89856    14.6575    25.3380    47.109   100
    ##                                    as(S, "matrix")[integer(0L)] 1338.978 1630.9390  2804.70381  1753.6315  3057.8415 49145.470   100
    ##                                                   S[c(0:6, NA)]    3.936   13.4480    23.07152    22.6730    30.6885    51.865   100
    ##                                     as(S, "matrix")[c(0:6, NA)] 1399.494 1632.7225  3501.08143  1732.6805  3527.0455 52206.202   100
    ##                                                       S[-(0:6)] 4647.227 4967.7035  6952.08464  5269.9350  6838.3490 53817.010   100
    ##                                         as(S, "matrix")[-(0:6)] 3553.470 3856.1935  5735.59619  5295.8265  5936.1850 49335.546   100
    ##                                                  S[logical(0L)]    2.665    8.2000    16.17614    14.1860    24.7230    48.462   100
    ##                                    as(S, "matrix")[logical(0L)] 1360.831 1660.8895  3502.28355  1784.3610  3267.9050 53634.683   100
    ##                                                         S[TRUE]  788.471 1085.5160  1984.08020  1176.8435  1326.7395 51734.046   100
    ##                                           as(S, "matrix")[TRUE] 3443.918 3966.2375  7363.07930  5436.9280  6295.1195 55704.896   100
    ##                                                        S[FALSE]    2.952    8.0360    17.27535    15.0470    25.1740    68.183   100
    ##                                          as(S, "matrix")[FALSE] 2395.507 2634.5985  3325.12009  2733.5110  4074.1290  6107.155   100
    ##                                                           S[NA]  626.275  973.5860  1443.53989  1006.7345  1126.7415  6053.855   100
    ##                                             as(S, "matrix")[NA] 3255.646 3768.9250  5903.40919  4326.3815  5923.5365 55901.122   100
    ##                                           S[c(TRUE, FALSE, NA)] 7653.347 8117.8565 10351.58406  9357.7990 10398.1330 59595.591   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)] 3224.035 3497.5870  4686.21226  4673.9180  5551.6870  8011.769   100
    ##                                                           S[, ]    1.681    4.1410    10.38694     9.5530    15.6210    26.896   100
    ##                                             as(S, "matrix")[, ] 2960.405 3312.6155  4651.66156  4729.2270  5400.7660  9078.138   100
    ##                                                       S[NULL, ]    5.084   16.6460    27.64056    27.7365    36.9615    71.832   100
    ##                                         as(S, "matrix")[NULL, ] 1390.351 1641.7425  3045.94330  1847.5215  3408.9040 51014.291   100
    ##                                                S[integer(0L), ]    6.601   17.5685    29.90458    29.1920    39.8315   102.746   100
    ##                                  as(S, "matrix")[integer(0L), ] 1421.224 1626.4495  2815.29944  1726.4280  3019.0965 52159.503   100
    ##                                                         S[1L, ]    6.437   16.6665    30.39002    30.0735    40.8565    70.151   100
    ##                                           as(S, "matrix")[1L, ] 1399.822 1643.8950  2280.02763  1721.6310  2829.1025  6825.885   100
    ##                                           S[1L, , drop = FALSE]    8.118   21.9145    33.62082    31.4265    43.0705   110.331   100
    ##                             as(S, "matrix")[1L, , drop = FALSE] 1353.984 1631.1850  2386.11595  1741.8440  3259.5205  5794.038   100
    ##                                                 S[c(0:6, NA), ]   13.161   31.5495    44.77692    45.1615    55.7805   124.927   100
    ##                                   as(S, "matrix")[c(0:6, NA), ] 1466.857 1659.2905  3900.60060  1755.4150  3162.1045 53611.190   100
    ##                                  S[-seq_len(nrow(S))[-(1:6)], ]   21.648   44.5465    56.58492    56.1905    68.8800   127.100   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:6)], ] 1229.549 1664.3745  2412.15628  1731.3480  3189.0005  7270.940   100
    ##                                                S[logical(0L), ]  153.627  283.8020   302.07365   300.6325   321.4195   491.590   100
    ##                                  as(S, "matrix")[logical(0L), ] 1398.715 1656.3385  2876.50793  1769.7855  3110.0755 51410.679   100
    ##                                                       S[TRUE, ]    6.150   15.5185    24.06331    23.5955    31.5495    73.390   100
    ##                                         as(S, "matrix")[TRUE, ] 2958.109 3273.3580  4942.44545  4434.1910  5381.1475 50971.241   100
    ##                                                      S[FALSE, ]  153.545  281.4035   294.75925   300.6120   318.1600   661.986   100
    ##                                        as(S, "matrix")[FALSE, ] 1493.056 1650.0655  2452.12267  1750.3515  3316.8385  5457.633   100
    ##                                                         S[NA, ] 1111.838 1271.6355  1665.58154  1320.7535  1419.2970  4203.525   100
    ##                                           as(S, "matrix")[NA, ] 2411.866 2776.2330  4399.72394  3764.3330  4763.8925 52036.093   100
    ##                                         S[c(TRUE, FALSE, NA), ]  581.421  788.7375  1069.06926   859.7700   934.7180  4595.854   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ] 2297.640 2585.9520  3772.04059  2785.0275  4084.6045 48244.454   100
    ##                                              S[character(0L), ]    6.396   14.3705    28.15511    27.4085    37.9865    91.430   100
    ##                                as(S, "matrix")[character(0L), ] 1251.197 1645.4120  2801.40454  1731.2865  3261.8985 50581.167   100
    ##                                           S[rownames(S)[1:6], ]   22.427   45.6330    80.00494    61.4795    75.0505  1925.155   100
    ##                             as(S, "matrix")[rownames(S)[1:6], ] 1476.697 1686.4940  4631.42232  1825.5455  3700.3935 52205.300   100
    ##                                           S[matrix(0L, 0L, 2L)]   15.990   41.8405    52.18480    54.1815    65.6205   122.385   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)] 1454.967 1626.8390  2837.61410  1736.5140  3185.1260 50963.861   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]   22.837   49.8560    61.23678    62.3815    71.0735   177.120   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))] 1459.846 1661.9145  2837.76170  1754.3490  3121.9040 51067.263   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]   90.159  147.8050   168.37552   166.3165   178.7395   364.080   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))] 1256.404 1734.8330  2383.26440  1832.9255  3110.2600  6031.838   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 8710.245 9460.0940 14437.38699 10872.8310 11635.9025 59630.892   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 3608.041 4127.5520  5262.68169  5294.4120  6037.9060  9587.604   100

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
    ## <bytecode: 0x12d07d040>
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
    ##                                    expr      min       lq      mean    median        uq      max neval
    ##                                    t(L)  584.086  609.096  732.6667  767.2945  815.3465 1221.431   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4490.976 5013.357 5232.0596 5137.9765 5264.5025 7134.615   100

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
    ## <bytecode: 0x129d37818>
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
    ## <bytecode: 0x129de4618>
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
