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
    ##                                                             S[]    1.394    4.4075    12.46482   11.6235    18.2450    37.351   100
    ##                                               as(S, "matrix")[] 1851.150 2398.3770  4716.48420 3701.8080  4257.7475 56355.074   100
    ##                                                         S[NULL]    5.289    8.0360    17.64763   14.8625    23.7595    70.766   100
    ##                                           as(S, "matrix")[NULL] 1480.510 1701.8895  3297.85386 1945.2655  3519.2350 59889.233   100
    ##                                                  S[integer(0L)]    2.911    9.2455    18.20113   18.7575    25.9530    44.034   100
    ##                                    as(S, "matrix")[integer(0L)] 1460.133 1732.8445  3185.92222 2042.3740  3531.3915 52831.862   100
    ##                                                   S[c(0:6, NA)]    4.387   14.7805    28.77462   27.2855    38.7860   147.395   100
    ##                                     as(S, "matrix")[c(0:6, NA)] 1458.821 1715.9115  3013.34010 1868.3085  3301.8940 53043.545   100
    ##                                                       S[-(0:6)] 4493.682 4731.4410  7334.41005 5492.9750  6745.7300 62563.745   100
    ##                                         as(S, "matrix")[-(0:6)] 3411.733 5148.8825  6726.14184 5554.0445  6271.2575 54458.332   100
    ##                                                  S[logical(0L)]    3.116    8.9790    19.35118   18.5935    27.6135    69.249   100
    ##                                    as(S, "matrix")[logical(0L)] 1484.528 1728.3550  3108.41295 1918.0415  3337.1745 52695.250   100
    ##                                                         S[TRUE]  874.366 1142.3010  1592.34447 1226.7815  1437.1525  5195.725   100
    ##                                           as(S, "matrix")[TRUE] 3771.918 4463.9160  6494.10726 5553.5320  6046.0240 56571.636   100
    ##                                                        S[FALSE]    3.772    9.0405    21.90138   19.7825    27.4700    90.446   100
    ##                                          as(S, "matrix")[FALSE] 2193.623 2687.5910  3635.78119 3164.4620  4397.5575  6919.447   100
    ##                                                           S[NA]  907.453 1001.9785  2745.69046 1082.0310  2543.8040 57134.812   100
    ##                                             as(S, "matrix")[NA] 3465.443 4103.0955  7481.94527 5266.9420  6117.1180 72372.462   100
    ##                                           S[c(TRUE, FALSE, NA)] 7788.524 8438.5585 10771.77666 9733.8100 10646.7160 58002.454   100
    ##                             as(S, "matrix")[c(TRUE, FALSE, NA)] 3314.809 3666.5480  6714.14155 4947.9005  5526.0005 54452.756   100
    ##                                                           S[, ]    2.665    5.9655    13.90843   15.2520    18.9420    42.394   100
    ##                                             as(S, "matrix")[, ] 3175.409 3448.5510  4685.53576 4857.0445  5434.0375  8798.026   100
    ##                                                       S[NULL, ]  146.452  284.5810   302.04823  307.1925   338.6190   690.686   100
    ##                                         as(S, "matrix")[NULL, ] 1409.785 1711.6885  3802.57657 2044.0345  3544.5935 70609.298   100
    ##                                                S[integer(0L), ]  136.940  287.1230   310.74802  308.1970   338.7625   544.439   100
    ##                                  as(S, "matrix")[integer(0L), ] 1471.941 1756.2555  3110.08083 2002.5015  3365.6490 53586.057   100
    ##                                                         S[1L, ]    8.159   19.7825    33.57572   34.0300    46.1455    75.809   100
    ##                                           as(S, "matrix")[1L, ] 1526.348 1732.1475  2729.73982 2041.7385  3664.2725  6364.881   100
    ##                                           S[1L, , drop = FALSE]    9.020   22.7960    36.63227   36.6130    46.4120   103.525   100
    ##                             as(S, "matrix")[1L, , drop = FALSE] 1428.809 1728.9700  3560.40023 1915.5200  3537.9720 54629.138   100
    ##                                                 S[c(0:6, NA), ]   14.432   31.6930    47.28612   46.7400    58.4250   120.581   100
    ##                                   as(S, "matrix")[c(0:6, NA), ] 1563.412 1746.2720  3727.40020 1999.8775  3431.8640 52217.395   100
    ##                                  S[-seq_len(nrow(S))[-(1:6)], ]   21.853   46.8425    61.05433   60.8235    75.0505   134.931   100
    ##                    as(S, "matrix")[-seq_len(nrow(S))[-(1:6)], ] 1491.990 1747.4815  3668.87639 2118.6545  3542.2360 52155.731   100
    ##                                                S[logical(0L), ]  140.261  283.0640   311.04076  317.9140   346.7370   621.724   100
    ##                                  as(S, "matrix")[logical(0L), ] 1535.696 1720.5035  3062.61923 1921.6290  3458.4935 51153.773   100
    ##                                                       S[TRUE, ]    6.150   16.1540    29.66555   32.3080    38.9705    94.136   100
    ##                                         as(S, "matrix")[TRUE, ] 3081.642 3467.7800  5630.68539 4407.6640  5303.3500 54913.883   100
    ##                                                      S[FALSE, ]  151.290  283.3305   315.57495  312.6250   342.5755   621.560   100
    ##                                        as(S, "matrix")[FALSE, ] 1490.924 1711.7500  3092.78703 1961.6860  3389.6955 52981.512   100
    ##                                                         S[NA, ] 1101.219 1288.3635  1605.84331 1340.0030  1452.0765  4544.235   100
    ##                                           as(S, "matrix")[NA, ] 2585.542 2947.7565  4705.06652 4141.5125  4902.7390 55521.708   100
    ##                                         S[c(TRUE, FALSE, NA), ]  647.513  817.0480  1552.65401  908.4575  1019.3625 48210.137   100
    ##                           as(S, "matrix")[c(TRUE, FALSE, NA), ] 2470.865 2762.0470  3700.20408 3076.2505  4494.5635  7216.779   100
    ##                                              S[character(0L), ]  146.329  280.9935   301.75918  308.0535   338.9265   561.331   100
    ##                                as(S, "matrix")[character(0L), ] 1411.958 1706.3380  3587.18061 2049.0980  3400.1710 53304.100   100
    ##                                           S[rownames(S)[1:6], ]   28.290   49.7945    69.49254   64.2060    87.0020   189.256   100
    ##                             as(S, "matrix")[rownames(S)[1:6], ] 1391.417 1762.3850  3150.16366 1953.5065  3499.6165 50907.240   100
    ##                                           S[matrix(0L, 0L, 2L)]   19.680   47.8060    58.52955   60.2700    69.9050    97.416   100
    ##                             as(S, "matrix")[matrix(0L, 0L, 2L)] 1506.176 1680.7130  2472.49352 1848.7515  3097.1195  6513.055   100
    ##                                S[cbind(c(0:6, NA), c(NA, 6:0))]   27.552   64.5545    86.72361   86.4895   102.3155   189.666   100
    ##                  as(S, "matrix")[cbind(c(0:6, NA), c(NA, 6:0))] 1227.745 1746.0465  3720.59584 2108.5480  3817.1000 52484.469   100
    ##                S[cbind(c(rownames(S), NA), c(NA, colnames(S)))]  128.371  198.4810   221.89979  222.1790   240.0140   363.219   100
    ##  as(S, "matrix")[cbind(c(rownames(S), NA), c(NA, colnames(S)))] 1384.160 1795.0825  3078.42514 2067.0560  3294.7805 53655.839   100
    ##                     S[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 7743.137 8409.7355 10245.89918 9617.3700 10683.3495 55113.389   100
    ##       as(S, "matrix")[matrix(c(TRUE, FALSE), nrow(S), ncol(S))] 3870.482 4632.0570  7291.84344 5604.8025  6724.1230 55563.159   100

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
    ## <bytecode: 0x10f905d98>
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
    ##                                    expr      min        lq      mean    median        uq      max neval
    ##                                    t(L)  536.198  649.3375  753.1409  739.0455  838.9625 1394.164   100
    ##  as(t(as(L, "dtrMatrix")), "dtpMatrix") 4455.101 5307.0195 5507.8465 5433.5660 5534.3850 7281.682   100

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
    ## <bytecode: 0x10c310fb8>
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
    ## <bytecode: 0x10c3bfdb8>
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
