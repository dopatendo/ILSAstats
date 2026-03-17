# Standard Error for Estimates with Replicate Weights and Plausible Values

Calculates the standard error given a vector or list of previous
estimations.

## Usage

``` r
repse(er, e0, setup = NULL, method)

repsecomp(se)

pvse(PVse, PVe0, df = FALSE)
```

## Arguments

- er:

  a vector or a list containing any statistic of interest (e.g.,
  percent, mean, variance, regression coefficient). If it is a vector or
  list of `length==1`, the function estimates standard errors without
  plausible values. If it is a list with `length>1`, it estimates
  standard errors with plausible values.

- e0:

  a numeric vector or a vector containing any statistic of interest
  (e.g., percent, mean, variance, regression coefficient), computed
  using total weights. For scenarios without plausible values, `e0`
  should be a single value. For scenarios with plausible values, `e0`
  should be a vector of the same length as `er`.

- setup:

  an optional list produced by
  [`repsetup`](https://dopatendo.github.io/ILSAstats/reference/repsetup.md).

- method:

  a string indicating the name of the replication method. Available
  options are: `"JK2-full"`, `"JK2-half"`, `"FAY-0.5"`, and
  `"JK2-half-1PV"`.  
    
  Additionally, ILSA names can be used, defaulting into:

  - `"TIMSS"`, `"PIRLS"`, or `"LANA"` for `"JK2-full"`;

  - `"ICILS"`, `"ICCS"`, or `"CIVED"` for `"JK2-half"`;

  - `"PISA"` or `"TALIS"` for `"FAY-0.5"`;

  - and `"oldTIMSS"`, `"oldPIRLS"`, or `"RLII"` for `"JK2-half-1PV"`.

  Note that `"oldTIMSS"` and `"oldPIRLS"` refer to the method used for
  TIMSS and PIRLS before 2015, where within imputation variance is
  estimated using only 1 plausible value.

- se:

  a numeric vector with standard errors, used by `repsecomp()` to
  estimate a composite standard error.

- PVse:

  a numeric vector containing the standard errors of the estimates of
  each plausible value.

- PVe0:

  a numeric vector containing the point estimates of each plausible
  value.

- df:

  a logical value indicating if degrees should be calculated.

## Value

the standard error.

## Details

The standard errors are calculated using a modifier \\m\\, for JK2-full:
\\m = 0.5\\; for JK2-half: \\m = 1\\; and for FAY-0.5:
\\\frac{1}{R(1-0.5)^2}\\. Depending on the statistic, one of the
following formulas is used.

The standard error not involving plausible values is calculated by:

\$\$\sqrt{m\times \sum\_{r=1}^{R}(\varepsilon_r-\varepsilon_0)^2}.\$\$

The standard error involving plausibles values and replicate weights is
calculated by:

\$\$\sqrt{\left\[ \sum\_{p=1}^{P} \left( m\times
\sum\_{r=1}^{R}(\varepsilon\_{rp}-\varepsilon\_{0p})^2 \right)
\dfrac{1}{P}\right\]+ \left\[ \left(1+ \dfrac{1}{P} \right)
\dfrac{\sum\_{p=1}^{P}
(\varepsilon\_{0p}-\overline{\varepsilon}\_{0p})^{2}}{P-1}
\right\]}.\$\$

The standard error involving plausibles values without replicate weights
is calculated by:

\$\$\sqrt{ \dfrac{\sum\_{p=1}^{P} SE^2\_{\varepsilon\_{0P}}}{P}+ \left\[
\left(1+ \dfrac{1}{P} \right) \dfrac{\sum\_{p=1}^{P}
(\varepsilon\_{0p}-\overline{\varepsilon}\_{0p})^{2}}{P-1}
\right\]}.\$\$

The standard error of the difference of two statistics (\\a\\ and \\b\\)
from independent samples is calculated by:

\$\$\sqrt{SE_a^{2}+SE_b^{2}}.\$\$

The standard error of the difference of two statistics (\\a\\ and \\b\\)
from dependent samples not involving plausible values is calculated by:

\$\$\sqrt{m\times \sum\_{r=1}^R((a_r-b_r)-(a_0-b_0))^2}.\$\$

The standard error of the difference of two statistics (\\a\\ and \\b\\)
from dependent samples involving plausible values is calculated by:

\$\$\sqrt{\left\[ \sum\_{p=1}^{P} \left( m\times
\sum\_{r=1}^{R}((a\_{rp}-b\_{rp})-(a\_{0p}-b\_{0p}))^2 \right)
\dfrac{1}{P}\right\]+ \left\[ \left(1+ \dfrac{1}{P} \right)
\dfrac{\sum\_{p=1}^{P} \left((a\_{0p}-b\_{0p})- (
\overline{a}\_{0p}-\overline{b}\_{0p}) \right)^{2}}{P-1} \right\]}.\$\$

The standard error of a composite estimate is calculated by:

\$\$\sqrt{\dfrac{\sum\_{c=1}^CSE^2\_{\varepsilon_c}}{C^{2}}}.\$\$

The standard error of the difference between an element (\\a\\) of the
composite and the composite is calculated by:

\$\$\sqrt{\dfrac{\sum\_{c=1}^CSE^2\_{\varepsilon_c}}{C^{2}}+\left(\dfrac{(C-1)^2-1}{C^2}\right)SE^2_a}.\$\$

Where \\\varepsilon\\ represents a statistic of interest, the subindex
\\0\\ indicates an estimate using the total weights, \\r\\ indicates a
replicate from a total of \\R\\, \\p\\ indicates a plausible value from
a total of \\P\\, and \\c\\ indicates an element in a composite estimate
from value a total of \\C\\.

## Examples

``` r
# Creation of replicate weights
RW <- repcreate(df = repdata, # the data frame with all the information
                wt = "wt", # the total weights column name
                jkzone = "jkzones", # the jkzones column name
                jkrep = "jkrep", # the jkreps column name
                repwtname = "REPWT", # the desired name for the rep weights
                reps = 50, # the number of replications
                method = "ICILS") # the name of the method aka the study name

# Non-PVs ----

## Mean with total weights
E0 <- stats::weighted.mean(x = repdata$item01, w = repdata$wt, na.rm = TRUE)
E0
#> [1] 3.616723

## Means by replication
ER <- as.vector(apply(RW,2,function(i){
  stats::weighted.mean(x = repdata$item01, w = i, na.rm = TRUE)
}))
ER
#>  [1] 3.616836 3.619077 3.619685 3.616048 3.617228 3.616271 3.615712 3.614548
#>  [9] 3.614198 3.616319 3.619552 3.615179 3.618895 3.615820 3.616441 3.618920
#> [17] 3.614744 3.616959 3.617163 3.616906 3.616051 3.616516 3.615775 3.616377
#> [25] 3.616088 3.616901 3.618240 3.616946 3.616905 3.614426 3.617217 3.615657
#> [33] 3.616894 3.617397 3.616456 3.616732 3.618875 3.615592 3.614090 3.617166
#> [41] 3.616536 3.617829 3.616730 3.617733 3.615143 3.616488 3.615361 3.617289
#> [49] 3.614539 3.617939

## Standard error by hand
repse(er = ER, e0 = E0, method = "ICILS")
#> [1] 0.009507831

## Standard error with repmean()
repmean(x = "item01",wt = "wt",repwt = RW,df = repdata, method = "ICILS")
#>      N    mean      se      sd    sdse     var   varse
#> 1 4483 3.61672 0.00951 0.66582 0.01106 0.44331 0.01474


# PVs ----

## Mean with total weights
E0 <- sapply(1:5,function(i){
  stats::weighted.mean(x = repdata[,paste0("Math",i)], w = repdata$wt,
                       na.rm = TRUE)
})
E0
#> [1] 0.001906542 0.003607742 0.004161708 0.004761400 0.011764149

## Means by replication
ER <- lapply(1:5, function(j){
  as.vector(apply(RW,2,function(i){
    stats::weighted.mean(x = repdata[,paste0("Math",j)], w = i, na.rm = TRUE)
  }))
})
ER
#> [[1]]
#>  [1] -0.0004463519 -0.0013915224 -0.0025850130  0.0002870367  0.0044730765
#>  [6] -0.0001805856  0.0002410582  0.0023785992  0.0034924355  0.0041688703
#> [11]  0.0027315969 -0.0013845991  0.0022268749  0.0049718156 -0.0012978324
#> [16] -0.0030897725  0.0048509156  0.0052540038 -0.0004898690  0.0033412055
#> [21]  0.0023642694  0.0013834984  0.0034001030  0.0014746144  0.0088756296
#> [26]  0.0002331920  0.0013730372  0.0022503022 -0.0003744138  0.0003263507
#> [31]  0.0002312919  0.0019272061  0.0014090183  0.0037021501 -0.0007417683
#> [36]  0.0012311237  0.0010829931 -0.0005750963 -0.0001365065  0.0039112943
#> [41]  0.0045666259  0.0013860257  0.0011836546  0.0063289710  0.0003306555
#> [46]  0.0057918641  0.0016523310  0.0020257912  0.0050599525 -0.0001661918
#> 
#> [[2]]
#>  [1]  9.151146e-04  6.510882e-04  7.228567e-04  1.952493e-03  6.565699e-03
#>  [6]  1.905492e-03  1.985875e-03  3.868356e-03  6.020139e-03  5.705783e-03
#> [11]  4.097872e-03  1.654258e-03  3.304069e-03  5.635072e-03  1.462933e-03
#> [16] -7.737843e-04  6.646233e-03  5.939021e-03  5.806423e-04  4.916165e-03
#> [21]  4.508211e-03  3.176515e-03  4.551143e-03  3.629559e-03  7.882345e-03
#> [26]  2.299348e-03  3.270834e-03  4.182122e-03  2.349030e-03  2.017711e-03
#> [31]  2.125472e-03  3.479761e-03  5.132273e-03  4.468423e-03  1.680023e-05
#> [36]  3.722541e-03  1.735126e-03  1.933815e-03  1.922831e-03  6.087988e-03
#> [41]  5.405040e-03  3.414345e-03  2.934067e-03  7.678973e-03  2.696717e-03
#> [46]  6.272218e-03  2.623887e-03  3.741456e-03  7.799882e-03  6.984930e-04
#> 
#> [[3]]
#>  [1]  0.0026583656  0.0021790740 -0.0006707041  0.0018652054  0.0063107986
#>  [6]  0.0020340825  0.0038384874  0.0034967078  0.0054091793  0.0043268897
#> [11]  0.0047223490  0.0041556744  0.0032789035  0.0075911736  0.0014657388
#> [16] -0.0007034002  0.0072839378  0.0068366028  0.0036396124  0.0062648908
#> [21]  0.0056344104  0.0024393639  0.0068828897  0.0028940098  0.0105812495
#> [26]  0.0021499642  0.0030532361  0.0046045186  0.0035065188  0.0019424457
#> [31]  0.0023749823  0.0039159275  0.0031923458  0.0052322946  0.0005640173
#> [36]  0.0043281569  0.0044479420  0.0024286577  0.0042483975  0.0060771971
#> [41]  0.0069031531  0.0044022787  0.0040208002  0.0065809468  0.0053812817
#> [46]  0.0065192490  0.0045817970  0.0035340821  0.0072481617  0.0009752336
#> 
#> [[4]]
#>  [1]  4.383773e-03  3.559327e-03  5.864619e-04  3.894075e-03  6.504664e-03
#>  [6]  3.224943e-03  5.393037e-03  2.205715e-03  4.873473e-03  7.407686e-03
#> [11]  6.041855e-03  5.396652e-03  3.951987e-03  8.699723e-03  3.577050e-03
#> [16]  4.747380e-04  5.714284e-03  3.859487e-03  2.529166e-03  3.959320e-03
#> [21]  4.899845e-03  5.009764e-03  5.974318e-03  6.217162e-03  9.901946e-03
#> [26]  3.034254e-03  4.747220e-03  5.384879e-03  3.686355e-03  1.340286e-04
#> [31]  4.390538e-03  3.616674e-03  4.270481e-03  4.722080e-03  2.615539e-03
#> [36]  4.282240e-03  3.838436e-03  3.266618e-03  4.082850e-03  6.256556e-03
#> [41]  9.591444e-03  4.163683e-03  3.908256e-03  9.121732e-03  5.388792e-03
#> [46]  7.359566e-03  4.376797e-03  6.944483e-03  8.257020e-03 -1.956457e-05
#> 
#> [[5]]
#>  [1] 0.009725467 0.008698491 0.008330467 0.009508291 0.014594119 0.010902192
#>  [7] 0.011272954 0.010641134 0.013145012 0.014314543 0.012040355 0.010514591
#> [13] 0.010580799 0.014829591 0.010897973 0.006959437 0.014830595 0.013798612
#> [19] 0.010111703 0.012577115 0.012962675 0.011133646 0.013250337 0.011176855
#> [25] 0.016405816 0.010905420 0.011298591 0.012859263 0.010062711 0.009442275
#> [31] 0.010239275 0.011848111 0.011033086 0.012323824 0.008713693 0.011506440
#> [37] 0.010754110 0.010767209 0.010248354 0.015029168 0.013636771 0.010704115
#> [43] 0.009393543 0.015193123 0.011210993 0.014307394 0.012240039 0.012568500
#> [49] 0.017340583 0.007800080
#> 

## Standard error by hand
repse(er = ER, e0 = E0, method = "ICILS")
#> [1] 0.01634417

## Standard error with repmean()
repmean(x = paste0("Math",1:5),wt = "wt",repwt = RW,df = repdata, method = "ICILS",PV = TRUE)
#>      N    mean      se      sd    sdse     var   varse
#> 1 5000 0.00524 0.01634 1.02282 0.01099 1.04619 0.02248

```
