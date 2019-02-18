Testing new features in the twang package
================

What's new
----------

The current working version of **twang** that implements xgboost, as well as other enhancements, is available in this repository as part of the branch named **fast**. You need to be sure to switch to the **fast** branch prior to installing the **twang** package. The changes to the package improve the speed and scalability of the main package functions.

### Major improvements

1.  Faster implementation of key components of the TWANG package to improve speed when calculating and optimizing the balance statistics.
2.  The integration of **xgboost** for estimating the propensity score model. This provides users with access to a cutting-edge implementation of gradient boosting, while also providing substantial speed improvements.
3.  A new option (*ks.exact*) was added to improve the speed when calculating the p-value for the weighted two-sample Kolmogorov–Smirnov (KS) statistic.
4.  By default, the standardized difference is calculated at every iteration when either es.mean or es.max is specified. Previous versions of **twang** only provided the user with a standardized differences for a grid of iterations.

### Minor updates

1.  Formula must be specified in the ps function.
2.  dots now pass parameters to booster
3.  Added the n.minobsinnode option to the ps function. Previously, it was always fixed at 10.
4.  Removed the direct optimizer, which was a hidden feature.
5.  Missing data in **xgboost** is handled differently than in **gbm**.
6.  Added a default print for the ps class

### Other changes

In the process of developing these improvements, several issues were discovered in the previous implementation of **twang**. These issues have been resolved in the new implementation of the **twang**, as well as any calls to the legacy code. See [Issue \#1](https://github.com/mattcefalu/twang/issues/1) for details.

1.  When there was missing data in a numeric variable, the denominator used for the standardized difference on the missing data indicator confused ATT and ATE.
2.  Irrespective of the user specified estimand, the standardized difference for factor variables always used the SD among treated as the denominator.
3.  The denominator used in the standardized difference when for the estimand ATE was inconsistent. We changed all instances to use the pooled standard deviation for the full sample based on the survey weights. The previous denominators are described here:
    1.  For numeric variables, the denominator was the standard deviation of the full sample based on the combined propensity score and sampling weights.
    2.  For missing values on numeric variables, ignoring the issue described in (1), the denominator was the standard deviation of the full sample based only on the sampling weights.
    3.  For factor variables, as noted in (2), the denominator was always the SD among treated.

Testing the updated balance calculations
----------------------------------------

The function *dx.wts* has not been updated and still relies on the older version of the balance calculations. Therefore, we will use *dx.wts* to check that the updated balance calculations match the expected output.

First, let's load the data.

``` r
library(twang)
data(lalonde)
```

Next, we will use the *ps* function with ks.exact=TRUE. This is necessary to ensure the same algorithm for estimating the p-value for the KS statistics is used in both calculations.

``` r
M1<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , booster="gbm", ks.exact=T , verbose=F)
w = get.weights(M1 , stop.method = "es.mean")
D1 = dx.wts(x = M1 , estimand="ATT" , vars=c("age","black","hispan","married","nodegree","re74","re78") , treat.var="treat")

# check that the balance summary is the same
(summary(M1)[,-c(8,10)]) - (D1$summary.tab[,-c(1,10)])
```

    ##   n.treat n.ctrl ess.treat ess.ctrl        max.es       mean.es
    ## 1       0      0         0        0  0.000000e+00  1.110223e-16
    ## 2       0      0         0        0 -4.996004e-16 -8.326673e-17
    ## 3       0      0         0        0 -9.714451e-16 -6.245005e-16
    ##          max.ks       mean.ks
    ## 1 -1.110223e-16 -5.551115e-17
    ## 2 -2.775558e-17 -6.938894e-18
    ## 3  1.387779e-17  0.000000e+00

``` r
# check that full balance table is the same
bal.table(D1)$es.mean - bal.table(M1)$es.mean
```

    ##          tx.mn tx.sd ct.mn ct.sd std.eff.sz stat p ks ks.pval
    ## age          0     0     0     0          0    0 0  0       0
    ## black        0     0     0     0          0    0 0  0       0
    ## hispan       0     0     0     0          0    0 0  0       0
    ## married      0     0     0     0          0    0 0  0       0
    ## nodegree     0     0     0     0          0    0 0  0       0
    ## re74         0     0     0     0          0    0 0  0       0
    ## re78         0     0     0     0          0    0 0  0       0

The same strategy can be used for testing when estimand="ATE".

``` r
M2<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATE" , stop.method = c("es.mean","ks.max") , booster="gbm" , ks.exact = T, verbose=F)
w = get.weights(M2 , stop.method = "es.mean")
D2 = dx.wts(x = M2 , estimand="ATE" , vars=c("age","black","hispan","married","nodegree","re74","re78") , treat.var="treat")

# check that the balance summary is the same
(summary(M2)[,-c(8,10)]) - (D2$summary.tab[,-c(1,10)])
```

    ##   n.treat n.ctrl ess.treat ess.ctrl        max.es       mean.es
    ## 1       0      0         0        0  2.220446e-16  1.110223e-16
    ## 2       0      0         0        0  6.661338e-16  2.220446e-16
    ## 3       0      0         0        0 -5.551115e-17 -2.220446e-16
    ##          max.ks       mean.ks
    ## 1 -1.110223e-16 -5.551115e-17
    ## 2 -5.551115e-17  0.000000e+00
    ## 3  5.551115e-17  0.000000e+00

``` r
# check that full balance table is the same
bal.table(D2)$es.mean - bal.table(M2)$es.mean
```

    ##          tx.mn tx.sd ct.mn ct.sd std.eff.sz stat p ks ks.pval
    ## age          0     0     0     0          0    0 0  0       0
    ## black        0     0     0     0          0    0 0  0       0
    ## hispan       0     0     0     0          0    0 0  0       0
    ## married      0     0     0     0          0    0 0  0       0
    ## nodegree     0     0     0     0          0    0 0  0       0
    ## re74         0     0     0     0          0    0 0  0       0
    ## re78         0     0     0     0          0    0 0  0       0

Testing other updates
---------------------

#### New default implementation that uses xgboost

``` r
M1<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max"))
```

#### gbm, but with new balance optimization

``` r
M2<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max"), booster = "gbm")
```

#### xgboost, but with a different tree method

See <https://xgboost.readthedocs.io/en/latest/parameter.html> for the different options.

``` r
M3<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , tree_method = "auto")
```

#### Legacy implementation of ps

``` r
M4<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , version="legacy")
```

#### Keep only every the n-th iteration

The option *n.keep* only keeps every n-th iteration of the propensity score and optimizes over this set instead of all iterations

``` r
M5<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , n.keep = 10)
```

#### Set grid size for an initial search for optimal balance

The *n.grid* option sets the grid size for an initial search of the region most likely to have the minimum. A value of n.grid=50 uses a 50 point grid from 1:n.trees. It finds the minimum, say at grid point 35. It then looks for the actual minimum between grid points 34 and 36.

``` r
M6<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , n.grid = 100)
```

#### Combining an initial grid search with saving only every n-th iteration

The options *n.grid* and *n.keep* can be combined. For now, *n.grid* corresponds to a grid on the kept iterations. For example, n.keep=10 with n.trees=5000 will keep values (10,20,...,5000).n.grid=10 then splits this vector into 10 points. Thus, n.grid\*n.keep must be less than or equal to n.trees.

``` r
M7<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000 , interaction.depth=3 , 
       shrinkage=.01 , estimand = "ATT" , stop.method = c("es.mean","ks.max") , n.grid = 10 , n.keep=10)
```

#### Control xgboost through params

You can use the **gbm** or **xgboost** parameter names directly in the ps function, or you can use the params option for xgboost.

``` r
M8<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , estimand = "ATT" ,  n.trees=5000 ,
       stop.method = c("es.mean","ks.max") , params=list(eta=0.01 , min_child_weight=10 , max_depth=3 ) )
M9<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , estimand = "ATT" ,  n.trees=5000 ,
       stop.method = c("es.mean","ks.max") , eta=0.01 , min_child_weight=10 , max_depth=3  )
```