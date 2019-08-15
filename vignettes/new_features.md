New features in twang 2.0
================

What's new
----------

The **twang** package has undergone extensive revisions to improve computational efficiency. Off the shelf, **twang** should behave as it has in the past, i.e. code written for prior versions of **twang** should still run without modification. However, behind the scenes, there has been substantial revisions to improve the computational efficiency of the package. In certain situations, results from the updated version of **twang** will no longer replicate prior implementations.

The most up-to-date version of **twang** is available at <https://github.com/mattcefalu/twang/tree/fast> as part of the branch named **fast**.

### Major improvements

1.  Faster implementation of key components of the **twang** package to improve speed when calculating and optimizing the balance statistics. The updated implementation is now default. The user can request that the prior version of **twang** be used with the `version="legacy"` option.
2.  The integration of **xgboost** for estimating the propensity score model. This provides users with access to a cutting-edge implementation of gradient boosting, while also providing substantial speed improvements in larger datasets. The user can request that **xgboost** be used for gradient boosting with the `version="xgboost"` option.
3.  A new option (`ks.exact`) was added to improve the speed when calculating the p-value for the weighted two-sample Kolmogorov–Smirnov (KS) statistic. The default behavior when calculating the weighted two-sample KS p-value differs from previous implementations of **twang**. More information is provided below.
4.  Two new options were added that control the iterations over which the algorithm optimizes balance (`n.keep` and `n.grid`). More information is provided below.
5.  The `n.minobsinnode` option is passed to the underlying gradient boosting algorithm. Previously, this parameter was always fixed at 10.
6.  By default, the specifeied stopping methods are directly optimized. Previous versions of **twang** relied on a optimization routine, which was not guaranteed to find the true minimum.

### Minor updates and notes

1.  Removed the direct optimizer, which was a hidden feature of **twang**. This option is now redudant with the default implementation.
2.  Missing data in **xgboost** is handled differently than in **gbm**. In **gbm**, missing values are placed in their own node, while in **xgboost** missing values are placed in the left or right node based on minimizing the objective function.
3.  Added default print methods for the `ps`, `mnps`, and `iptw` classes that print the basic summaries.

### Other changes

In the process of developing these improvements, several issues were discovered in the previous implementation of **twang**. These issues have been resolved in the new implementation of the **twang**, as well as any calls to the legacy code. See [Issue \#1](https://github.com/mattcefalu/twang/issues/1) for details.

1.  When there was missing data in a numeric variable, the denominator used for the standardized difference on the missing data indicator confused ATT and ATE.
2.  Irrespective of the user specified estimand, the standardized difference for factor variables always used the standard deviation among treated as the denominator.
3.  The denominator used in the standardized difference when for the estimand ATE was inconsistent. We changed all instances to use the pooled standard deviation for the full sample based on the survey weights. The previous denominators are described here:
    1.  For numeric variables, the denominator was the standard deviation of the full sample based on the combined propensity score and sampling weights.
    2.  For missing values on numeric variables, ignoring the issue described in (1), the denominator was the standard deviation of the full sample based only on the sampling weights.
    3.  For factor variables, as noted in (2), the denominator was always the SD among treated.

Testing the updated balance calculations
----------------------------------------

The function `dx.wts` has not been updated and still relies on the older version of the balance calculations. Therefore, we will use `dx.wts` to check that the updated balance calculations match the expected output.

First, let's load some data.

``` r
library(twang)
data(lalonde)
```

Next, we will fit a propensity score model using the `ps` function with `ks.exact=TRUE`. Setting `ks.exact=TRUE` is necessary to ensure the same algorithm for estimating the p-value of the KS statistics is used in all calculations. More details on the `ks.exact` option are provided below.

``` r
M1<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , estimand = "ATT" ,
       n.trees=5000 , stop.method = c("es.mean","ks.max") , ks.exact=T , verbose=F)
w = get.weights(M1 , stop.method = "es.mean")
D1 = dx.wts(x = M1 , estimand="ATT" , vars=c("age","black","hispan","married","nodegree","re74","re78") , treat.var="treat")

# check that the balance summary is the same
all( abs((summary(M1)[,-c(8,10)]) - (D1$summary.tab[,-c(1,10)])) < 10^{-14} )
```

    ## [1] TRUE

``` r
# check that full balance table is the same
all( abs(bal.table(D1)$es.mean - bal.table(M1)$es.mean) < 10^{-14} )
```

    ## [1] TRUE

The same strategy can be used for testing when `estimand="ATE"`.

``` r
M2<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , estimand = "ATE" ,
       n.trees=5000 , stop.method = c("es.mean","ks.max") , ks.exact=T , verbose=F)
w = get.weights(M2 , stop.method = "es.mean")
D2 = dx.wts(x = M2 , estimand="ATE" , vars=c("age","black","hispan","married","nodegree","re74","re78") , treat.var="treat")

# check that the balance summary is the same
all( abs((summary(M2)[,-c(8,10)]) - (D2$summary.tab[,-c(1,10)])) < 10^{-14} )
```

    ## [1] TRUE

``` r
# check that full balance table is the same
all( abs(bal.table(D2)$es.mean - bal.table(M2)$es.mean) < 10^{-14} )
```

    ## [1] TRUE

Description and usage of new options
------------------------------------

All options for the `ps`, `mnps`, and `iptw` are fully described in their corresponding help files. This document provides a brief overview of the new options' behaviors and implementation.

#### Specifying the version of **twang**

Version 2.0 of **twang** defaults to using updated balance and optimization routines. Once installed, the user does not need to modify any existing code to use the update functions that improves computational efficiency. However, if the user wishes to use the prior implementation of **twang**, they can request this by specifying the `version="legacy"` option.

``` r
M1<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , version="legacy" , verbose)
```

The new options of **twang** are not allowed when the user requests the legacy version of the code. If specified together, the user should receive an error message.

``` r
M1<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , version="legacy" , ks.exact=T)
```

    ## Error in ps(treat ~ age + black + hispan + married + nodegree + re74 + : Option ks.exact is not allowed with version='legacy'

#### Specifying the algorithm used to calculated the KS p-value

Version 2.0 of **twang** updated the calculation of the weighted KS p-value to follow the behavior of `ks.test`. Specifically, if `ks.exact=NULL` and the product of the effective sample sizes is less than 10,000, then an approximation based on the exact distribution of the unweighted KS statistic is used. This approximation via the exact distribution can also be requested directly `ks.exact=TRUE`. Otherwise, an approximation based on the asymptotic distribution of the unweighted KS statistic is used.

#### Controlling the optimization

Two new options were added that control how **twang** determines the iterations over which the algorithm optimizes balance. Before describing the options, we first describe how **twang** determines the optimal iteration. The algorithm first creates a grid of points, which is a 25 point grid by default, and assesses the balance at each iteration of the grid. The algorithm finds the optimal iteration based on this grid, and then does a finer search for the optimum between this grid point and its neighbors. For example, say we find that grid point 15 achieves the optimal balance, then the algorithm focuses it search for the optimal iteration between grid points 14 and 16 .

The option `n.keep` is a numeric variable indicating that the algorithm should only consider every `n.keep`-th iteration of the propensity score model and optimize balance over this set instead of all iterations. The default value is 1, which is to optimize over all iterations. Here, we specify `n.keep=50` and see that the optimal iteration is a multiple of 50.

``` r
ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.keep=50 , verbose=F)
```

    ##             n.treat n.ctrl ess.treat ess.ctrl    max.es   mean.es
    ## unw             185    429 185.00000 429.0000 1.3085999 0.4721851
    ## ks.mean.ATE     185    429  74.04940 226.5065 0.5921766 0.2198502
    ## es.mean.ATE     185    429  74.46888 222.1981 0.5996123 0.2194846
    ##                max.ks max.ks.p   mean.ks iter
    ## unw         0.6404460       NA 0.2659359   NA
    ## ks.mean.ATE 0.2898190       NA 0.1235886 1050
    ## es.mean.ATE 0.2934582       NA 0.1238387 1100

The option `n.grid` is numeric variable that sets the grid size for an initial search of the region most likely to minimize the stop.method. When `n.keep=1`, a value of `n.grid=50` uses a 50 point grid from `1:n.trees`.

``` r
M <- ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.grid=50 , verbose=F)
# how many points are in the initial grid?
length(M$iters)
```

    ## [1] 50

The options `n.grid` and `n.keep` can be combined. The option `n.grid` corresponds to a grid of points on the kept iterations as defined by `n.keep`. For example, `n.keep=10` with `n.trees=5000` will keep the iterations `(10,20,...,5000)`. The option `n.grid=10` then splits this vector into 10 points. Thus, `n.grid*n.keep` must be less than or equal to `n.trees`.

``` r
M <- ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000, n.keep=10, n.grid=10 , verbose=F)
summary(M)
```

    ##             n.treat n.ctrl ess.treat ess.ctrl    max.es   mean.es   max.ks
    ## unw             185    429 185.00000 429.0000 1.3085999 0.4721851 0.640446
    ## ks.mean.ATE     185    429  73.24226 224.5819 0.5932635 0.2190305 0.290351
    ## es.mean.ATE     185    429  73.24226 224.5819 0.5932635 0.2190305 0.290351
    ##             max.ks.p   mean.ks iter
    ## unw               NA 0.2659359   NA
    ## ks.mean.ATE       NA 0.1233264 1070
    ## es.mean.ATE       NA 0.1233264 1070

``` r
M <- ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.trees=5000, n.keep=100, n.grid=100 , verbose=F)
```

    ## Error in ps.fast(formula = formula, data = data, n.trees = nrounds, interaction.depth = max_depth, : n.tress must be at least n.grid times n.keep

#### Additional boosting options

This version of **twang** allows the user to specify the **gbm** option `n.minobsinnode`, which controls the minimum number of observations in the terminal nodes of the trees.

``` r
M <- ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , n.minobsinnode = 100 , verbose=F)
M$gbm.obj$n.minobsinnode
```

    ## [1] 100

#### Using xgboost

This version of **twang** supports the use of **xgboost** to perform the gradient boosting for the estimation of the propensity scores. The user request the use of **xgboost** through the `version="xgboost"` option. The behavior of **xgboost** can be controlled using the option names of **gbm**, but can also be controlled directly using the options of **xgboost**. This includes the use of `nrounds`, `max_depth`, `eta`, `subsample`, and `tree_method`. In addition, the list of parameters passed to **xgboost** can be specified with `params`. See <https://xgboost.readthedocs.io/en/latest/parameter.html> for a description of the different options.

``` r
## the following model specifications are equivalent
M3<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , verbose=F,
       n.trees=5000 , shrinkage = 0.05 , interaction.depth = 3, n.minobsinnode=15,  version="xgboost")
M4<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , verbose=F,
       nrounds=5000 , eta = 0.05 , max_depth = 3, min_child_weight=15, version="xgboost")
M5<-ps(treat~age+black+hispan+married+nodegree+re74+re78 , data=lalonde , verbose=F, 
       nrounds=5000 , params=list(eta = 0.05 , max_depth = 3 , min_child_weight=15) , version="xgboost")

# verify this is true by checking balance. 
all.equal(bal.table(M3),bal.table(M4))
```

    ## [1] TRUE

``` r
all.equal(bal.table(M4),bal.table(M5))
```

    ## [1] TRUE
