#' Weighted mediation analysis.
#'
#' Estimate causal mediation mechanism of a treatment
#' using propensity score weighting.
#'
#' For users comfortable with [ps], any options prefaced with
#' `ps_` are passed directly to the `ps()` function.
#'
#' @param formula.med 
#'   A object of class [formula] relating the mediatior(s)
#'   to the covariates (potential confounding variables).
#' @param data 
#'   A dataset of class [data.frame] that includes the treatment indicator, mediator(s), and covariates. 
#' @param a_treatment 
#'   The (character) name of the treatment variable, which must be
#'   dichotomous (0, 1).
#' @param y_outcome 
#'   The (character) name of the outcome variable, y. If this is not provided, then
#'   no effects will be calculated and a warning will be raised. Default : `NULL`.
#' @param total_effect_weights 
#'   A vector of total effect weights, which if left `NULL`
#'   then total_effect_ps must be supplied. Default : `NULL`.
#' @param total_effect_ps 
#'   A ps object that contains the total effect weights,
#;   which if left `NULL` then total_effect_weights must be supplied. Default : `NULL`.
#' @param total_effect_stop_rule
#'   The stopping rule (`ks.mean`, `ks.max`, `es.mean`, `es.max`) for the total effect weights, which 
#'   only needs to be specified if total_effect_ps is provided. Default : `NULL`.
#' @param method
#'   The method for getting weights ("ps", "logistic", or "crossval"). Default : `"ps"`.
#' @param sampw 
#'   Optional sampling weights Default : `NULL`.
#' @param ps_n.trees 
#'   Number of gbm iterations passed on to [gbm]. Default: 10000.
#' @param ps_interaction.depth 
#'   A positive integer denoting the tree depth used in
#'   gradient boosting. Default: 3.
#' @param ps_shrinkage
#'   A numeric value between 0 and 1 denoting the learning rate.
#'   See [gbm] for more details. Default: 0.01.
#' @param ps_bag.fraction 
#'   A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See [gbm] for
#'   more details. Default: 1.0.
#' @param ps_n.minobsinnode An integer specifying the minimum number of observations 
#'   in the terminal nodes of the trees used in the gradient boosting.  See [gbm] for
#'   more details. Default: 10.
#' @param ps_perm.test.iters 
#'   A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If `perm.test.iters=0`
#'   then the function returns an analytic approximation to the p-value. Setting
#'   `perm.test.iters=200` will yield precision to within 3% if the true
#'   p-value is 0.05. Use `perm.test.iters=500` to be within 2%. Default: 0.
#' @param ps_verbose 
#'   If `TRUE`, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: `FALSE`.
#' @param ps_stop.method 
#'   A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are `ks.mean`, `ks.max`, `es.mean`, and `es.max`. `ks` refers to the
#'   Kolmogorov-Smirnov statistic and es refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (`.max`) or the mean (`.mean`). 
#'   Default: `c("ks.mean", "ks.max")`.
#' @param ps_version 
#'  "gbm", "xgboost", or "legacy", indicating which version of the twang package to use.
#'   * `"gbm"` uses gradient boosting from the [gbm] package.
#'   * `"xgboost"` uses gradient boosting from the [xgboost] package.
#'   * `"legacy"` uses the prior implementation of the `ps` function.
#' @param ps_ks.exact `NULL` or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   `NULL`, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting `ks.exact = TRUE` will add substantial
#'   computation time for larger sample sizes. Default: `NULL`.
#' @param ps_n.keep  
#'   A numeric variable indicating the algorithm should only
#'   consider every `n.keep`-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default : 1.
#' @param ps_n.grid 
#'   A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the `stop.method`. A
#'   value of `n.grid=50` uses a 50 point grid from `1:n.trees`. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36.If specified with `n.keep>1`, `n.grid` 
#'   corresponds to a grid of points on the kept iterations as defined by ```n.keep```. Default: 25.
#' @param ps_cv.folds 
#'   A numeric variable that sets the number of cross-validation folds if 
#'   using method='crossval'. Default: 10. 
#' @param ps_keep.data 
#'   A logical variable that determines if the dataset should be saved 
#'   in the resulting `ps` model objects. Default: `FALSE`. 
#' @return mediation object
#'   The `mediation` object includes the following:
#'   - `model_a_res` The model A `ps()` results.
#'   - `model_m1_res` The model M1 `ps()` results.
#'   - `model_m0_res` The model M0 `ps()` results.
#'   - `data` The data set used to compute models
#'   - `stopping_methods` The stopping methods passed to `stop.method`.
#'   - `datestamp` The date when the analysis was run.
#'   - For each `stop.method`, a list with the following:
#'     * `TE` The total effect.
#'     * `NDE_0` The natural direct effect, holding the mediator constant at 0.
#'     * `NIE_1` The natural indirect effect, holding the exposure constant at 1.
#'     * `NDE_1` The natural direct effect, holding the mediator constant at 1.
#'     * `NIE_0` The natural indirect effect, holding the exposure constant at 0.
#'     * `expected_treatment0_mediator0` E(Y(0, M(0)))
#'     * `expected_treatment1_mediator1` E(Y(1, M(1)))
#'     * `expected_treatment1_mediator0` E(Y(1, M(0)))
#'     * `expected_treatment0_mediator1` E(Y(0, M(1)))
#'
#' @seealso [ps]
#' @keywords models, multivariate
#'
#' @export
wgtmed <- function(formula.med,
                               data,
                               a_treatment,
                               y_outcome = NULL,
                               total_effect_wts = NULL,
                               total_effect_ps = NULL,
                               total_effect_stop_rule = NULL,
                               method="ps",
                               sampw = NULL,
                               ps_n.trees = 10000,
                               ps_interaction.depth = 3,
                               ps_shrinkage = 0.01,
                               ps_bag.fraction = 1.0,
                               ps_n.minobsinnode = 10,
                               ps_perm.test.iters = 0,
                               ps_verbose = FALSE,
                               ps_stop.method = c("ks.mean", "ks.max"),
                               ps_version = "gbm",
                               ps_ks.exact = NULL,
                               ps_n.keep = 1,
                               ps_n.grid = 25,
                               ps_cv.folds=10,
                               ps_keep.data=FALSE) {

   # check that data is not a tibble or data.table
   if( class(data)[1] %in% c("tbl_df","tbl","data.table")){
      stop("The wgtmed function currently does not support data.table or tibble. Please convert your data object to a data.frame.")
   }else{
      if (class(data)[1] != "data.frame"){
         warning("Data classes other than data.frame may cause errors." , call.=FALSE)
      }
   }

  # Check the specification of total effect weights 
  # Set total_effect_covars to NULL and set to value later if possible
  total_effect_covars <- NULL

  if(is.null(total_effect_wts) & is.null(total_effect_ps)) stop("Either total effects weights or a total effect ps object must be provided")
  if(!is.null(total_effect_wts)) { 
       if(!is.null(total_effect_ps)) {
		warning("Both total weights and total effects ps object provided. Weights are used.") 
		total_effect_ps	<- NULL
       }
       if(!is.vector(total_effect_wts)) stop("total_effect_wts must be a vector")
       if(!is.null(total_effect_wts) & length(unique(total_effect_wts))>1 ) {warning("Reminder to check that all confounders used for treatment (to obtain supplied\ntotal effect weights) were included in confounders for the mediation model")} 
       }	else{
    			if(total_effect_ps$estimand != "ATE") stop("Total effect must be ATE. Estimand in total_effect_ps != 'ATE'")
    			if(length(total_effect_stop_rule) > 1) warning("Multiple stopping rules provided for total_effect_ps. Only the first is used")
    			if(length(total_effect_stop_rule) == 0) stop("A stopping rule must be specified for the total effects. total_effect_stop_rule must be specified.") 
    			total_effect_wts <- get.weights(total_effect_ps, total_effect_stop_rule[1])
    			var.names.tx <- total_effect_ps$gbm.obj$var.names
    }

  # Check mediator and covariates are in the data
  form.vars	<- trimws(unlist(strsplit(Reduce(paste,deparse(formula.med)),"[~+]")))
  if(!all(form.vars %in% names(data))) stop("All variables in mediation model (formula.med) are not in the dataset")

  # Get the mediators and covariates 
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula.med", "data"), names(mf), 0) 
  mf.med <- mf[c(1, m)]  
  mf.med[[1]] <- as.name("model.frame")
  names(mf.med)[2] <- "formula"
  mf.med$na.action <- na.pass
  mf.med$subset <- rep(FALSE, nrow(data))
  mf.med <- eval(mf.med, parent.frame())
  Terms.med <- attr(mf.med, "terms")
  var.names.med <- attr(Terms.med, "term.labels")
  m_mediators <- trimws(unlist(strsplit(Reduce(paste,deparse(formula.med[[2]])),split="\\+"))) 

  # Check for errors in specification of mediators and treatment and covariates
  #* Treatment must be in the data
  if(!(a_treatment %in% names(data))) stop("Treatment variable is not in the dataset")
  if(!all(data[,a_treatment] %in% c(0:1))) stop("Treatment must be a dichotomous 0,1 variable")

  #* Treatment cannot be a mediator
  if(a_treatment %in% m_mediators){stop("Treatment variable is listed as mediator")}

  #* Treatment should not be specified as a covariate in the mediator model
  if(a_treatment %in% var.names.med){warning("Treatment should not be specified as a predictor in mediator model\nIt has been excluded from formula.med so the function could proceed")
                                    var.names.med <- setdiff(var.names.med, a_treatment)}

  #* The mediator covariates must include any total effect treatment model covariates 
  #* and the variables used to get total effect weights must be in the data
  if(!is.null(total_effect_ps)) {
      if(!identical(var.names.tx, "1") & (length(setdiff(var.names.tx, var.names.med)) > 0)){
        warning("Confounders for treatment must also be included as confounders for mediator \n Omitted confounders added to mediation model")
        var.names.med <- union(var.names.tx, var.names.med) }
      if(!all(var.names.tx %in% names(data))) stop("The variables in the provided total_effect_ps object are not in the dataset")
    } 
  
    #* The outcome cannot equal treatment, the mediator, or any of the covariates and it must be in the data.frame
    if(!is.null(y_outcome)){
       if(y_outcome == a_treatment){stop("The outcome variables equals the treatment variable")}
       if(y_outcome %in% m_mediators){stop("The outcome variables equals a mediator variable")}
       if(y_outcome %in% var.names.med){stop("The outcome variables equals a covariate")}
       if(!(y_outcome %in% names(data))){stop("The outcome variable is not in the dataset")}
	}

  # Create formula for tx without mediator
  #* Model to weights for total effect
  form <- as.formula(paste(a_treatment, "~", paste(var.names.med, collapse="+"))) 
  
  twang:::check_missing(data[,a_treatment])
  twang:::check_missing(data[,m_mediators])
  if (is.null(y_outcome)) {
      warning(paste("The `y_outcome` parameter is NULL. Therefore, only", 
          "weights will be returned; no effects will be calculated.\n", 
          sep = " "))
  }    
  else {
      twang:::check_missing(data[,y_outcome])
  }

  # Generates weights for estimating four population mean
  # E[Y(1.M(1)] -- w_11 Standard IPTW for tx group
  # E[Y(0.M(0)] -- w_00 Standard IPTW for ctrl group
  # E[Y(1.M(0)] -- w_10 Counterfactual mean. p(A = 0 | M, X)/p(A = 1 | M, X) * 1/(1-p(A=1 | X))
  #                     wts only for the tx group
  # E[Y(0.M(0)] -- w_00 Counterfactual mean. p(A = 1 | M, X)/p(A = 0 | M, X) * 1/p(A=1 | X)
  #                     wts only for the ctrl group

  #* Pull w_11 and w_00 from the total
    total_effect_wts <- as.matrix(total_effect_wts)
    w_11 <- total_effect_wts
    w_11[data[,a_treatment] == 0,] <- NA
    w_00 <- total_effect_wts
    w_00[data[,a_treatment] == 1,] <- NA

  if(method=="ps") {
    ps_args <- list(formula = form, n.trees = ps_n.trees, interaction.depth = ps_interaction.depth, 
      shrinkage = ps_shrinkage, bag.fraction = ps_bag.fraction, 
      perm.test.iters = ps_perm.test.iters, verbose = ps_verbose, 
      stop.method = ps_stop.method, version = ps_version, sampw = sampw, 
      ks.exact = ps_ks.exact,keep.data=ps_keep.data)
    if (ps_version != "legacy") {
      append(ps_args, list(n.minobsinnode = ps_n.minobsinnode, 
          n.keep = ps_n.keep, n.grid = ps_n.grid))
    }
    #* Get p(A|X) using same covariates and stopping rules as for mediation model 
    model_a_res <- do.call(ps, c(list(data = data, estimand = "ATE"), 
          ps_args))

    #* Calculate w_10 weights
    #* Note p(A = 0 | M, X)/p(A = 1 | M, X) are the "ATT" weights for the tx group when estimate the average treatment on the control
    #* So we run ps with "tx" = 1- tx and ATT and the estimand
    ps_args$formula <- as.formula(paste("a_treatment0~", paste(c(m_mediators, var.names.med), collapse="+"))) 
    data[,"a_treatment0"] <- 1 - data[,a_treatment]
    model_m0_res <- do.call(ps, c(list(data = data, estimand = "ATT"), 
      ps_args))
    
    #~ Remove a_treatment0 variable in data
    data <- data[,-which(colnames(data)=="a_treatment0")]

    #~  1/(1-p(A=1 | X))
    w_1 <- 1/(1 - model_a_res$ps)
    w_1[data[,a_treatment] == 0, ] <- NA
    
    #~ p(A = 0 | M, X)/p(A = 1 | M, X)
    w_2 <- model_m0_res$w
    w_2[data[,a_treatment] == 0, ] <- NA
    w_10 <- as.matrix(w_2 * w_1)

    #* Calculate w_01 weights
    #* Note p(A = 1 | M, X)/p(A = 0 | M, X) are the "ATT" weights for the crtl group when estimate = ATT
    #* So we run ps with ATT and the estimand
    ps_args$formula <- as.formula(paste(a_treatment,"~", paste(c(m_mediators, var.names.med), collapse="+"))) 
    model_m1_res <- do.call(ps, c(list(data = data, estimand = "ATT"), 
      ps_args))

    #~  1/p(A=1 | X)
    w_1 <- 1/model_a_res$ps
    w_1[data[,a_treatment] == 1, ] <- NA

    #~ p(A = 1 | M, X)/p(A = 0 | M, X)
    w_2 <- model_m1_res$w
    w_2[data[,a_treatment] == 1, ] <- NA
    w_01 <- as.matrix(w_2 * w_1)

    #~  Make w_11 and w_00 matrices same dimensions as w_10 and w_01
    if(ncol(w_10)>ncol(w_11)) {
      w_11 <- matrix(as.vector(w_11),ncol=ncol(w_10),nrow=nrow(w_11))
      w_00 <- matrix(as.vector(w_00),ncol=ncol(w_10),nrow=nrow(w_00))
    }
    if(!is.null(total_effect_ps)) {
      colnames(w_11) <- colnames(w_00) <- rep(paste(total_effect_stop_rule[1],total_effect_ps$estimand,sep="."),ncol(w_11))
    } 
  }
  
  if(method!="ps") {
    form_m <- as.formula(paste(a_treatment,"~", paste(c(m_mediators, var.names.med), collapse="+"))) 
    if(method=="crossval") {
      #* Fit total effects model to get p(A|X)  
      model_a_res <- gbm(formula=form, data = data, weights = sampw, 
          distribution = "bernoulli", n.trees = ps_n.trees, 
          interaction.depth = ps_interaction.depth, n.minobsinnode = ps_n.minobsinnode, 
          shrinkage = ps_shrinkage, bag.fraction = ps_bag.fraction, train.fraction = 1, 
          verbose = ps_verbose, keep.data = FALSE, cv.folds=ps_cv.folds)
      best.iter <- gbm.perf(model_a_res, method="cv",plot.it=FALSE)
      model_a_preds <- predict(model_a_res, n.trees=best.iter, newdata=data, type="response")

    #* Fit mediation model
      model_m0_res <- gbm(formula=form_m, data = data, weights = sampw, 
          distribution = "bernoulli", n.trees = ps_n.trees, 
          interaction.depth = ps_interaction.depth, n.minobsinnode = ps_n.minobsinnode, 
          shrinkage = ps_shrinkage, bag.fraction = ps_bag.fraction, train.fraction = 1, 
          verbose = ps_verbose, keep.data = FALSE, cv.folds=ps_cv.folds)
      best.iter <- gbm.perf(model_m0_res, method="cv",plot.it=FALSE)
      model_m0_preds <- predict(model_m0_res, n.trees=best.iter, newdata=data, type="link")
     }
    if(method=="logistic") {
      #* Fit total effects model to get p(A|X)  
      #* Suppress warnings to suppress warning that arises if specify sampling weights 
      suppressWarnings( model_a_res <- glm(form,data=data,family="binomial",weights=sampw))
      model_a_preds <- predict(model_a_res,type="response")

      #* Fit mediation model 
      suppressWarnings(model_m0_res <- glm(form_m,data=data,family="binomial",weights=sampw))
      model_m0_preds <- predict(model_m0_res,type="link")
    }
    #~  1/(1-p(A=1 | X))
    w_1 <- 1/(1 - model_a_preds)
    w_1[data[,a_treatment] == 0] <- NA  

    #~ p(A = 0 | M, X)/p(A = 1 | M, X)
    w_2 <- 1/exp(model_m0_preds)
    w_2[data[,a_treatment] == 0] <- NA
    w_10 <- as.matrix(w_2 * w_1)

    #* Calculate w_01 weights: p(A = 1 | M, X)/p(A = 0 | M, X) * 1/p(A=1 | X)
    model_m1_res <- model_m0_res

    #~  1/p(A=1 | X)
    w_1 <- 1/model_a_preds
    w_1[data[,a_treatment] == 1] <- NA

    #~ p(A = 1 | M, X)/p(A = 0 | M, X)
    w_2 <- exp(model_m0_preds)
    w_2[data[,a_treatment] == 1] <- NA
    w_01 <- as.matrix(w_2 * w_1)
    
    ps_stop.method <- method

  }
  results <- list(method=method,model_m0 = model_m0_res, model_m1 = model_m1_res, 
      model_a = model_a_res, mediator_names = m_mediators, covariate_names = var.names.med, 
      a_treatment=a_treatment, y_outcome=y_outcome, 
      stopping_methods = ps_stop.method, data = data, datestamp = date())
  class(results) <- "mediation"
  attr(results, "w_11") <- w_11
  attr(results, "w_00") <- w_00
  attr(results, "w_10") <- w_10
  attr(results, "w_01") <- w_01
  if (is.null(y_outcome)) {
      return(results)
  }
  if(method=="ps") {
    stop_methods <- c(ps_stop.method)
    for (i in 1:length(stop_methods)) {
        stop_method <- stop_methods[i]
        effects_name = paste(stop_method, "effects", sep = "_")
        w_11_temp <- w_11[, i]
        w_00_temp <- w_00[, i]
        w_10_temp <- w_10[, i]
        w_01_temp <- w_01[, i]
        results[[effects_name]] <- twang:::calculate_effects(w_11_temp, 
            w_00_temp, w_10_temp, w_01_temp, data[,y_outcome],sampw=sampw)
  } }
  else {
  results[[paste0(method,"_effects")]] <- twang:::calculate_effects(w_11, 
            w_00, w_10, w_01, data[,y_outcome],sampw=sampw)
  }
  return(results)
}
