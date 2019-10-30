#' Describe the effects from a mediation object
#' 
#' Describe the effects, and calculate standard
#' errors and confidence intervals from a mediation object
#' 
#' @param x A mediation object
#' @param y_outcome The outcome; if `NULL`,
#'   then Y must have been provided to the
#'   original mediation function.
#'
#' @method desc_effects mediation
#' @export
desc_effects.mediation <- function(x, y_outcome = NULL) {
  
  # this is just a helper function to calcualte CI and SE for TE and NDE
  get_ci_and_se <- function(eff, w, y_outcome, a_treatment) {
    dsgn <- svydesign(id=~1, weights=~w, data=data.frame(Y=y_outcome, A=a_treatment))
    stderr <- svyglm(Y ~ A, design = dsgn)
    stderr <- summary(stderr)$coeff["A", "Std. Error"]
    ci <- eff + stderr * qnorm(.975) * c(-1, 1)
    return(c(eff, stderr, ci[1], ci[2]))
  }
  
  # this is just a helper function to calcualte CI and SE for NIE
  get_ci_and_se_nie <- function(eff, w, y_outcome, a_treatment, ind) {
    dsgn <- svydesign(id=~1, weights=~wts, data=subset(data.frame(Y=y_outcome, wts=w),
                                                       subset=(a_treatment==ind)))
    stderr <- svymean(~Y, design=dsgn)
    stderr <- SE(stderr)
    ci <- eff + c(stderr) * qnorm(.975) * c(-1, 1)
    return(c(eff, stderr, ci[1], ci[2]))
  }
  
  # Get the data from the original object
  data <- x$data
  
  # Check to make sure that either `y_outcome` has been provided, or that
  # `y_outcome` was provided in the original weighted mediation` calculation
  if (is.null(y_outcome) & !('Y' %in% names(data))) {
    stop(paste("Either a `y_outcome` must be provided, or the original",
               "mediation object must include an outcome in the data set", sep=" "))
  }
  
  # If no `y_outcome` was provided, use the 'Y' from the original data set
  if (is.null(y_outcome)) {
    y_outcome <- data[, 'Y']
  }
  
  # Also, grab the treatment, 'A', from the original data set
  a_treatment <- x$data[, 'A']
  
  # Get the names of the stopping methods
  stop_methods <- c(x$stopping_methods)
  
  # Get the weights (which were saved as attributes)
  w_11 <- attr(x, 'w_11')
  w_00 <- attr(x, 'w_00')
  w_10 <- attr(x, 'w_10')
  w_01 <- attr(x, 'w_01') 
  
  # Now loop through each of the stopping methods to calculate the
  # confidence intervals and standard errors for each effect
  results <- list()
  for (i in 1:length(stop_methods)) {
    
    # Grab the name of the stopping method
    stop_method <- stop_methods[i]
    effects_name = paste(stop_method, "effects", sep = "_")
    
    # Let's get the effects for this stopping method; if there are no
    # effects in the actual mediation object, then we will calculate them
    if (!(effects_name %in% x)) {
      effects <- twang:::calculate_effects(w_11[,i], w_00[,i], w_10[,i], w_01[,i], y_outcome)
    } else {
      effects <- x[[effects_name]]
    }
    
    # We want to calculate the confidence intervals and standard errors,
    # so we collect the weights that we need to calculate these things
    w_te   <- ifelse(!is.na(w_11[, i]), w_11[, i], w_00[, i])
    w_nde0 <- ifelse(!is.na(w_10[, i]), w_10[, i], w_00[, i])
    w_nie1 <- w_11[, i] - w_10[, i]
    w_nde1 <- ifelse(!is.na(w_11[, i]), w_11[, i], w_01[, i])
    w_nie0 <- w_01[, i] - w_00[, i]
    
    # First, we calculate the TE standard error and confidence intervals
    te_res <- get_ci_and_se(effects$total_effect, w_te, y_outcome, a_treatment)
    
    # Second, we calculate the NDE standard errors and confidence intervals
    nde0_res <- get_ci_and_se(effects$natural_direct_effect0, w_nde0, y_outcome, a_treatment)
    nde1_res <- get_ci_and_se(effects$natural_direct_effect1, w_nde1, y_outcome, a_treatment)
    
    # Third, we calculate the NIE standard errors and confidence intervals
    nie0_res <- get_ci_and_se_nie(effects$natural_indirect_effect0, w_nie0, y_outcome, a_treatment, 0)
    nie1_res <- get_ci_and_se_nie(effects$natural_indirect_effect1, w_nie1, y_outcome, a_treatment, 1)
    
    # Package all of these into a data frame
    df_effects <- data.frame(te_res,
                             nde0_res,
                             nie1_res,
                             nde1_res,
                             nie0_res)
    
    # Clean up the structure of the data frame, and append it the results
    df_effects <- transpose(df_effects)
    rownames(df_effects) <- c('TE', 'NDE_0', 'NIE_1', 'NDE_1', 'NIE_0')
    colnames(df_effects) <- c('effect', 'std.err', 'ci.min', 'ci.max')
    results[[effects_name]] <- df_effects
    
  }
  return(results)
}
