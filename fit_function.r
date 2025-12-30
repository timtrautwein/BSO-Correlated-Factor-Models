################## fit function #############################
#' @title Model Evaluation Function for BSO
#' @description 
#' This function evaluates the fit of correlated factor models created during the Bee Swarm Optimization process.
#' It takes a vector of item assignments (representing a candidate model), fits the corresponding 
#' structural equation model using lavaan, and calculates various fit indices.
#' 
#' The function is specifically designed for evaluating correlated factor models in the BSO algorithm, 
#' transforming various fit indices into a common "nectar" scale for optimization.
#' Adaptions are needed for other optimization problems.
#' 
#' @param item_assignment A numeric vector where each position represents an item and each value 
#'        indicates which factor the item is assigned to (-1 = not included in model)
#' @param dat A data frame containing the dataset with all variables to be analyzed
#' @param item_names A character vector containing the names of all items (length must match item_assignment)
#' @param fit_crit A character vector specifying which fit criteria to extract and evaluate:
#'        - "cfi": Comparative Fit Index
#'        - "cfi.scaled": Scaled/Robust Comparative Fit Index...
#'        - "rmsea": Root Mean Square Error of Approximation
#'        - "min_omega2": Minimum reliability (omega) across all factors
#'        - "min_loading": Minimum standardized factor loading
#'        - "expl_var": Explained variance ratio
#'        - "n_items": Proportion of items included in the model
#'        - "n_items_min": Whether minimum required items per factor is reached
#'        - "max_fac_cor1": Maximum factor correlation
#'        - "residual": Residual-based fit index using top 10% of residuals
#'        - "residual_mean": Average of largest residuals
#' @param logistic_weights A list of parameters for transforming each fit criterion
#'        using a logistic function, where each element is a named vector with:
#'        - d: horizontal shift parameter (Optimization Threshold)
#'        - a: steepness parameter (positive = higher is better, negative = lower is better)
#' @param nu_weights A numeric vector of weights for combining the transformed fit criteria
#' @param nu_min Threshold value below which any single criterion will cause the overall
#'        fit value to be set to zero (prevents optimization focusing only on some criteria)
#' @param verbose Logical; if TRUE, prints detailed information to console
#' @param return_fit Logical; if TRUE, returns the lavaan fit object for debugging
#' @param debug_fit_mode Logical; if TRUE, propagates warning and error messages
#' @param ignore_warnings Logical; if TRUE, warning messages during model fitting are ignored
#'        (execution continues). If FALSE (default), any warning will be treated as a failure and penalized.
#' @param keep_warnings Logical; if TRUE, collected warning messages are attached as an attribute
#'        ('warnings') to the returned numeric vector (useful for diagnostics).
#' @param ... Additional arguments passed to lavaan's CFA function
#'
#' @return A named numeric vector containing:
#'   - fit_overall: The weighted combination of all fit criteria
#'   - nu: Individual transformed fit values for each criterion
#'   - fit_detailed: Detailed fit criteria
#'   - factor reliabilities: Reliability coefficients for each factor
#'
#' @details
#' The function works in several stages:
#' 1. It builds a structural equation model (CFA) based on the item assignment
#' 2. It fits both a null model and the candidate model using lavaan
#' 3. It calculates various fit indices as specified in fit_crit
#' 4. It transforms all indices to a common scale using logistic functions
#' 5. It combines these transformed values using nu_weights to produce an overall fit value
#'
#' If model fitting fails or produces errors, the function returns a vector of NAs
#' for the fit indices and 0 for the overall fit to prevent BSO from crashing.
#'
#' @note
#' The function is highly customized for correlated factor models and would need substantial 
#' modification to be applied to other types of optimization problems.

fit.function = function(item_assignment, # a vector of item assignments as detailed in the section "Model Specification Search in Bifactor Modeling"
                        dat, # a data frame containing the data set
                        item_names, # a vector of item names (length must equal the length of item assignments)
                        fit_crit = c("cfi", "rmsea", "min_omega2", "min_loading"), # a vector of fit criteria to be extracted
                        # list (length equal to fit_crit) with weights for the logistic transformation (see section The Optimization Function)
                        logistic_weights = list(c(d = 0.9, a = 70),
                                                c(d = 0.06, a = -70),
                                                c(d = 0.4, a = 70),
                                                c(d = 0.33, a = 70)),
                        nu_weights = c(1,1,1,1), # a vector of weights for the aggregation function (see section The Optimization Function)
                        nu_min = 10^(-4), # threshold value, if the nu value of any fit criterion is below this threshold, nu
                        # will be returned as zero to prevent BSO from optimizing a subset of the fit criteria while neglecting
                        # one criterion completely
                        max_nest_fac,
                        verbose = FALSE, # return information in console
                        return_fit = FALSE, # return fit object for debugging
                        debug_fit_mode = FALSE, # option to pass through error messages and warnings to help debugging the fit function
                        ignore_warnings = FALSE, # option to ignore warning messages during model fitting (default: FALSE)
                        keep_warnings = FALSE, # attach warning messages as attribute 'warnings' to the output
                        ... # arguments to be passed to lavaan
){
  
  
  
  #### Warning / Error handling ####
  # The optimization process will frequently propose models that do not converge or yield warnings.
  # By default, any warning is treated as a failure and penalized (fit_overall = 0).
  # If ignore_warnings = TRUE, warning messages are suppressed and evaluation continues.
  # Errors are always penalized to prevent BSO from crashing due to problematic model specifications.

  warn_msgs <- character(0)

  # Create a vector of appropriate length with NAs but set fit_overall to 0
  penalty_out <- function() {
    out_names <- c("fit_overall", paste0("nu.",fit_crit), fit_crit, c(paste0("f", 1:max_nest_fac)))
    out <- c(0, rep(NA, length(out_names) - 1))
    names(out) <- out_names
    out
  }

  # Core model evaluation (separated to allow warning handling via withCallingHandlers)
  eval_fit <- function() {
 
  # show item assignments in verbose mode
    if (verbose) print(item_assignment)

    # the highest number among the item assignments represents the number of nested factors
    n_nest_fac = max(item_assignment)

    # Create a null model that explicitly models all variances
    # This ensures the covariance matrix inside lavaan is complete
    # regardless of which items are selected
    null_model <- paste0(item_names,"~~",item_names, collapse = "\n")


    # #### Build correlated factor model definition ####
    # Create lavaan syntax for each factor, restricting loadings to be positiv
    model=c()
    for(fac in 1:n_nest_fac) {
      model <- paste0(model, "\n","f", fac, " =~ ",
                      paste0("lower(0) * ",
                             item_names[item_assignment==fac], collapse = " + "))

    # Note that we restricted all factor loadings to be positive. 
    }
    #### Model fitting ####

    # write current model in a text file in verbose mode
    if (verbose) writeLines(model, "Model_tmp.txt")

    # print current model to console in verbose mode
    if (verbose) print(model)

    # Fit Null Model (only if needed)
    if ("expl_var" %in% fit_crit) {
    fit_null <- cfa(model = null_model, 
                    data = dat, 
                    se = "none", # do not compute standard errors to speed up computation
                    ...) # additional arguments can be passed to lavaan from the main function (see lavOptions)
    }



    # Fit current candidate model
    fit_model <- cfa(model = model, 
                     data = dat, 
                     std.lv = T, # standardize latent variables
                     se = "none", # do not compute standard errors to speed up computation
                     ...) # additional arguments can be passed to lavaan from the main function (see lavOptions)


    #### Fit criterion evaluation ####
    ## In this section, we compute a number of fit criteria
    ## depending on the input to the fit criterion arguments

    # Initialize vector to store fit indices
    fit_inds <- c()

    # All fit indices in fit_crit which are computed in lavaan will be taken from lavaan
    fit_inds <- fitMeasures(fit_model, fit_crit)

    # -------------------------------------------------------------------
    # Check: Ensure that all requested fit criteria are available (custom must be added if another function is included)
    custom_crit <- c("min_omega2", "min_loading", "expl_var", "n_items", "n_items_min", "max_fac_cor1", "residual", "residual_mean")
    standard_crit <- names(fitMeasures(fit_model))
    for(crit in fit_crit){
      if(!(crit %in% standard_crit) && !(crit %in% custom_crit)){
        stop(paste("Error: The fit criterion", crit, "is not available in the fit function."))
      }
    }
    # -------------------------------------------------------------------

    ### Calculate reliability (omega) for each correlated factor

    rel_tmp <-  suppressMessages( # Suppress semTools informational messages during reliability computation.
	withCallingHandlers( #supress unnecessary warnings...
	  semTools::reliability(fit_model)["omega2",],
	  warning = function(w) {
	  if (grepl("reliability\\(\\) function was deprecated", conditionMessage(w))) {
            invokeRestart("muffleWarning")   # only this warning
        }
        # keep the rest
  }
    )
	)
    # Ensure consistent naming of reliability results 
    # If rel_tmp has only one entry, name it 'f1'
    if(length(rel_tmp) == 1) {
      names(rel_tmp) <- "f1"
    } else {
      names(rel_tmp) <- tolower(gsub("F", "f", names(rel_tmp)))
    }

    # Prepare a fixed-length reliability vector to ensure consistent output structure
    rel <- rep(NA, max_nest_fac ) # length of the vector should be the same irrespective of the number of factors
    names(rel) <- c(paste0("f", 1:max_nest_fac))

    # Store reliability values in the results vector
    rel[match(names(rel_tmp),names(rel))] <- rel_tmp


    ### Calculate and store minimum reliability across all factors if requested
    if ("min_omega2" %in% fit_crit){
           fit_inds <- c(fit_inds, min_omega2 = min(rel, na.rm = TRUE)) 
    }


    ### Extract minimal (standardized) factor loading across all factors
    if ("min_loading" %in% fit_crit){
      # Extract all standardized loadings and their estimation status
      lambda_std <- lavInspect(fit_model, what = "std")$lambda
      lambda_free <- (lavInspect(fit_model, what = "free")$lambda != 0)

      # Only consider parameters that were actually estimated
      lambda_std <- lambda_std[lambda_free]
      fit_inds <- c(fit_inds, min_loading = min(abs(lambda_std)))
    }



    ### Calculate explained variance proportion if requested
    if ("expl_var" %in% fit_crit){

      item_cov <-  lavInspect(fit_null, what = "sampstat")$cov

      item_resid_var <- diag(lavInspect(fit_model, what = "theta"))

      expl_var <- sum(diag(item_cov) - item_resid_var)/ sum(diag(item_cov))

      fit_inds <- c(fit_inds, expl_var = expl_var)
    }


    ### Calculate proportion of items included in the model if requested
    if ("n_items" %in% fit_crit) {
      fit_inds <- c(fit_inds, n_items = sum(item_assignment > 0)/length(item_assignment))
    }


    ### Check if minimum items per factor requirement is met if requested
    if ("n_items_min" %in% fit_crit) {
      # Check if the condition is met and set the value accordingly
      nu_n_items_min <- ifelse(sum(item_assignment > 0) >= n_items_min, 1, 0) # Set logistic weigths on NA
      #print(n_items_min)
      # Add the new indicator to the fitness indicators
      fit_inds <- c(fit_inds, n_items_min = nu_n_items_min)
     }


    ### Calculate maximum factor correlation if requested
    if ("max_fac_cor1" %in% fit_crit) {
      result <- tryCatch({
        # Extract factor correlation matrix
        fac_cor <- lavInspect(fit_model, "cor.lv")

        # Handle special cases based on number of factors
        if(ncol(fac_cor) == 1) {
          max_fac_cor <- 0 # No correlation with only one factor
        } else {
          # Extract lower triangle of correlation matrix (excluding diagonal)
          fac_cor <- fac_cor[lower.tri(fac_cor, diag = FALSE)]
          if(length(fac_cor) == 0) {
            max_fac_cor <- NA # No correlations to extract
          } else {
            # Find maximum correlation value
            max_fac_cor <- max(fac_cor, na.rm = TRUE)
          }
        }
        max_fac_cor
      },
      # Handle potential warnings/errors gracefully
      warning = function(w) {NA},
      error = function(e) {NA})

      fit_inds["max_fac_cor1"] <- result
    }


    ### Calculate residual-based fit index if requested
    if ("residual" %in% fit_crit) {
      # Extract standardized residuals
      residual= resid(fit_model,type="standardized")$cov

      # Focus on the top 10% largest residuals
      num_values_top_10_percent <- ceiling(length(residual) * 0.1)
      top_10_percent_residual <- sort(abs(residual), decreasing = TRUE)[1:num_values_top_10_percent]

      # Calculate sum of top residuals
      sum_residuals_top_10_percent <- sum(top_10_percent_residual,na.rm=T)

      # Compare to ideal threshold (2.58 standard deviations)
      dif_std_residual=length(top_10_percent_residual)*2.58-sum_residuals_top_10_percent
      names(dif_std_residual) <- "dif_std_residual"

      # Apply logistic penalty function - higher values are worse
      pun_log <- sapply(sum_residuals_top_10_percent, function(x) {
       f_x = 1- 1 / (1 + exp(.04*(length(top_10_percent_residual)*2.58-x))) 

        return(f_x)
      })

      # Calculate mean penalty value
      pun_log_mean <- mean(pun_log)

      # Check for invalid values and set to 0 if needed
      if(!exists("pun_log_mean") || is.na(pun_log_mean) || is.null(pun_log_mean)) {
        pun_log_mean <- 0
      }


      fit_inds <- c(fit_inds, residual = pun_log_mean)
    }

    ### Calculate mean of largest residuals if requested    
    if ("residual_mean" %in% fit_crit) {
      # Extract standardized residuals
      residual= resid(fit_model,type="standardized")$cov
      # Get top 10% of residuals by magnitude
      num_values_top_10_percent <- ceiling(length(residual) * 0.1)
      top_10_percent_residual <- sort(abs(residual), decreasing = TRUE)[1:num_values_top_10_percent]

      # Also get top 10 residuals (alternative approach)
      #top_10_residual <- sort(abs(residual), decreasing = TRUE)[1:10]

      # Apply logistic penalty function to each residual
      pun_log <- sapply(top_10_percent_residual, function(x) {
        f_x = 1- 1 / (1 + exp(2*(2.58-x))) # Log function


        return(f_x)
      })

      # Calculate mean penalty value
       pun_log_mean <- mean(pun_log)
      fit_inds <- c(fit_inds, residual_mean = pun_log_mean)
    }



    #### Output preparation ####


    ### Order fit indices to match the order in fit_crit for correct weighting

    fit_inds <- fit_inds[fit_crit]


    # Combine all fit information with factor reliabilities and round to 3 decimals
    fit_detailed <- round(c(fit_inds, rel), 3)

    # Print detailed fit information to console if in verbose mode
    if (verbose) print(fit_detailed)    


    # Transform all fit criteria to a common scale (0-1) using logistic functions
    nu <- sapply(fit_crit, function(i_crit) {
      index <- which(fit_crit == i_crit)

      # Special handling for max_fac_cor1 criterion
      if (i_crit == "max_fac_cor1") {
        fac_cor <- lavInspect(fit_model, "cor.lv") 
        if (ncol(fac_cor) == 1) {
          # For single-factor models, use a default value of 0.5
          out <- 0.5
        } else {
          # Extract relevant correlations and get maximum
          fac_cor <- fac_cor[lower.tri(fac_cor, diag = FALSE)]
          max_fac_cor <- if (length(fac_cor) == 0) NA else max(fac_cor, na.rm = TRUE)
          # Apply logistic transformation
          out <- logistic(max_fac_cor, d = logistic_weights[[index]]["d"], a = logistic_weights[[index]]["a"])
        } 
      } else if (i_crit == "residual" || i_crit == "residual_mean" || i_crit == "n_items_min") {
        # These criteria are already on a 0-1 scale, so no transformation needed
        out <- fit_inds[index]
      } else {
        # Apply standard logistic transformation to other criteria
        out <- logistic(fit_inds[index], d = logistic_weights[[index]]["d"], a = logistic_weights[[index]]["a"])
      }

      names(out) <- NULL
      out

    })

    # Calculate overall fitness ("nectar") by weighted sum of transformed indices
    # If any criterion falls below nu_min threshold, set overall fit to 0
    fit_overall <- ifelse(all(nu >= nu_min), yes = nu %*% nu_weights, no = 0)


    # Print overall fit value if in verbose mode
    if (verbose) print(fit_overall)

    # Return either the full fit object (for debugging) or the fit values
    if (return_fit) {
      fit_model 
    } else {
      c(fit_overall = fit_overall, nu = nu, fit_detailed)

    }

    #### Error handling section ####
    # Any warning or error will result in the fit function returning NA values
    # This prevents BSO from crashing due to problematic model specifications
  }

  res <- tryCatch(
    withCallingHandlers(
      eval_fit(),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))

        # In debug mode: re-throw warnings to inspect the underlying issue (unless warnings are ignored)
        if (debug_fit_mode && !ignore_warnings) stop(w)

        # Default behavior: treat warnings as failures (penalize)
        # Optional behavior: ignore warnings and continue evaluation
        if (ignore_warnings) {
          invokeRestart("muffleWarning")
        } else {
          stop(conditionMessage(w)) # escalate warning to error so it is handled by the error handler below
        }
      }
    ),
    error = function(e) {
      # handle error if fit cirterion isn´t included
      if (grepl("Error: The fit criterion", e$message)) {
        stop(e)  # rethrow the error – so that it isn´t catched by the trycatch
      }

      # In debug mode: re-throw error to inspect the underlying issue
      if (debug_fit_mode) stop(e)

      penalty_out()
    }
  )

  # Optionally attach warning messages for diagnostics
  if (keep_warnings) attr(res, "warnings") <- warn_msgs

  res

  
  
} #end fit function
