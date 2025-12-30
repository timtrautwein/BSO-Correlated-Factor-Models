# BSOs Main-Function

# Install packages if required 
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Load all required packages
install_and_load("psych")        # For psychometric functions and factor analysis
install_and_load("lavaan")       # For structural equation modeling 
install_and_load("parallel")     # For parallel computing capabilities
install_and_load("GPArotation")  
install_and_load("OpenMx")       # For extended structural equation modeling (and the Holzinger-Dataset for Testing Purposes)
install_and_load("semTools")     # For additional SEM tools and utilities
install_and_load("ggplot2")      # For visualization
install_and_load("scales")       # For scale transformations and formatting
install_and_load("tidyverse")    # For data manipulation and visualization

# Load necessary helper functions from external files
source("tool_functions.R")  # Contains subfunctions for BSO operations
source("fit_function.R")    # Contains the fitness evaluation function

#' @title Bee Swarm Optimization (BSO) for Correlated Factor Analysis
#' @description 
#' Implements the Bee Swarm Optimization algorithm for discovering optimal factor structures
#' in psychometric data. The algorithm mimics the behavior of bee colonies, where scout bees
#' explore new solutions (major model changes) and onlooker bees refine existing solutions
#' (minor model changes). The BSO approach avoids local optima through its exploratory nature.
#'
#' @param item_names Vector of item names to be used in the factor model
#' @param data Data frame containing the variables
#' @param max_iter Maximum number of iterations without improvement before stopping
#' @param n_bees Total number of bees (Number of Solutions to explore in one Iteration) to use in the optimization
#' @param nest_site_bees Number of nest site scout bees for initial exploration (default = n_bees). 
#' @param use_scout_nest_init Logical. If TRUE, use Scout-Nest-Bees initialization; if FALSE, use original random initialization (default = TRUE)
#' @param min_efa_loading Minimum EFA loading required for items to be retained in CFA (default = 0, i.e., no minimum)
#' @param n_start_items Number of items in start model (for original initialization)
#' @param percent_scouts Proportion of bees assigned as scouts (major model changes)
#' @param top_best Number of best solutions to focus modification efforts on
#' @param min_fac Minimum number of factors (preferred alias for min_nest_fac)
#' @param max_fac Maximum number of factors (preferred alias for max_nest_fac)
#' @param min_nest_fac Minimum number of factors (legacy name; backward compatibility)
#' @param max_nest_fac Maximum number of factors (legacy name; backward compatibility)
#' @param depletion Maximum "age" of a solution before it's considered depleted (exhausted)
#' @param summaryfile File path to save detailed results of all explored models
#' @param summaryfile_fin File path to save only the best final model
#' @param write_solutions Logical; if TRUE, write solutions to files (default is TRUE)
#' @param scouts Vector specifying number of scout bees per top solution (calculated from percent_scouts)
#' @param onlookers Vector specifying number of onlooker bees per top solution (calculated from percent_scouts)
#' @param fit_crit Vector of fit criteria to extract and optimize
#' @param logistic_weights List of parameters for logistic transformation of fit criteria
#' @param nu_weights Weights for the aggregation function that combines fit criteria
#' @param nu_min Threshold value for minimum acceptability of individual fit criteria
#' @param n_items_min Minimum number of items required in a factor
#' @param verbose Logical; if TRUE, print detailed progress information
#' @param debug_fit_mode Logical; if TRUE, return fit object for debugging
#' @param parallel Logical; if TRUE, use parallel processing
#' @param nCores Number of CPU cores to use when parallel=TRUE
#' @param seed Random seed for reproducibility (in case of parallel processing the number of cores must be the same as before for reproducibility)
#' @param fun Logical; if TRUE, display ASCII art (Easter egg feature)
#' @param plot_nectar Logical; if TRUE, create plots of optimization progress
#' @param plot_list List of aesthetic parameters for the nectar plot
#' @param cluster_mode Logical; if TRUE, disable plotting (for use on clusters)
#' @param balance_n_fac Logical; if TRUE, start with uniform factor distribution
#' @param ignore_warnings Logical; if TRUE, warnings during model fitting are muffled and the model is evaluated; if FALSE, any warning is treated as a failed evaluation (fit_overall = 0).
#' @param ... Additional arguments passed to lavaan
#' 
#' @return A ggplot object showing optimization progress (if plot_nectar=TRUE)
#' @details
#' The BSO algorithm works as follows:
#' 1. Generate initial solutions via nest scout bees using explanatory factor models as starting points
#' 2. Evaluate each model using fit criteria
#' 3. Select top solutions for further exploration
#' 4. Apply scout bees to make major changes (add/remove/merge/split factors)
#' 5. Apply onlooker bees to make minor changes (add/remove/swap items)
#' 6. Re-evaluate new models and keep the best solutions
#' 7. Continue as long as a better model is found within max_iter iterations
#'
#' The "correlated factors" terminology refers to factor models where factors are allowed to
#' correlate with each other; nest is used for correlated factors for consistency with the original function.

BSO = function(item_names,      # Vector of item names to be analyzed
               data,            # Data frame containing the data set
               max_iter,        # Maximum iterations without improvement before stopping
               n_bees,          # Total number of bees (computational agents)
               nest_site_bees = n_bees, # Number of nest site scout bees for initial exploration
               n_start_bees = NULL, # Deprecated parameter for backward compatibility
               use_scout_nest_init = TRUE, # Use Scout-Nest-Bees initialization (TRUE) or original random initialization (FALSE)
               min_efa_loading = 0, # Minimum EFA loading for items to be retained (default = 0, no minimum)
               n_start_items = length(item_names),  # Number of items in starting models
               percent_scouts,  # Proportion of scout bees (major model changes)
               top_best,        # Number of best solutions to focus on
               min_fac = NULL , # minimum number of correlated factors allowed (preferred alias; overrides min_nest_fac if provided)
               max_fac = NULL, # maximum number of correlated factors allowed (preferred alias; overrides max_nest_fac if provided)
               min_nest_fac = NULL, # minimum number of factors (legacy name; kept for backward compatibility)
               max_nest_fac = NULL, # maximum number of factors (legacy name; kept for backward compatibility)
               depletion,       # Maximum "age" of a solution before it's considered exhausted and won´t be longer used as top solution
               summaryfile,     # File path for detailed results output
               summaryfile_fin, # File path for final best model output
               write_solutions = TRUE,  # Whether to write solutions to files
               # Calculate scout and onlooker distribution based on solution quality (weighted by rank)
               scouts = round(prop.table(top_best:1) * n_bees * percent_scouts),     
               onlookers = round(prop.table(top_best:1) * n_bees * (1 - percent_scouts)),  
               # Note: Distribution weighted by rank – top solutions receive more attention
			   fit_crit = c("cfi", "rmsea", "min_omega2", "min_loading"),  # Fit criteria to optimize (for other implemented criteria see fit_function)
               # Parameters for logistic transformation of fit indices (shape and scale parameters)
               logistic_weights = list(c(d = 0.9, a = 70),    # CFI transformation (higher is better)
                                       c(d = 0.06, a = -70),  # RMSEA transformation (lower is better)
                                       c(d = 0.4, a = 70),    # min_omega2 transformation (higher is better)
                                       c(d = 0.33, a = 70)),  # min_loading transformation (higher is better)
               nu_weights = c(1,1,1,1),  # Equal weights for combining fit criteria by default
               nu_min = 0,               # Minimum threshold for any single fit criterion
               n_items_min = 3,          # Minimum items per factor
               verbose = FALSE,          # Whether to print detailed information (not compatible with parallel mode)
               debug_fit_mode = FALSE,   # Whether to return full fit object for debugging
               parallel = TRUE,          # Whether to use parallel processing
               nCores = parallel::detectCores() - 2,  # Number of cores to use (leaves 2 for system)
               seed = 1,                 # Random seed for reproducibility
               fun = TRUE,               # Easter egg: display ASCII art
               plot_nectar = TRUE,       # Whether to create convergence plot
               # Aesthetic parameters for the nectar (fitness) plot
               plot_list = list(xlim = c(0,max_iter*2), 
                                ylim = c(0,sum(nu_weights)),
                                ylab = "Overall Nectar Value",
                                xlab = "Iteration",
                                jitter_width = 0.5,
                                alpha = 0.2,
                                size = 1,
                                width = 10,
                                height = 8,
                                dpi = 450),
               cluster_mode = FALSE,     # For use on computing clusters (disables plots)
               balance_n_fac = FALSE,    # Whether to balance factor counts in initial models
               ignore_warnings  = FALSE, # If TRUE, ignore/muffle warnings in fit.function(); if FALSE, treat warnings as evaluation failures (fit_overall =0).
               ...){                     # Additional arguments passed to lavaan
  
  # Handle backward compatibility: if n_start_bees is provided, use it for nest_site_bees
  if (!is.null(n_start_bees)) {
    warning("Parameter 'n_start_bees' is deprecated. Please use 'nest_site_bees' instead.")
    nest_site_bees <- n_start_bees
  }
  # Resolve factor-range arguments with backward compatibility (min/max factors)
	if (!is.null(min_fac)) {
	min_nest_fac <- min_fac
	}
	if (!is.null(max_fac)) {
	max_nest_fac <- max_fac
	}

	# Backward compatibility: warn when legacy argument names are used
	if (!is.null(min_nest_fac) && is.null(min_fac)) {
	warning("Parameter 'min_nest_fac' is deprecated. Please use 'min_fac' instead.")
	}
	if (!is.null(max_nest_fac) && is.null(max_fac)) {
	warning("Parameter 'max_nest_fac' is deprecated. Please use 'max_fac' instead.")
	}
	
	#### Input checks ####
	# fit_crit, logistic_weights, and nu_weights must match in length because each fit criterion
	# is transformed and weighted element-wise. Mismatched lengths would misalign indices, so we stop early with a clear error.
	if (!is.list(logistic_weights)) {
	stop("Error: 'logistic_weights' must be a list with one (d,a) entry per element of 'fit_crit'.")
	}
	if (length(logistic_weights) != length(fit_crit)) {
	stop(paste0("Error: Length mismatch. 'logistic_weights' has length ", length(logistic_weights),
				" but 'fit_crit' has length ", length(fit_crit), ". They must be equal."))
	}
	if (length(nu_weights) != length(fit_crit)) {
	stop(paste0("Error: Length mismatch. 'nu_weights' has length ", length(nu_weights),
				" but 'fit_crit' has length ", length(fit_crit), ". They must be equal."))
	}
  
  # see for yourself
  if(fun){
    cat(readLines("fun1.txt", warn = F), sep = "\n")
  }
  
  # Initialize the nectar (fitness) plot if enabled
  if (plot_nectar){
    # Create base plot with proper axis limits and labels
    conv_plot <- ggplot() + 
      xlim(plot_list$xlim[1], plot_list$xlim[2]) +
      ylim(plot_list$ylim[1], plot_list$ylim[2]) +
      xlab(plot_list$xlab) +
      ylab(plot_list$ylab) 
    # Display empty plot at start (unless in cluster mode)
    if(!cluster_mode) plot(conv_plot)
  }
  
  # Set up parallel processing environment if enabled
  if (parallel) {
    # Create PSOCK cluster for Windows compatibility
    myCL <- makePSOCKcluster(names = nCores, outfile = ifelse(verbose, paste0("parallel_verbose",Sys.time(),".txt"), ""))
    
    # Set up random number generation for reproducibility in parallel mode
    RNGkind("L'Ecuyer-CMRG") # set seed for reproducibility to mirror parallels default
    clusterSetRNGStream(cl = myCL, iseed = seed)
    
    # Export all variables and functions to the cluster nodes
    clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
    clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    
    # Load necessary packages on each cluster node
    clusterEvalQ(cl = myCL, expr = {
      library(psych)
      library(lavaan)
      library(semTools)
      library(ggplot2)
      library(scales)
      library(GPArotation)
    })
  } else {
    set.seed(seed, kind = "L'Ecuyer-CMRG")  # set seed for reproducibility to mirror parallels default
  }
  
  # Store total number of items
  n_items <- length(item_names) #number of items in total
  

  ########################### Phase 1: Generate Initial Solutions ##############################
  # Choose initialization method based on use_scout_nest_init parameter
  
  if (use_scout_nest_init) {
    # Use Scout-Nest-Bees initialization with EFA
    if (verbose) cat("Using Scout-Nest-Bees initialization with EFA...\n")
    
    if (!parallel) {
      solutions <- initialize.scout.nest.bees(
        item_names = item_names,
        data = data,
        nest_site_bees = nest_site_bees,
        min_nest_fac = min_nest_fac,
        max_nest_fac = max_nest_fac,
        n_items_min = n_items_min,
        min_efa_loading = min_efa_loading,
        fit_crit = fit_crit,
        logistic_weights = logistic_weights,
        nu_weights = nu_weights,
        nu_min = nu_min,
        debug_fit_mode = debug_fit_mode,
		ignore_warnings = ignore_warnings,
        verbose = verbose,
        seed = seed,
        ...
      )
    } else {
      # Export additional functions needed for Scout-Nest-Bees
      clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
      clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
      clusterEvalQ(cl = myCL, expr = {
        library(psych)
        library(lavaan)
        library(semTools)
        library(ggplot2)
        library(scales)
        library(GPArotation)
      })
      
      # Scout-Nest-Bees already returns a complete dataframe, so we just run it once
      solutions <- initialize.scout.nest.bees(
        item_names = item_names,
        data = data,
        nest_site_bees = nest_site_bees,
        min_nest_fac = min_nest_fac,
        max_nest_fac = max_nest_fac,
        n_items_min = n_items_min,
        min_efa_loading = min_efa_loading,
        fit_crit = fit_crit,
        logistic_weights = logistic_weights,
        nu_weights = nu_weights,
        nu_min = nu_min,
        debug_fit_mode = debug_fit_mode,
		ignore_warnings = ignore_warnings,
        verbose = verbose,
        seed = seed,
        ...
      )
    }
    
  } else {
    # Use original random initialization
    if (verbose) cat("Using original random initialization...\n")
    
    if (!parallel){
      solutions <- sapply(1:nest_site_bees, function(i_scout){ 
        if (verbose) print(i_scout)
        
        initialize.bees(item_names = item_names,
                        data = data,
                        max_iter = max_iter,
                        n_bees = n_bees,
                        n_start_items = n_start_items,
                        fit_crit = fit_crit,
                        logistic_weights = logistic_weights,
                        nu_weights = nu_weights,
                        nu_min = nu_min,
                        balance_n_fac = balance_n_fac,
                        debug_fit_mode = debug_fit_mode,
					    ignore_warnings = ignore_warnings,
                        verbose = verbose,
                        ...)
      }, simplify = "matrix")
    } else if (parallel){
      solutions <- parSapply(myCL, 1:nest_site_bees, function(i_scout){ 
        if (verbose) print(i_scout)
        initialize.bees(item_names = item_names,
                        data = data,
                        max_iter = max_iter,
                        n_bees = n_bees,
                        n_start_items = n_start_items,
                        fit_crit = fit_crit,
                        logistic_weights = logistic_weights,
                        nu_weights = nu_weights,
                        nu_min = nu_min,
                        balance_n_fac = balance_n_fac,
                        debug_fit_mode = debug_fit_mode,
						ignore_warnings = ignore_warnings,
                        verbose = verbose,
                        ...
        )
      }, simplify = "matrix")
    }
    
    # Transpose solutions object for convenience purposes
    solutions <- data.frame(t(solutions))
  }

  ########################### Dynamic adjustment of top_best ##############################
  
  # Check how many valid solutions we actually have (not all NA)
  valid_solutions <- !is.na(solutions$fit_overall) & solutions$fit_overall >= 0 # mark solutions with a valid evaluation (fit_overall not NA and >= 0)
  n_valid_solutions <- sum(valid_solutions)
  
  if (verbose) {
    cat("\n=== DYNAMIC ADJUSTMENT ===\n")
    cat("Total solutions generated:", nrow(solutions), "\n")
    cat("Valid solutions (non-NA, >0):", n_valid_solutions, "\n")
    cat("Originally requested top_best:", top_best, "\n")
  }
  
  # Store original top_best for reference
  original_top_best <- top_best
  
  # Adjust top_best if we have fewer valid solutions than requested
  if (n_valid_solutions < top_best) {
    top_best <- max(1, n_valid_solutions)  # At least 1, but not more than available
    
    if (verbose) {
      cat("ADJUSTMENT: Reducing top_best to:", top_best, "\n")
    }
    
    warning(paste0("Reduced top_best from ", original_top_best, " to ", top_best, 
                   " because only ", n_valid_solutions, " valid solutions were generated."))
  }
  
  # Filter solutions to only keep valid ones for the initial sorting
  if (n_valid_solutions > 0) {
    solutions_valid <- solutions[valid_solutions, ]
    solutions_invalid <- solutions[!valid_solutions, ]
    
    # Sort valid solutions by fit_overall
    solutions_valid <- solutions_valid[order(solutions_valid$fit_overall, decreasing = TRUE), ]
    
    # Recombine: valid solutions first, then invalid ones
    solutions <- rbind(solutions_valid, solutions_invalid)
  } else {
    stop("No valid solutions generated. All fit_overall values are NA or <= 0. Check your fit function and data.")
  }
  
  # Recalculate scouts and onlookers based on the adjusted top_best
  if (is.null(scouts) || is.null(onlookers)) {
    scouts <- round(prop.table(top_best:1) * n_bees * percent_scouts)
    onlookers <- round(prop.table(top_best:1) * n_bees * (1 - percent_scouts))
    
    if (verbose) {
      cat("Recalculated scouts:", scouts, "\n")
      cat("Recalculated onlookers:", onlookers, "\n")
    }
  }
  
  # Validate bee distribution: total scouts + onlookers must equal n_bees
  if (sum(scouts + onlookers) != n_bees) {
    n_bees <- sum(scouts + onlookers)
    warning(paste0("The current combination of n_bees, percent_scouts and top_best does not work out, I changed n_bees to "), n_bees)
  }
  
  if (verbose) {
    cat("Final adjusted parameters:\n")
    cat("- top_best:", top_best, "\n")
    cat("- n_bees:", n_bees, "\n")
    cat("- scouts:", paste(scouts, collapse = ", "), "\n")
    cat("- onlookers:", paste(onlookers, collapse = ", "), "\n")
    cat("==========================\n\n")
  }
  
  # Add initial solutions to plot if plot_nectar is TRUE
  if (plot_nectar && nrow(solutions) > 0) {
    # Prepare initial solutions for plotting
    init_plot_data <- solutions
    init_plot_data$iter <- 0  # Mark as iteration 0
    
    # Determine role based on initialization method
    if (use_scout_nest_init) {
      init_plot_data$role <- as.factor("nest_site_scouts")
    } else {
      init_plot_data$role <- as.factor("random_init")
    }
    
    init_plot_data$n_nest_fac <- as.factor(init_plot_data$n_nest_fac)
    
    # Add to plot with different shape/color to distinguish from regular iterations
    suppressWarnings(expr = {
      conv_plot <- conv_plot + 
        geom_jitter(data = init_plot_data, 
                    aes(x = iter,
                        y = fit_overall,
                        color = n_nest_fac,
                        shape = role),
                    width = plot_list$jitter_width * 0.5,  # Less jitter for initial points
                    alpha = plot_list$alpha * 1.5,  # Slightly more prominent
                    size = plot_list$size * 1.2) +  # Slightly larger
        scale_shape_manual(values = c("nest_site_scouts" = 17, "random_init" = 18, 
                                     "scout" = 16, "onlooker" = 15))  # Different shapes
      
      if(!cluster_mode) plot(conv_plot)
    })
  }
  
  ################################################ Phase 2: Main Optimization Loop ##############################################
  # In this phase, scout bees perform major modifications and onlooker bees perform minor modifications
  # to the top solutions in each iteration until convergence
  
  
  iter <- 1        # Overall iteration counter
  counter <- 1     # Counter for iterations without improvement (resets when improvement found)
  best_fit <- solutions$fit_overall[1]  # Track best fitness seen so far
  best_solution <- solutions[1,]        # Store best solution seen so far

  # Main optimization loop
  while (counter <= max_iter){ # Continue until max_iter iterations without improvement
    start <- Sys.time() # Measure runtime
    
    # Adaptive top_best: if the initial phase yielded too few valid solutions, gradually increase top_best as more valid solutions accumulate (up to original_top_best).
    # Rationale: exploring more top solutions becomes meaningful once enough distinct, valid solutions exist.
    current_valid_solutions <- sum(!is.na(solutions$fit_overall) & solutions$fit_overall > 0)
    
    # Gradually increase top_best as we get more solutions, but don't exceed original request
    if (current_valid_solutions > top_best && top_best < original_top_best) {
      new_top_best <- min(original_top_best, current_valid_solutions)
      
      if (new_top_best > top_best) {
        if (verbose) {
          cat("ITERATION", iter, ": Increasing top_best from", top_best, "to", new_top_best, "\n")
        }
        
        top_best <- new_top_best
        
        # Recalculate scouts and onlookers for the new top_best
        scouts <- round(prop.table(top_best:1) * n_bees * percent_scouts)
        onlookers <- round(prop.table(top_best:1) * n_bees * (1 - percent_scouts))
        
        # Adjust n_bees if necessary
        if (sum(scouts + onlookers) != n_bees) {
          n_bees <- sum(scouts + onlookers)
          if (verbose) {
            cat("Adjusted n_bees to", n_bees, "\n")
          }
        }
      }
    }
    # Select top solutions to explore in this iteration
    tmp_solutions <- solutions[1:top_best,] 
    
    # Check which solutions are still "fresh" (not depleted)
    fresh <- solutions$age <= depletion 
    
    # Increase age of best/top solutions that are still fresh
    solutions$age[fresh][1:top_best] <- as.numeric(solutions$age[fresh][1:top_best]) + 1 
    
    # Prepare a data frame that defines which bees will modify which solutions
    allBees <- data.frame(i_bee_ind = 1:n_bees)
    # Assign bees to solutions based on scouts and onlookers vectors
    allBees$solution_focus = c(rep(c(1:top_best), times = scouts), rep(c(1:top_best), times = onlookers)) 
    # Assign bee roles - scouts first, then onlookers
    allBees$role = c(rep("scout", sum(scouts)), rep("onlooker", sum(onlookers)))
    # Add the current model parameters for each solution
    allBees <- data.frame(allBees, tmp_solutions[allBees$solution_focus,item_names]) 
   
    # Process each bee's operation (either in serial or parallel)
    if (!parallel){
      # Single-core implementation
      new_solutions <- sapply(1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        # Get the current item assignment for this bee to modify
        item_assignment_orig <- as.numeric(as.character(unlist(allBees[i_bee, item_names])))
        names(item_assignment_orig) <- item_names
      
        # Apply appropriate bee operation based on the bee's role
        if (allBees$role[i_bee] == "scout"){
          # Scout bees make major model changes (add/remove/merge/split factors)
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          # Ensure consistent factor numbering
          item_assignment <-reorder_factors(item_assignment = item_assignment)
          
        } else if (allBees$role[i_bee] == "onlooker") {
          # Onlooker bees make minor model changes (add/remove/swap items)
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <-reorder_factors(item_assignment = item_assignment)
        } else {stop("Fatal error: Bee swarm out of control.")} # if an unknown role occurs...
        
        # Check if this model configuration has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) identical(item_assignment, i_solution))

        if (any(compare_item_assignment)){ # if you saw this item configuration before
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          # the new age is the maximal age of this solution + 1
          # (because solutions is sorted by fit)
          new_sol$age <- max(solutions$age[which(compare_item_assignment)]) + 1 
        } else { 
          # If this is a new model, evaluate it

          # Calculate fit metrics for this model
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
							  ignore_warnings = ignore_warnings,
							  max_nest_fac = max_nest_fac,
                              ...) #evaluate model
          
          
          # Save model results: item-factor allocation, seed, iteration, 
          # number of factors, age, and fit statistics
          n_nest_fac <- max(item_assignment)
          
          # include metadata for consistency
          new_sol <- c(unlist(item_assignment), 
                      seed = seed, 
                      iteration = iter, 
                      n_nest_fac = n_nest_fac, 
                      age = 0, 
                      fit,
                      # Metadata for evolved solutions
                      init_method = 2,  # 2 = evolved (not initial)
                      efa_rotation_code = 0,  # Not applicable
                      efa_n_factors = 0,  # Not applicable
                      efa_n_items_sampled = sum(item_assignment > 0))  # Current items in model
          
          new_sol
        }
        #new_sol
      }, simplify = "matrix")
      
    } else if (parallel){
      # Parallel implementation - similar to single-core but distributed across nodes
      
      # Export updated environment variables to all cluster nodes
      clusterExport(cl = myCL, varlist = ls(envir = environment()), envir = environment())
      clusterExport(cl = myCL, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
      
      # Process each bee's operation in parallel
      new_solutions <- parSapply(myCL, 1:nrow(allBees), function(i_bee){
        
        if (verbose) print(allBees[i_bee,])
        
        # Get the current item assignment for this bee to modify
        item_assignment_orig <- as.numeric(as.character(unlist(allBees[i_bee, item_names])))
        names(item_assignment_orig) <- item_names
        
        # Apply appropriate bee operation based on the bee's role
        if (allBees$role[i_bee] == "scout"){
          item_assignment <- scout.bee(item_assignment = item_assignment_orig,
                                       max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <-reorder_factors(item_assignment = item_assignment)
          
        } else if (allBees$role[i_bee] == "onlooker") {
          # Onlooker bees make minor model changes
          item_assignment <- onlooker.bee(item_assignment = item_assignment_orig, 
                                          max_nest_fac = max_nest_fac)
          
          # Ensure consistent factor numbering
          item_assignment <-reorder_factors(item_assignment = item_assignment)
        } else {stop("Fatal error: Bee swarm out of control.")} # if an unknown role occurs...
        
        # Check if this model configuration has been tested before
        compare_item_assignment <- apply(solutions[, item_names], 1, function(i_solution) identical(item_assignment, i_solution))

        if (any(compare_item_assignment)){ # if you saw this item configuration before
          if (verbose) print ("Empty Flower")
          new_sol <- solutions[max(which(compare_item_assignment)),]
          new_sol$age <- new_sol$age + 1 
        } else { # if you did not see this model before
          
          # Calculate fit metrics for this model
          fit <- fit.function(item_assignment = item_assignment, 
                              item_names = item_names,
                              dat = data,
                              verbose = verbose,
                              fit_crit = fit_crit,
                              logistic_weights = logistic_weights,
                              nu_weights = nu_weights,
                              nu_min = nu_min,
                              debug_fit_mode = debug_fit_mode,
							  ignore_warnings = ignore_warnings,
							  max_nest_fac = max_nest_fac,
                              ...
                              ) #evaluate model
          
          
          # Store model information and results
          n_nest_fac <- max(item_assignment)
          
          # include metadata for consistency
          new_sol <- c(unlist(item_assignment), 
                      seed = seed, 
                      iteration = iter, 
                      n_nest_fac = n_nest_fac, 
                      age = 0, 
                      fit,
                      # Metadata for evolved solutions
                      init_method = 2,  # 2 = evolved (not initial)
                      efa_rotation_code = 0,  # Not applicable
                      efa_n_factors = 0,  # Not applicable
                      efa_n_items_sampled = sum(item_assignment > 0))  # Current items in model
    
          new_sol
        }
      }, simplify = "matrix")
    }

    
    # Transform results to a data frame for easier processing
    new_solutions <- data.frame(t(new_solutions))
    

  #  new_solutions2 <- 
    # Update visualization if plotting is enabled
    if (plot_nectar) {
      # Create a copy of solutions for plotting
      new_solutions2 <- new_solutions
      # Add information needed for the plot
      new_solutions2$iter <- iter 
      new_solutions2$role <- as.factor(allBees$role)
      new_solutions2$n_nest_fac <- as.factor(new_solutions2$n_nest_fac)
      # Update plot and suppress NA-warnings (e.g., when fit is not converged)
      suppressWarnings(expr = {
        conv_plot <- conv_plot + 
          geom_jitter(data = new_solutions2, 
                      aes(x = iter,
                          y = fit_overall,
                          color = n_nest_fac,
                          shape = role),
                      width = plot_list$jitter_width,
                      alpha = plot_list$alpha,
                      size = plot_list$size) 
        
        # Ensure consistent shape scale across all plot layers
        if (!"shape" %in% names(conv_plot$scales$scales)) {
          conv_plot <- conv_plot + 
            scale_shape_manual(values = c("nest_site_scouts" = 17, "random_init" = 18, 
                                         "scout" = 16, "onlooker" = 15))
        }
        
       if(!cluster_mode) plot(conv_plot)
      })
      
    } #plot quality of new solution
    

    # Verify that at least some models have valid fit values
    
    if (all(is.na(new_solutions$fit_overall))) stop("All initial fit values are NA. Please check fit-function. Try setting verbose = TRUE or debug_fit_mode = TRUE.")

    # Sort new solutions by fit quality
    new_solutions <- new_solutions[order(new_solutions$fit_overall, decreasing = TRUE),]
    
    # Check if we've found a new best solution
    if (max(new_solutions$fit_overall, na.rm = TRUE) > best_fit) { # if we have a new winner...
      best_fit <- max(new_solutions$fit_overall, na.rm = TRUE)
      best_solution <- new_solutions[1,]
      best_solution$nCores <- nCores # add number of cores to best solution
      counter <- 0 # reset counter to zero because it will be counted +1 below
    } 
  
    # Combine old and new solutions, then sort by quality
    solutions <- robust_bind_rows(solutions, new_solutions)
    solutions <- solutions[order(solutions$fit_overall, decreasing = T),] #sort by quality
    
    # Update iteration counters
    counter <- counter + 1 #increase count timer
    iter <- iter + 1   #increase overall iteration counter
    end <- Sys.time() 

  } # End of main optimization loop
  
  # Clean up: stop cluster if parallel processing was used
  if (parallel) stopCluster(myCL)
  
  # Save results to files if requested
  if (write_solutions) {
    write.table(solutions, summaryfile, sep=";",row.names = F, col.names = T, append = F)
  }
  
  # Capture session information for reproducibility
  session_info <- capture.output(sessionInfo())
  
  # Write session information and best solution to final summary file
  writeLines(c("# Session Information:", session_info, "# End of Session Information", ""), con = summaryfile_fin)
  write.table(best_solution, summaryfile_fin, sep=";",row.names = F, col.names = T, append = F)
  
  # Return convergence plot if created
  if(plot_nectar) {
    conv_plot + xlim(0,iter)
  } else {
      NULL
    }
  
}
############################################# End BSO #################################################