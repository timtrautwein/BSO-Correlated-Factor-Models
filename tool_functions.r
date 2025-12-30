################## BSO subfunctions #############################
#' @title Bee Swarm Optimization (BSO) Helper Functions
#' @description This file contains a collection of sub-functions which are utilized during BSO algorithm execution.
#' These functions handle model initialization, scout bee operations, onlooker bee operations, and various helper utilities.

#' @title Initialize Bee Solutions
#' @description Generates the random start models in the first iteration of the BSO algorithm.
#' 
#' @param item_names Vector of item names.
#' @param data Data frame containing the data set.
#' @param max_iter Maximum number of iterations.
#' @param n_bees Number of bees.
#' @param n_items Number of items (inferred from name vector by default).
#' @param verbose Logical; if TRUE, print detailed information.
#' @param debug_fit_mode Logical; if TRUE, return fit object for debugging.
#' @param fit_crit Vector of fit criteria to extract.
#' @param logistic_weights Weights for logistic transformation.
#' @param nu_weights Weights for the aggregation function.
#' @param nu_min Threshold value for fit criteria.
#' @param n_start_items How many items should be in the start model (allows starting with a subset of items).
#' @param balance_n_fac Logical; if TRUE, start models are sampled such that the number of factors is uniformly distributed, if FALSE, random item assignments are used resulting in more candidate models with fewer factors 
#' @param ... Additional arguments passed to lavaan.
#'
#' @return A vector containing item assignments, seed, iteration number, factor count, age, and fit statistics.
#'
#' @details
#' Item assignment is coded as follows:
#' - -1: Item is not in the model
#' -  1: max_nest_fac: Item loads on this correlated factor (nest is the name for consistency with the original function)
#'
#' The function ensures that each factor has at least 3 items.
#### Original Initialization function ####
# Generates the random start models in the first iteration
initialize.bees <- function(item_names, # a vector of item names (length must equal the length of item assignments)
                            data, # a data frame containing the data set
                            max_iter, # maximum number of iterations
                            n_bees, # number of bees
                            n_items = length(item_names), # number of items (inferred from name vector)
                            verbose = FALSE, # verbose mode (passed to fit function)
                            debug_fit_mode = FALSE, # debug mode (passed to fit function)
							ignore_warnings = FALSE, # whether warnings should be ignored in fit.function
                            fit_crit = c("cfi", "rmsea", "min_omega2", "min_loading"), # see fit-function()
                            logistic_weights, # see fit-function()
                            nu_weights, # see fit-function()
                            nu_min, # see fit-function()
                            n_start_items = length(item_names), # How many items should be in the start model? (allows start with a subselection of items)
                            balance_n_fac = FALSE, # if TRUE, the start models are sampled such that the number of factors is uniformly distributed
                            ...){ # further arguments to be passed to lavaan
  #####
  # Legacy initialization: generate random start models (pre Scout-Nest initialization).
  #####
  # item_assignment is a vector with numbers with the following coding scheme:
  # -1                ... Item is not in the model
  
  # 1:max_nest_fac    ... Item loads on this correlated factor (nest is the name for consistency with the original function)
  
  if (!balance_n_fac){ 
    # If not balancing factor counts, create random assignment
    # if we start with a subset of the items, compute how many are excluded
    n_items_excluded <- n_items - n_start_items
    
    # return error for impossible argument combinations
    if(n_items_excluded < 0 | (n_items - n_items_excluded) < 0) stop("Invalid value for n_start_items")
    
    # assign random item allocations
    item_assignment <- rep(-1, n_items_excluded)
    item_assignment <- c(item_assignment, resamp(1:max_nest_fac, n_items - n_items_excluded, replace = TRUE)) # generate models
    
    # Note that across all start models this procedure will result in a higher number of models
    # with fewer factors (because these are more frequent among all possible random assignments)

  } else {
    
    # This approach tries to balance the number of factors across all start models. 
    
    # draw a random number of factors first
    n_fac_tmp <- sample(min_nest_fac:max_nest_fac, size = 1)
    item_assignment <- c(rep(1:n_fac_tmp, each = 3)) # place at least 3 items for each factor
    
    # the other items can be assigned to any of the factors 
    n_items_excluded <- n_items - n_start_items
    n_items_remaining <- n_items - length(item_assignment) - n_items_excluded
    # return error for impossible argument combinations
    if(n_items_excluded < 0 | n_items_remaining < 0) stop("Invalid value for n_start_items")
    item_assignment <- c(item_assignment, rep(-1, n_items_excluded), sample(1:n_fac_tmp, size = n_items_remaining, replace = TRUE))
      
    # randomize order
    item_assignment <- sample(item_assignment, size = length(item_assignment))
  }
  
  # remove nested factors with less than 3 items
  item_counts <- table(item_assignment[item_assignment %in% 1:max_nest_fac])
  fac_2b_removed <- as.numeric(names(item_counts)[item_counts < 3])
  # remove the factor by removing these items 
  item_assignment[item_assignment %in% fac_2b_removed] <- -1 #Items will be removed instead of moving them to g factor
  
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)
  
  # Reorder Factors here
  item_assignment <-reorder_factors(item_assignment = item_assignment)
  
  # Name assignment so that the names are exported
  names(item_assignment) <- item_names
  
  # if there are items loading on g only, the max fac is the general factor
  # otherwise it is 
  n_nest_fac <- max(item_assignment)
  
  # fit the the model
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
  
  # save items - factor allocation, seed, iteration, number of factors, "age" and quality in solution object
  # Use consistent numeric metadata columns that match Scout-Nest-Bees format
  new_sol <- c(item_assignment, 
               seed = seed, 
               iteration = 0, 
               n_nest_fac = n_nest_fac, 
               age = 0, 
               fit,
               # Metadata for traceability - consistent with Scout-Nest-Bees format
               init_method = 0,  # 0 = random initialization (numeric)
               efa_rotation_code = 0,  # Not applicable for random init
               efa_n_factors = 0,  # Not applicable for random init  
               efa_n_items_sampled = sum(item_assignment > 0))  # Current items in model
  new_sol 
}

#' Calculate Geomin epsilon based on Mplus defaults
#' 
#' @param n_factors Number of factors
#' @return Numeric epsilon value
calculate_geomin_epsilon <- function(n_factors) {
  if (n_factors == 2) 1e-4      # 0.0001
  else if (n_factors == 3) 1e-3  # 0.001
  else 1e-2                      # 0.01 (≥4 factors)
}

#' Perform EFA for BSO initialization
#' 
#' @param data Standardized data frame
#' @param n_factors Number of factors to extract
#' @param rotation Rotation method: "geominQ", "geominQ_05", "equamax", or "promax"
#' @param item_names Vector of item names
#' @return psych::fa object or NULL if failed
perform.efa.for.bso <- function(data, n_factors, rotation, item_names) {
  tryCatch({
    if (rotation == "geominQ") {
      # Use Mplus default epsilon
      fa(data[, item_names],
         nfactors = n_factors,
         rotate = "geominQ",
         fm = "ml",
         delta = calculate_geomin_epsilon(n_factors),
         max.iter = 1000)
    } else if (rotation == "geominQ_05") {
      # Use epsilon = 0.5
      fa(data[, item_names],
         nfactors = n_factors,
         rotate = "geominQ",
         fm = "ml",
         delta = 0.5,
         max.iter = 1000)
    } else {
      # Equamax or Promax
      fa(data[, item_names],
         nfactors = n_factors,
         rotate = rotation,
         fm = "ml",
         max.iter = 1000)
    }
  }, error = function(e) {
    NULL  # Return NULL if EFA fails
  })
}

#' Convert EFA loadings to CFA item assignments
#' 
#' @param efa_result Result from perform.efa.for.bso
#' @param item_names All item names (including those not in EFA)
#' @param min_loading Minimum loading required for item to be retained
#' @return Named vector of item assignments (-1 = excluded, 1:k = factor)
efa.to.cfa.assignment <- function(efa_result, item_names, min_loading = 0) {
  if (is.null(efa_result)) {
    # If EFA failed, return all items excluded
    assignment <- rep(-1, length(item_names))
    names(assignment) <- item_names
    return(assignment)
  }
  
  # Get loadings matrix
  loadings <- efa_result$loadings[]
  used_items <- rownames(loadings)
  
  # Initialize assignment with -1 (excluded) for all items
  assignment <- rep(-1, length(item_names))
  names(assignment) <- item_names
  
  # For each item in the EFA, find its highest loading
  for (item in used_items) {
    item_loadings <- abs(loadings[item, ])
    max_loading <- max(item_loadings)
    
    # Check if item meets minimum loading requirement
    if (max_loading >= min_loading) {
      # Assign to factor with highest loading
      # If tie, randomly choose
      max_factors <- which(item_loadings == max_loading)
      if (length(max_factors) > 1) {
        assigned_factor <- sample(max_factors, 1)
      } else {
        assigned_factor <- max_factors
      }
      assignment[item] <- assigned_factor
    }
  }
  
  # Check that each factor has at least 3 items
  factor_counts <- table(assignment[assignment > 0])
  for (fac in names(factor_counts)) {
    if (factor_counts[fac] < 3) {
      # Remove this factor by setting items to -1
      assignment[assignment == as.numeric(fac)] <- -1
    }
  }
  
  # Renumber factors to be consecutive
  assignment <- renew_factor_numbers(assignment)
  assignment <- reorder_factors(assignment)
  
  return(assignment)
}

#' Scout-Nest-Bees Initialization using EFA
#' 
#' @description Initializes BSO using EFA with different rotation methods and item sampling strategies.
#' These nest site scout bees explore the solution space to find promising initial models (nest sites)
#' before the main optimization begins.
#' 
#' @param item_names Vector of item names
#' @param data Data frame containing the dataset
#' @param nest_site_bees Number of nest site scout bees to generate for initial exploration
#' @param min_nest_fac Minimum number of factors
#' @param max_nest_fac Maximum number of factors
#' @param n_items_min Minimum number of items in the itemset (0 = no minimum)
#' @param min_efa_loading Minimum EFA loading for items to be retained
#' @param fit_crit Vector of fit criteria
#' @param logistic_weights Weights for logistic transformation
#' @param nu_weights Weights for aggregation function
#' @param nu_min Threshold value for nu
#' @param debug_fit_mode Debug mode flag
#' @param verbose Verbose output flag
#' @param seed Random seed
#' @param ... Additional arguments for lavaan
#' 
#' @return Data frame of solutions with failed attempts information for fallback
initialize.scout.nest.bees <- function(item_names,
                                      data,
                                      nest_site_bees,
                                      min_nest_fac,
                                      max_nest_fac,
                                      n_items_min = 0,
                                      min_efa_loading = 0,
                                      fit_crit,
                                      logistic_weights,
                                      nu_weights,
                                      nu_min,
                                      debug_fit_mode = FALSE,
                                      verbose = FALSE,
									  ignore_warnings = FALSE, # whether warnings should be ignored in fit.function
                                      seed = 1,
                                      ...) {
  
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  n_items <- length(item_names)
  
  # Define rotation methods
  rotation_methods <- c("geominQ", "geominQ_05", "equamax", "promax")
  n_rotations <- length(rotation_methods)
  
  # Define factor range
  factor_range <- min_nest_fac:max_nest_fac
  n_factor_levels <- length(factor_range)
  
  # Define item sampling percentiles
  if (n_items_min > 0 && n_items_min < n_items) {
    # Use provided minimum
    item_percentiles <- c(
      min = n_items_min,
      p25 = round(n_items_min + 0.25 * (n_items - n_items_min)),
      p50 = round(n_items_min + 0.50 * (n_items - n_items_min)),
      p75 = round(n_items_min + 0.75 * (n_items - n_items_min)),
      p100 = n_items
    )
  } else {
    # Use default minimum of 4-5 items for EFA
    min_efa_items <- max(5, min_nest_fac * 3)  # At least 3 items per factor
    item_percentiles <- c(
      min = min_efa_items,
      p25 = round(min_efa_items + 0.25 * (n_items - min_efa_items)),
      p50 = round(min_efa_items + 0.50 * (n_items - min_efa_items)),
      p75 = round(min_efa_items + 0.75 * (n_items - min_efa_items)),
      p100 = n_items
    )
  }
  item_percentiles <- unique(item_percentiles)
  
  # Generate all possible combinations
  all_combinations <- expand.grid(
    n_factors = factor_range,
    rotation = rotation_methods,
    n_items_sample = item_percentiles,
    stringsAsFactors = FALSE
  )
  
  # Calculate how many bees we need for basic coverage
  n_basic_combinations <- nrow(all_combinations)
  
  if (verbose) {
    cat("Scout-Nest-Bees initialization:\n")
    cat("Total possible combinations:", n_basic_combinations, "\n")
    cat("Requested nest site scout bees:", nest_site_bees, "\n")
  }
  
  # Initialize solutions list and failed attempts tracker
  solutions_list <- list()
  failed_attempts <- list()  # Track failed EFA attempts for fallback
  
  # For 100% item sampling, we only need one solution per factor/rotation combo
  # because the results will be identical
  unique_100_percent <- all_combinations[all_combinations$n_items_sample == n_items, ]
  other_combos <- all_combinations[all_combinations$n_items_sample < n_items, ]
  
  # Track which combinations have been used and their results
  efa_cache <- list()
  
  # Strategy for allocating bees to combinations
  if (nest_site_bees <= nrow(unique_100_percent)) {
    # If we have fewer bees than unique 100% combinations, just use those
    sampled_combinations <- unique_100_percent[1:nest_site_bees, ]
  } else if (nest_site_bees <= n_basic_combinations) {
    # Use all unique 100% combinations, then sample from others
    sampled_combinations <- rbind(
      unique_100_percent,
      other_combos[1:min(nest_site_bees - nrow(unique_100_percent), nrow(other_combos)), ]
    )
  } else {
    # We have more bees than basic combinations
    # First use all unique combinations
    sampled_combinations <- all_combinations
    
    # Then add replicates from non-100% sampling (which will have random item selection)
    remaining_bees <- nest_site_bees - nrow(sampled_combinations)
    if (remaining_bees > 0 && nrow(other_combos) > 0) {
      # Sample with replacement from non-100% combinations
      additional_indices <- sample(1:nrow(other_combos), remaining_bees, replace = TRUE)
      sampled_combinations <- rbind(sampled_combinations, other_combos[additional_indices, ])
    }
  }
  
  # Ensure we have exactly nest_site_bees
  sampled_combinations <- sampled_combinations[1:min(nest_site_bees, nrow(sampled_combinations)), ]
  
  if (verbose) {
    cat("Generating", nrow(sampled_combinations), "nest site scout bees...\n")
  }
  
  # Generate solutions for each sampled combination
  for (i in 1:nrow(sampled_combinations)) {
    combo <- sampled_combinations[i, ]
    
    # Sample items if not using all
    if (combo$n_items_sample < n_items) {
      sampled_items <- sample(item_names, combo$n_items_sample, replace = FALSE)
    } else {
      sampled_items <- item_names
    }
    
    if (verbose) {
      cat(sprintf("Nest site bee %d: %d factors, %s rotation, %d items\n", 
                  i, combo$n_factors, combo$rotation, combo$n_items_sample))
    }
    
    # Create cache key for EFA results
    cache_key <- paste(combo$n_factors, combo$rotation, combo$n_items_sample, 
                      paste(sort(sampled_items), collapse = "_"), sep = "_")
    
    # Check if we've already computed this exact EFA
    if (cache_key %in% names(efa_cache)) {
      # Use cached assignment
      item_assignment <- efa_cache[[cache_key]]
      if (verbose) cat("  -> Using cached EFA result\n")
    } else {
      # Perform EFA
      efa_result <- perform.efa.for.bso(data, combo$n_factors, combo$rotation, sampled_items)
      
      # Convert to item assignment
      item_assignment <- efa.to.cfa.assignment(efa_result, item_names, min_efa_loading)
      
      # Cache the result
      efa_cache[[cache_key]] <- item_assignment
    }
    
    # Skip if no valid assignment
    if (all(item_assignment == -1)) {
      if (verbose) cat("  -> EFA failed or no items met criteria\n")
      
      # Track failed attempt for potential fallback use
      failed_attempts[[length(failed_attempts) + 1]] <- list(
        n_factors = combo$n_factors,
        rotation = combo$rotation,
        n_items_sample = combo$n_items_sample,
        sampled_items = if (combo$n_items_sample < n_items) sampled_items else item_names
      )
      
      next
    }
    
    # Fit the model
    n_nest_fac <- max(item_assignment)
    
    fit <- fit.function(item_assignment = item_assignment,
                       item_names = item_names,
                       dat = data,
                       verbose = FALSE,
                       fit_crit = fit_crit,
                       logistic_weights = logistic_weights,
                       nu_weights = nu_weights,
                       nu_min = nu_min,
                       debug_fit_mode = debug_fit_mode,
					   ignore_warnings = ignore_warnings,
					   max_nest_fac = max_nest_fac,
                       ...)
    
    # Create solution vector - ensure all numeric types and consistent metadata
    # Encode rotation method as numeric
    rotation_code <- match(combo$rotation, rotation_methods)  # 1-4 for the four methods
    
    new_sol <- c(item_assignment,
                seed = seed,
                iteration = 0,
                n_nest_fac = n_nest_fac,
                age = 0,
                fit,
                # Metadata for Scout-Nest-Bees initialization (all numeric)
                init_method = 1,  # 1 = scout_nest_efa initialization
                efa_rotation_code = rotation_code,  # 1-4 for rotation methods
                efa_n_factors = combo$n_factors,
                efa_n_items_sampled = combo$n_items_sample)
    
    solutions_list[[length(solutions_list) + 1]] <- new_sol
  }
  
  # Convert to data frame - ensure numeric types are preserved
  if (length(solutions_list) == 0) {
    stop("No valid solutions generated from Scout-Nest-Bees initialization")
  }
  
  # Combine solutions into matrix first (preserves numeric types)
  solutions_matrix <- do.call(rbind, solutions_list)
  
  # Convert to data frame
  solutions <- as.data.frame(solutions_matrix)
  
  # Ensure all columns that should be numeric are numeric
  # All columns except row names should be numeric at this point
  solutions[] <- lapply(solutions, function(x) as.numeric(as.character(x)))
  
  if (verbose) {
    cat("Generated", nrow(solutions), "valid solutions\n")
    if (length(failed_attempts) > 0) {
      cat("Failed attempts:", length(failed_attempts), "\n")
    }
  }
  
  # Return solutions with failed attempts info for fallback
  attr(solutions, "failed_attempts") <- failed_attempts
  return(solutions)
}

#' Robust data frame binding that handles mismatched columns
#' 
#' @param df1 First data frame
#' @param df2 Second data frame  
#' @param fill_value Value to use for missing columns (default NA)
#' @return Combined data frame
robust_bind_rows <- function(df1, df2, fill_value = NA) {
  # Get all unique column names
  all_cols <- unique(c(names(df1), names(df2)))
  
  # Add missing columns to df1
  missing_in_df1 <- setdiff(all_cols, names(df1))
  for (col in missing_in_df1) {
    # Use appropriate fill values for metadata columns
    if (col %in% c("init_method", "efa_rotation_code", "efa_n_factors", "efa_n_items_sampled")) {
      df1[[col]] <- 0  # Use 0 for metadata columns that should be numeric
    } else {
      df1[[col]] <- fill_value
    }
  }
  
  # Add missing columns to df2  
  missing_in_df2 <- setdiff(all_cols, names(df2))
  for (col in missing_in_df2) {
    # Use appropriate fill values for metadata columns
    if (col %in% c("init_method", "efa_rotation_code", "efa_n_factors", "efa_n_items_sampled")) {
      df2[[col]] <- 0  # Use 0 for metadata columns that should be numeric
    } else {
      df2[[col]] <- fill_value
    }
  }
  
  # Reorder columns to match
  df1 <- df1[, all_cols]
  df2 <- df2[, all_cols]
  
  # Now rbind should work
  return(rbind(df1, df2))
}

#' @title Scout Bee Function
#' @description Performs major model modifications as a scout bee in the BSO algorithm.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A modified item assignment vector after applying a major operation.
#'
#' @details Scout bees perform major operations on models:
#' 1. Add a new factor
#' 2. Split an existing factor
#' 3. Remove a factor
#' 4. Merge two factors
#' 
#' The function randomly selects one of these operations based on what's possible 
#' with the current model configuration.
scout.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment, levels = 1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  # Check which major changes are allowed and choose one randomly 
  allowed_operations <- check.operations.scouts(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add new factor
            free_items <- which(item_assignment <= 0) # Which items are available?
            n_new_items <- resamp(3:length(free_items), size = 1) # How many items should the new factor have?
            new_fac_items <- resamp(free_items, size = n_new_items) # Choose items randomly
            item_assignment[new_fac_items] <- n_nest_fac + 1}, # Assign these items to the new factor
          { #split factors
            big_fac <- which(n_items_per_fac >= 6) # Which factors could be split?
            split_fac <- resamp(big_fac, size = 1) # Decide which factor should be split
            split_fac_items <- which(item_assignment == split_fac) # Which items belong to the factors?
            n_new_fac_items <- resamp(3:(length(split_fac_items) - 3), size = 1) # Split the factor randomly such that each factor has at least 3 items
            new_fac_items <- resamp(split_fac_items, size = n_new_fac_items) # Sample items for new factor
            
            item_assignment[new_fac_items]<- n_nest_fac + 1},# Assign these items to the new factor
          { #remove factor
            old_fac <- resamp(1:n_nest_fac, size = 1)
            item_assignment[item_assignment == old_fac] <- -1 }, # Remove the Factor including Items
          { #merge factors
            old_fac <- resamp(1:n_nest_fac, size = 2)
            item_assignment[item_assignment == old_fac[1]] <- old_fac[2] }
  )
  
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  item_assignment <- renew_factor_numbers(item_assignment = item_assignment)   
  item_assignment
}

#' @title Onlooker Bee Function
#' @description Performs minor model modifications as an onlooker bee in the BSO algorithm.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A modified item assignment vector after applying a minor operation.
#'
#' @details Onlooker bees perform minor operations on models:
#' 1. Add an item to a factor
#' 2. Remove an item from a factor
#' 3. Swap items between factors
#' 4. Delete an item from the item pool (redundant for correlated factors with operation 2)
#' 
#' The function randomly selects one of these operations based on what's possible 
#' with the current model configuration.
onlooker.bee <- function(item_assignment, max_nest_fac){
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)

  # Check which major changes are allowed and choose one randomly
  allowed_operations <- check.operations.onlookers(item_assignment = item_assignment,
                                                max_nest_fac = max_nest_fac)

  chosen_operation <- resamp(allowed_operations, size = 1) 

  switch (chosen_operation,
          { # add item to nested factor
            which_item    <- resamp(which(item_assignment <= 0), size = 1) # look for items not assigned to a nested factor
            to_which_factor   <- resamp(1:n_nest_fac, 1) # sample a new factor for this item
            item_assignment[which_item] <- to_which_factor}, # change assignment
          { # remove item from nested factor
            which_factor   <- resamp(which(n_items_per_fac > 3), size = 1) # look for a factor with more than the minimal item count
            which_item   <- resamp(which(item_assignment == which_factor), size = 1) # sample any item from this factor
            item_assignment[which_item] <- -1}, # remove item from factor
          { # swap item between nested factors
            item_1 <- resamp(which(item_assignment > 0), size = 1) # sample any item assigned to a nested factor
            item_2 <- resamp(which((item_assignment > 0) & (item_assignment != item_assignment[[item_1]])), size = 1) # sample any item from another nested factor
            item_assignment[c(item_1, item_2)] <- item_assignment[c(item_2, item_1)]},
          { # delete item from item pool #redundant with 2
            candidate_factors   <- c(which(n_items_per_fac > 3)) # look for a factor with more than the minimal item count 
            which_item   <- resamp(which(item_assignment %in% candidate_factors), size = 1) # sample any item from this factor or pick an item from the general factor
            item_assignment[which_item] <- -1} # remove item from item pool completely
  )
  item_assignment # return
}

#' @title Scout Bee Operation Selection
#' @description Decides which scout operations are possible based on an item assignment and the number of factors.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A vector of indices corresponding to allowed operations.
#'
#' @details Checks which of these major operations can be performed:
#' 1. Add a new factor
#' 2. Split a factor
#' 3. Remove a factor
#' 4. Merge factors
check.operations.scouts <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_fac_allowed = (n_nest_fac < max_nest_fac) & (sum(item_assignment <= 0) >= 3),     #can scouts add factors?
    spl_fac_allowed = n_nest_fac < max_nest_fac & any(n_items_per_fac >= 6),  #can scouts split factors?
    rmv_fac_allowed = n_nest_fac > min_nest_fac,                      #can scouts remove factors?
    mer_fac_allowed = n_nest_fac > min_nest_fac                      #can scouts merge factors?
  )
   which(operations)
}

#' @title Onlooker Bee Operation Selection
#' @description Decides which onlooker operations are possible based on an item assignment and the number of factors.
#' 
#' @param item_assignment Current item assignment vector.
#' @param max_nest_fac Maximum number of correlated factors allowed.
#'
#' @return A vector of indices corresponding to allowed operations.
#'
#' @details Checks which of these minor operations can be performed:
#' 1. Add an item to a factor
#' 2. Remove an item from a factor
#' 3. Swap items between factors
check.operations.onlookers <- function(item_assignment, max_nest_fac){
  
  n_items_per_fac <- table(factor(item_assignment,levels=1:max_nest_fac))
  n_nest_fac   <- max(item_assignment)
  
  operations <- c(
    add_item_allowed = any(item_assignment <= 0),         # can onlookers add item-factor allocations?
    rmv_item_allowed = any(n_items_per_fac > 3),   # can onlookers remove item-nested-factor allocations?
    swap_item_allowed = n_nest_fac > 1         # can onlookers swap item-nested-factor allocations?
  )
  
   which(operations)
}

#### Minor helper functions ####
#' @title Renew Factor Numbers
#' @description Ensures correct sequential factor numbering.
#' 
#' @param item_assignment Current item assignment vector.
#'
#' @return A vector with corrected sequential factor numbering.
#'
#' @details
#' Example: if a factor 5 occurs but factor 2 was removed,
#' the factors should be relabeled as 1:4 for consistency.
renew_factor_numbers <- function(item_assignment){
  # example: it is possible that a factor 5 occurs but factor 2 was removed
  # in that case, the factors should be relabeled 1:4
  # which items load on nested factors?
  items_nested <- item_assignment > 0
  item_assignment[items_nested] <- match(item_assignment[items_nested], sort(unique(item_assignment[items_nested])))
  item_assignment  
}

#' @title Robust Sample Function
#' @description Sample function that will also work correctly with vectors of length 1.
#' 
#' @param y Vector to sample from.
#' @param ... Additional arguments passed to sample().
#'
#' @return Sampled value(s) from the input vector.
resamp <- function(y,...){if(length(y)==1) y else sample(y,...)} 

#' Calculate maximum plausible factors for given number of items
#' 
#' @param n_items Number of items available
#' @param min_items_per_factor Minimum items per factor (default = 3 for CFA)
#' @return Maximum number of factors that can be reasonably supported
calculate.max.plausible.factors <- function(n_items, min_items_per_factor = 3) {
  return(floor(n_items / min_items_per_factor))
}

#' Generate random assignment matching failed EFA characteristics with plausible factor count
#' 
#' This function creates random item assignments that use the same items as failed EFA attempts
#' but calculates a realistic number of factors based on CFA requirements (≥3 items per factor).
#' This prevents mathematically impossible assignments that would fail in CFA fitting.
#' 
#' @param failed_attempt List containing information about failed EFA attempt
#' @param item_names Vector of all item names
#' @param min_nest_fac Minimum number of factors
#' @param max_nest_fac Maximum number of factors
#' @return Named vector of item assignments with realistic factor structure
generate.random.assignment.from.failed.efa <- function(failed_attempt, item_names, min_nest_fac, max_nest_fac) {
  # Use the same items that were sampled for the failed EFA
  target_items <- failed_attempt$sampled_items
  n_items_to_use <- length(target_items)
  
  # Calculate PLAUSIBLE number of factors based on available items
  # Rule: At least 3 items per factor for stable CFA
  max_plausible_factors <- calculate.max.plausible.factors(n_items_to_use, min_items_per_factor = 3)
  
  # Choose a realistic number of factors (not the failed EFA's unrealistic count)
  if (max_plausible_factors < min_nest_fac) {
    # Not enough items even for minimum factors - use minimum and add more items if possible
    n_factors <- min_nest_fac
    min_items_needed <- n_factors * 3
    
    if (length(item_names) >= min_items_needed) {
      # Add more items from the full set to make it work
      additional_items_needed <- min_items_needed - n_items_to_use
      additional_items <- sample(setdiff(item_names, target_items), additional_items_needed)
      target_items <- c(target_items, additional_items)
      n_items_to_use <- length(target_items)
    } else {
      # Still not enough items in total - use all available
      target_items <- item_names
      n_items_to_use <- length(target_items)
      n_factors <- max(1, floor(n_items_to_use / 3))
    }
  } else {
    # We can support some factors - choose a reasonable number
    n_factors <- min(max_plausible_factors, max_nest_fac)
    n_factors <- max(n_factors, min_nest_fac)
  }
  
  # Initialize assignment with all items excluded
  item_assignment <- rep(-1, length(item_names))
  names(item_assignment) <- item_names
  
  # Assign items to factors - ensure at least 3 items per factor
  items_to_assign <- target_items[1:min(length(target_items), n_items_to_use)]
  
  # First, assign exactly 3 items to each factor
  for (f in 1:n_factors) {
    start_idx <- (f - 1) * 3 + 1
    end_idx <- f * 3
    if (end_idx <= length(items_to_assign)) {
      factor_items <- items_to_assign[start_idx:end_idx]
      item_assignment[factor_items] <- f
    }
  }
  
  # Assign remaining items randomly to factors
  remaining_items <- items_to_assign[-(1:min(n_factors * 3, length(items_to_assign)))]
  if (length(remaining_items) > 0) {
    random_factors <- sample(1:n_factors, length(remaining_items), replace = TRUE)
    item_assignment[remaining_items] <- random_factors
  }
  
  # Clean up assignment
  item_assignment <- renew_factor_numbers(item_assignment)
  item_assignment <- reorder_factors(item_assignment)
  
  return(item_assignment)
}

#' @title Debug Helper
#' @description Debugging function to copy all temporary objects into the global environment.
#' 
#' @param tmp_env The environment to copy from (default is the current environment).
#'
#' @return Nothing, but copies all objects from the environment to the global environment.
allglobal <- function(tmp_env = environment()) {
  lss <- ls(envir = tmp_env)
  for (i in lss) {
    assign(i, get(i, envir = tmp_env), envir = globalenv())
  }
}

#' @title Reorder Factors
#' @description Reorders factor numbers based on their average positions in the item vector.
#' 
#' @param item_assignment Current item assignment vector.
#'
#' @return A vector with reordered factor numbering.
#'
#' @details
#' This improved function reorders factor numbers so that factors with items 
#' appearing earlier in the item vector tend to get lower factor numbers.
#' This creates more consistent labeling across multiple runs.
reorder_factors <- function(item_assignment) {
  # Remove -1 and other special values from consideration
  factors <- unique(item_assignment[item_assignment != -1])
  
  # Calculate average indices for each factor
  index_means <- sapply(factors, function(f) {
    mean(which(item_assignment == f))
  }, simplify = TRUE, USE.NAMES = TRUE)
  
  # Order factors by their average indices
  ordered_factors <- factors[order(index_means)]
  
  # Create a mapping from old factor numbers to new ordered factor numbers
  factor_mapping <- setNames(seq_along(ordered_factors), ordered_factors)
  
  # Apply the mapping to item_assignment to reorder factor labels
  item_assignment <- vapply(item_assignment, function(x) {
    if (x %in% names(factor_mapping)) {
      as.integer(factor_mapping[[as.character(x)]])
    } else {
      as.integer(x)  # Convert special values to integer explicitly
    }
  }, FUN.VALUE = integer(1))
  
  return(item_assignment)
}

##' Resolve factor-range arguments (min/max factors) with backward compatibility
##'
##' The code base historically used `min_nest_fac`/`max_nest_fac` to denote the allowed factor range.
##' Reason: Users often prefer the shorter naming `min_fac`/`max_fac`.
##'
##' This helper:
##' - accepts either naming scheme,
##' - checks for conflicts if both are provided, and
##' - returns a validated, consistent factor range.
##'
##' @param min_fac Optional integer. Minimum number of factors (preferred alias).
##' @param max_fac Optional integer. Maximum number of factors (preferred alias).
##' @param min_nest_fac Optional integer. Minimum number of factors (legacy name).
##' @param max_nest_fac Optional integer. Maximum number of factors (legacy name).
##' @return List with elements `min_nest_fac` and `max_nest_fac` (both integers).
#resolve_factor_range <- function(min_fac = NULL,
#                                 max_fac = NULL,
#                                 min_nest_fac = NULL,
#                                 max_nest_fac = NULL) {
#  # Detect conflicts if both naming schemes are used simultaneously
#  if (!is.null(min_fac) && !is.null(min_nest_fac) && (min_fac != min_nest_fac)) {
#    stop("Conflicting arguments: 'min_fac' and 'min_nest_fac' differ. Please provide only one or ensure they match.")
#  }
#  if (!is.null(max_fac) && !is.null(max_nest_fac) && (max_fac != max_nest_fac)) {
#    stop("Conflicting arguments: 'max_fac' and 'max_nest_fac' differ. Please provide only one or ensure they match.")
#  }
#  
#  min_out <- if (!is.null(min_fac)) min_fac else min_nest_fac
#  max_out <- if (!is.null(max_fac)) max_fac else max_nest_fac
#  
#  if (is.null(min_out) || is.null(max_out)) {
#    stop("You must provide a factor range via either (min_fac, max_fac) or (min_nest_fac, max_nest_fac).")
#  }
#  
#  min_out <- as.integer(min_out)
#  max_out <- as.integer(max_out)
#  
#  if (is.na(min_out) || is.na(max_out) || min_out < 1 || max_out < 1) {
#    stop("Factor range values must be positive integers.")
#  }
#  if (min_out > max_out) {
#    stop("Invalid factor range: min > max.")
#  }
#  
#  list(min_nest_fac = min_out, max_nest_fac = max_out)
#}#