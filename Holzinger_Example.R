  rm(list = ls()) # emtpy workspace
  
  source("BSO_main_function.R")
  # load Holzinger & Swineford data set
  data("HS.ability.data", package = "OpenMx")
  HS_data <- HS.ability.data
  item_names <- colnames(HS_data[,7:30])
  HS_data <- as.data.frame(apply(HS_data[, item_names], 2, scale))
  
  set.seed(456)
  
  # This is a small working example which will run several minutes (< 10 min.)
  # For more realistic settings of the hyperparameters please see text
  max_iter <- 5     		 # number of iterations without improvement after which to abort search
  n_bees <- 200     		 # number of n_bees
  nest_site_bees <- 200  # number of nest site bees for inital evaluation
  top_best <- n_bees*0.2 # how many solutions to use for onlookers (for example 20% of n_bees)
  percent_scouts <- .25  # percentage scouts
  min_nest_fac <- 1    	 # how many correlated factors minimum
  max_nest_fac <- 5    	 # how many correlated factors maximum
  depletion <- 5 			   # update best solutions after how many iterations
  cluster_mode <- F  # Are you running on cluster? If so, do not repeat plotting every iteration
  
 
  for (seed in 1) {
  
    res.name <- paste0("BSO_Holz_s" , seed , "_b" , n_bees, "_s", percent_scouts*n_bees, "_t", top_best, Sys.Date())
    summaryfile <- paste0(res.name, ".csv")
    summaryfile_fin = paste0(res.name, "_final.csv")
    
    conv_plot <- 
      BSO(item_names = item_names, 
          data = HS_data,
          max_iter = max_iter, 
          n_bees = n_bees,
          #n_start_bees = nest_site_bees,
          nest_site_bees = nest_site_bees,
          percent_scouts = percent_scouts, 
          top_best = top_best, 
          min_nest_fac = min_nest_fac, 
          max_nest_fac = max_nest_fac, 
          depletion = depletion, 
          summaryfile = summaryfile,
          summaryfile_fin = summaryfile_fin,
          seed = seed,
          verbose = F,
          parallel = T,
          nCores = 8, # Define the number of available cores 
          plot_nectar = TRUE,
          fit_crit = c("cfi", "rmsea",  "max_fac_cor1", "n_items","min_loading", "residual"),
          logistic_weights = list(c(d = 0.90, a = 30),
                                  c(d = 0.06, a = -50),
                                  c(d = 0.80, a = -30),     # Max factor correlation
                                  c(d = 0.95, a = 50),
                                  c(d = 0.30, a = 50),
                                  c(d = NA, a = NA)), #residual (local fit) wont be logit transformed a second time
         # nu_weights = c(1,.5,1,1,1,1),
          nu_weights = c(2,2,3,1,1,1), # Weighting matters, look at the comment examples below: This is close to the one in the paper
          #nu_weights = c(2,2,3,2,1,1), # stronger weights n_items (less item elimination) --> weights influence results
          nu_min = 0.00001,
          balance_n_fac = F,
          ignore_warnings  = F, # If True models which produce warnings won´t be penalized
          bounds = "pos.var", # additional lavaan agrument 
          cluster_mode = cluster_mode,
          plot_list = list(xlim = c(0, max_iter*10), 
                           #ylim = c(0, sum(c(1,1,.5,1,1,1))),
                           ylim = c(0, sum(c(2,2,3,1,1,1))),
                           ylab = "Overall Nectar Value",
                           xlab = "Iteration",
                           jitter_width = 0.5,
                           alpha = 0.2,
                           size = 1)
      )
    
    if (!is.null(conv_plot)){
      suppressWarnings(expr = {
        factor_cols <- hue_pal()(5)
        if(!cluster_mode) plot(conv_plot)
        ggsave(plot = conv_plot, filename = paste0(res.name,".pdf"), device = "pdf")
      })
    }
  }
  
