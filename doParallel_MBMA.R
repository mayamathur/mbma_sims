
# IMPORTANT NOTES -----------------------------

# Important things to remember: 
#
# - The returned Vhat is an estimate of T2 + t2w, *not* T2 itself

# Debugging help:
# 
# - The jobs may fail before fitting modAll with no apparent errors if 
#   k is too large for rma.mv. In that case, try setting p$k < 500 for modAll
#  and modPub small to prevent those models from being fit. 


# for interactive Sherlock:
# path = "/home/groups/manishad/MBMA"
# setwd(path)
# source("doParallel_MBMA.R")


# because Sherlock 2.0 restores previous workspace
rm( list = ls() )


# are we running locally?
run.local = TRUE

# should we set scen params interactively on cluster?
interactive.cluster.run = FALSE

# should lots of output be printed for each sim rep?
verbose = FALSE

# ~~ Packages -----------------------------------------------
toLoad = c("crayon",
           "dplyr",
           "foreach",
           "doParallel",
           "boot",
           "metafor",
           "robumeta",
           "data.table",
           "purrr",
           "metRology",
           "fansi",
           "MetaUtility",
           "ICC",
           "cfdecomp",
           "tidyr",
           "truncdist",
           "tibble",
           "tmvtnorm",
           "testthat",
           "truncreg",
           "truncnorm",
           "rstan",
           "optimx",
           "weightr")

if ( run.local == TRUE | interactive.cluster.run == TRUE ) toLoad = c(toLoad, "here")


# SET UP FOR CLUSTER OR LOCAL RUN ------------------------------

# ~~ Cluster Run ----------------------------------------
if (run.local == FALSE) {
  
  # load command line arguments
  args = commandArgs(trailingOnly = TRUE)
  
  cat("\n\n args received from sbatch file:", args)
  
  jobname = args[1]
  scen = args[2]  # this will be a number
  
  # load packages with informative messages if one can't be installed
  # **Common reason to get the "could not library" error: You did ml load R/XXX using an old version
  any.failed = FALSE
  for (pkg in toLoad) {
    
    cat( paste("\nAbout to try loading package", pkg) )
    
    tryCatch({
      # eval below needed because library() will otherwise be confused
      # https://www.mitchelloharawild.com/blog/loading-r-packages-in-a-loop/
      eval( bquote( library( .(pkg) ) ) )
    }, error = function(err) {
      cat( paste("\n*** COULD NOT LIBRARYIZE PACKAGE:", pkg) )
      any.failed <<- TRUE
    })
    
  }
  if ( any.failed == TRUE ) stop("Some packages couldn't be installed. See outfile for details of which ones.")
  
  # helper code
  path = "/home/groups/manishad/MBMA"
  setwd(path)
  source("helper_MBMA.R")
  
  # ~~ Cluster Run ----------------------------------------
  
  if ( interactive.cluster.run == FALSE ) {
    # get scen parameters (made by genSbatch.R)
    setwd(path)
    scen.params = read.csv( "scen_params.csv" )
    p <<- scen.params[ scen.params$scen == scen, ]
    
    cat("\n\nHEAD OF ENTIRE SCEN.PARAMS:\n")
    print(p)
  }
  
  # ~~ Interactive Cluster Run ----------------------------------------
  # alternatively, generate a simple scen.params in order to run doParallel manually in
  # Sherlock as a test
  if ( interactive.cluster.run == TRUE ) {
    path = "/home/groups/manishad/MBMA"
    setwd(path)
    source("helper_MBMA.R")
    
    scen.params = tidyr::expand_grid(
      # full list (save):
      # rep.methods = "naive ; gold-std ; pcurve ; maon ; 2psm ; rtma ; jeffreys-sd ; jeffreys-var ; mle-sd ; mle-var ; MBMA-mle-sd ; 2psm-MBMA-dataset ; prereg-naive",
      #rep.methods = "naive ; rtma ; 2psm",
      rep.methods = "naive ; 2psm ; sapb ; rtma",
      
      # args from sim_meta_2
      Nmax = 10,
      true.dist = "expo",
      Mu = c(0.5),
      t2a = c(0.2^2),
      t2w = c(0.2^2),
      m = 50,
      
      hack = c("affirm2"),
      rho = c(0),
      k.pub.nonaffirm = c(15),
      prob.hacked = c(0),
      
      eta = 5,
      
      true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
      
      muB = log(2),
      sig2B = 0.5,
      prob.conf = 0.5,
      
      # Stan control args
      stan.maxtreedepth = 20,
      stan.adapt_delta = 0.80,
      
      get.CIs = TRUE,
      run.optimx = FALSE )
    
    
    scen.params$scen = 1:nrow(scen.params)
    
    scen = 1
  }  # end "if ( interactive.cluster.run == TRUE )"
  
  
  
  # locally, with total k = 100, Nmax = 10, and sim.reps = 250, took 93 min total
  # for that I did sim.reps = 100 per doParallel
  
  # simulation reps to run within this job
  # **this need to match n.reps.in.doParallel in the genSbatch script
  if ( interactive.cluster.run == FALSE ) sim.reps = 2000
  if ( interactive.cluster.run == TRUE ) sim.reps = 1  
  
  # set the number of cores
  registerDoParallel(cores=16)
  
}



# FOR LOCAL USE  ------------------------------
if ( run.local == TRUE ) {
  #rm(list=ls())
  
  lapply( toLoad,
          require,
          character.only = TRUE)
  
  
  # helper fns
  code.dir = here()
  setwd(code.dir)
  source("helper_MBMA.R")
  
  
  # ~~ ********** Set Local Sim Params -----------------------------
  
  # for testing selection on SEs
  scen.params = tidyr::expand_grid(
    
    # for local runs, note that beta-sm is quite slow
    rep.methods = "naive ; mbma-MhatB ; mbma-muB ; 2psm ; beta-sm ; mbma-MhatB-gamma",
    
    
    # args from sim_meta_2
    Nmax = 1, 
    true.dist = "expo",
    Mu = c(0.25),
    t2a = c(0.25), 
    t2w = c(0),
    m = 50,
    
    #bm: try running existing hacking methods? affirm or affirm2?
    # remember that affirm will only work if prob.hacked < 1 else will never have nonaffirms
    #hack = c("favor-lowest-p"),
    #hack = "favor-gamma-ratio",
    hack = "affirm2",
    rho = c(0),
    k.pub.nonaffirm = c(30),
    prob.hacked = c(0), 
    
    eta = c(3),
    # within-study selection ratio; only used when hack=favor-gamma-ratio
    # important: other hacking methods will IGNORE gamma because can't be directly specified
    gamma = c(2), 
    SAS.type = "2psm", # "2psm" (original) or "carter"
    
    true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"), 
    
    # confounding parameters
    muB = 0.25,
    sig2B = 0.5,
    prob.conf = c(0.5),
    # muB = 0,
    # sig2B = 0,
    # prob.conf = 0,
    
    # Stan control args - only relevant if running RTMA
    stan.maxtreedepth = 25,
    stan.adapt_delta = 0.995,
    
    get.CIs = TRUE,
    run.optimx = FALSE )
  
  
  # more scens
  # scen.params = tidyr::expand_grid(
  #   
  #   rep.methods = "naive ; mbma-MhatB ; mbma-MhatB-true-t2 ; maon-adj-MhatB ; 2psm",
  #   
  #   
  #   # args from sim_meta_2
  #   Nmax = 1,
  #   Mu = c(0.5),
  #   t2a = c(0, 0.25^2, 0.5^2), 
  #   t2w = c(0),
  #   m = 50,
  #   
  #   hack = c("affirm"),  # but there will not be any hacking since prob.hacked is 0
  #   rho = c(0),
  #   k.pub.nonaffirm = c(5, 10, 15, 30, 50),
  #   prob.hacked = c(0),
  #   
  #   eta = c(1, 5, 10),
  #   
  #   true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
  #   
  #   # confounding parameters
  #   muB = log(1.5, 3),
  #   sig2B = 0.5,
  #   prob.conf = c(0.5), 
  #   
  #   # Stan control args - only relevant if running RTMA
  #   stan.maxtreedepth = 25,
  #   stan.adapt_delta = 0.995,
  #   
  #   get.CIs = TRUE,
  #   run.optimx = FALSE )
  
  
  
  scen.params$scen = 1:nrow(scen.params)
  
  
  sim.reps = 100  # reps to run in this iterate
  
  # set the number of local cores
  registerDoParallel(cores=8)
  
  scen = 1
  # data.frame(scen.params %>% filter(scen.name == scen))
  
  # just to avoid errors in doParallel script below
  jobname = "job_1"
  i = 1
}



# COMPILE STAN MODEL ONCE AT BEGINNING------------------------------

# only needed if running RTMA
# if ( run.local == TRUE ) setwd(code.dir)
# 
# if ( run.local == FALSE ) setwd(path)
# 
# source("init_stan_model_MBMA.R")




# RUN SIMULATION ------------------------------

if ( exists("rs") ) rm(rs)

# ~ ********** Beginning of ForEach Loop -----------------------------

# system.time is in seconds
doParallel.seconds = system.time({
  rs = foreach( i = 1:sim.reps, .combine = bind_rows ) %dopar% {
    #for debugging (out file will contain all printed things):
    #for ( i in 1:1000 ) {
    
    # only print info for first sim rep for visual clarity
    if ( i == 1 ) cat("\n\n~~~~~~~~~~~~~~~~ BEGIN SIM REP", i, "~~~~~~~~~~~~~~~~")
    
    # results for just this simulation rep
    if ( exists("rep.res") ) suppressWarnings( rm(rep.res) )
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    if ( verbose == TRUE ) {
      cat("\n\n scen variable:\n")
      print(scen)
      
      cat("\n\n scen.params again:\n")
      print(scen.params)
    }
    
    p = scen.params[ scen.params$scen == scen, names(scen.params) != "scen"]
    
    # calculate TOTAL heterogeneity
    p$V = p$t2a + p$t2w
    p$S = sqrt(p$V)
    
    if ( i == 1 & verbose == TRUE) cat("\n\nDIM AND HEAD OF P (SINGLE ROW OF SCEN.PARAMS):\n")
    if ( i == 1  & verbose == TRUE) print(dim(p)); print(p); print(p$Mu)
    
    # parse methods string
    all.methods = unlist( strsplit( x = p$rep.methods,
                                    split = " ; " ) )
    
    # ~ Simulate Dataset ------------------------------
    # includes unpublished studies
    d = sim_meta_2( Nmax = p$Nmax,
                    true.dist = p$true.dist,
                    true.sei.expr = p$true.sei.expr,
                    Mu = p$Mu,
                    t2a = p$t2a,
                    t2w = p$t2w,
                    m = p$m,
                    
                    hack = p$hack,
                    rho = p$rho,
                    prob.hacked = p$prob.hacked,
                    
                    k.pub.nonaffirm = p$k.pub.nonaffirm,
                    eta = p$eta,
                    gamma = p$gamma,
                    SAS.type = p$SAS.type,
                    
                    muB = p$muB,
                    sig2B = p$sig2B,
                    prob.conf = p$prob.conf,
                    
                    return.only.published = FALSE)
    
    #mean(d$Ci)
    
    
    # sanity
    #View(d %>% select(yi, vi, pval, affirm, Di.across.prob))
    
    # # TEMP    
    # d %>% group_by(affirm) %>%
    #   summarise( mean(Fi), mean(Di.across))
    # 
    # 
    # # using base-R filtering because filter() not working here
    # d[ d$Fi == 1 & d$Di.across == 1, ] %>%
    #   group_by(Ci, affirm) %>%
    #   summarise(n(),
    #             mean(Bi),
    #             mean(mui),
    #             mean(yi))
    
    d$Zi = d$yi / sqrt(d$vi)
    
    
    # ~ For MBMA: Gold-Standard Confounding Adjustment (TRUE muB, sigB) ------------------------------
    
    
    if ( any(d$Ci == 1) ) {
      # this section of code is fine even if d$Ci = 1 always
      d$yi.adj.true = d$yi - d$Ci * p$muB
      d$vi.adj.true = d$vi + d$Ci * p$sig2B
      
      d$tcrit.adj.true = d$tcrit
      d$tcrit.adj.true[ d$Ci == 1 ] = ( d$tcrit[ d$Ci == 1 ] * sqrt(d$vi[ d$Ci == 1 ]) - p$muB ) / sqrt(d$vi.adj.true[ d$Ci == 1 ])
      
      # d %>% group_by(Ci) %>%
      #   summarise( mean(yi),
      #              mean(yi.adj.true),
      #              mean(vi.adj.true),
      #              mean(tcrit),
      #              mean(tcrit.adj.true))
      
      expect_equal( d$affirm,
                    d$yi > d$tcrit * sqrt(d$vi) )
      # confirm that adjusted affirmative threshold is equivalent to the old one
      expect_equal( d$affirm,
                    d$yi.adj.true > d$tcrit.adj.true * sqrt(d$vi.adj.true) )
      
    } else {
      # if Ci = 0 always
      d$yi.adj.true = d$yi
      d$vi.adj.true = d$vi
      d$tcrit.adj.true = d$tcrit
    }
    
    
    
    # ***** For MBMA: Identifiable Confounding Adjustment (ESTIMATED muB, sigB) -----------------
    
    # dataset of only favored AND published results
    # (used in this section)
    #dp = d %>% filter(Fi == 1 & Di.across == 1)  # throwing weird error now?
    dp = d[ d$Fi == 1 & d$Di.across == 1, ]
    
    if ( any(dp$Ci == 1) ) {
      
      # P(A^*_i = a | C^*_i = 1, D^*_i = 1)
      ( P.affirm.pub = mean( dp$affirm[ dp$Ci == 1 ] ) )
      # and P(A^*_i = 0 | C^*_i = 1):
      P.nonaffirm.pub = 1 - P.affirm.pub
      
      # mean bias in published studies
      #  assumed to be correctly specified, so using sample average in underlying studies
      
      # avoid NAs in MhatB.affirm.obs if there are no affirms
      if ( P.affirm.pub > 0 ) {
        MhatB.affirm.obs = mean( dp$Bi[ dp$Ci == 1 & dp$affirm == 1] )
      } else {
        # doesn't matter what value is assigned since will be multiplied
        #  by probability of 0
        MhatB.affirm.obs = 0
      }
      
      if ( P.nonaffirm.pub > 0 ) {
        MhatB.nonaffirm.obs = mean( dp$Bi[ dp$Ci == 1 & dp$affirm == 0] )
      } else {
        # doesn't matter what value is assigned since will be multiplied
        #  by probability of 0
        MhatB.nonaffirm.obs = 0
      }
      
      denom = P.affirm.pub + p$eta * P.nonaffirm.pub
      
      # ~ Sample estimate of muB -------------------
      
      ( MhatB = (1/denom) * ( P.nonaffirm.pub * p$eta * MhatB.nonaffirm.obs +
                                P.affirm.pub * MhatB.affirm.obs ) )
      
      
      # # save: unweighted sample estimate one (just for comparison)
      # # will be way too high
      # mean( dp$Bi[ dp$Ci == 1 ] )
      # 
      # # underlying sample estimate of truth
      # mean( d$Bi[ d$Ci == 1 & d$Di == 1] )
      # 
      # # also should be close to...
      # p$muB
      
      # ~ Sample estimate of sig2B (only used for RTMA, not MBMA) -------------------
      
      # from 2022-7-4 theory
      # estimate *underlying* P(A^*_i = a | C^*_i = 1) from P(A^*_i = a | C^*_i = 1, D^*_i = 1)
      ( Pstar.affirm = P.affirm.pub/denom )
      Pstar.nonaffirm = P.nonaffirm.pub*p$eta/denom
      
      # variance of bias in published studies, assumed known
      # avoid NAs if there are no confounded affirms
      if ( any( dp$Ci == 1 & dp$affirm == 1 ) ) {
        ( shat2B.affirm.obs = var( dp$Bi[ dp$Ci == 1 & dp$affirm == 1 ] ) )
      } else {
        # doesn't matter what value is assigned since will be multiplied
        #  by probability of 0
        shat2B.affirm.obs = 0
      }
      
      if ( any(dp$Ci == 1 & dp$affirm == 0) ) {
        ( shat2B.nonaffirm.obs = var( dp$Bi[ dp$Ci == 1 & dp$affirm == 0 ] ) )
      } else {
        # doesn't matter what value is assigned since will be multiplied
        #  by probability of 0
        shat2B.nonaffirm.obs = 0
      }
      
      # IMPORTANT: at this point, shat2B.affirm.obs could be NA if there is only 1 study
      #  with C^*i = 1 and affirm = 1 (or same for affirm = 0)
      # then shat2B will be NA, and also vi.adj.est
      
      termA = Pstar.affirm*shat2B.affirm.obs + Pstar.nonaffirm*shat2B.nonaffirm.obs
      termB = ( Pstar.affirm * Pstar.nonaffirm ) * ( MhatB.affirm.obs^2 + MhatB.nonaffirm.obs^2 )
      termC = 2 * ( Pstar.affirm * Pstar.nonaffirm ) * ( MhatB.affirm.obs * MhatB.nonaffirm.obs )
      
      ( shat2B = termA + termB - termC )
      
      # # sanity checks
      # # close match
      # Pstar.affirm
      # mean(d$affirm[d$Ci==1])
      # 
      # # sample should be close to underlying because no selection on Bi^*
      # #  though requires quite large k to work
      # shat2B.nonaffirm.obs
      # var(d$Bi[d$Ci == 1 & d$affirm == 0])
      # 
      # # also matches pretty closely
      # shat2B
      # var(d$Bi[d$Ci==1])
      # p$sig2B
      
      
      # ~ ******* Adjusted yi, vi, tcrit in the dataset -------------------
      
      # adjusted yi's using the estimated MhatB
      # using the underlying dataset rather than dp so that we can do sanity checks later
      d$yi.adj.est = d$yi - d$Ci * MhatB
      
      d$vi.adj.est = d$vi + d$Ci * shat2B
      
      d$tcrit.adj.est = d$tcrit
      d$tcrit.adj.est[ d$Ci == 1 ] = ( d$tcrit[ d$Ci == 1 ] *
                                         sqrt(d$vi[ d$Ci == 1 ]) - MhatB ) / sqrt(d$vi.adj.est[ d$Ci == 1 ])
      
      
      # ~ Sanity checks on RTMA reparametrization -------------------
      
      # look for problems calculating the adjusted estimates
      # if ( any(is.na(d$yi.adj.est)) | any(is.na(d$vi.adj.est)) ) {
      # 
      #   cat("\n\n TEMP FLAG yi.adj.est:\n")
      #   print(head(d$yi.adj.est))
      # 
      #   cat("\n\n vi.adj.est:\n")
      #   print(head(d$vi.adj.est))
      # 
      #   cat("\n\n mean(dp$Ci):\n")
      #   print( mean(dp$Ci) )
      # 
      #   cat("\n\n mean(dp$affirm):\n")
      #   print( mean(dp$affirm) )
      # 
      #   cat("\n\n shat2B:\n")
      #   print( shat2B )
      # 
      #   cat("\n\n shat2B.nonaffirm.obs:\n")
      #   print( shat2B.nonaffirm.obs )
      # 
      #   cat("\n\n shat2B.affirm.obs:\n")
      #   print( shat2B.affirm.obs )
      # }
      # 
      # if ( length(d$yi.adj.est) == 0 | length(d$vi.adj.est) == 0 ) {
      # 
      #   cat("\n\n TEMP FLAG yi.adj.est:\n")
      #   print(head(d$yi.adj.est))
      # 
      #   cat("\n\n vi.adj.est:\n")
      #   print(head(d$vi.adj.est))
      # 
      #   cat("\n\n mean(dp$Ci):\n")
      #   print( mean(dp$Ci) )
      # 
      #   cat("\n\n mean(dp$affirm):\n")
      #   print( mean(dp$affirm) )
      # 
      #   cat("\n\n shat2B:\n")
      #   print( shat2B )
      # 
      #   cat("\n\n shat2B.nonaffirm.obs:\n")
      #   print( shat2B.nonaffirm.obs )
      # 
      #   cat("\n\n shat2B.affirm.obs:\n")
      #   print( shat2B.affirm.obs )
      # }
      
      
      expect_equal( d$affirm,
                    d$yi > d$tcrit * sqrt(d$vi) )
      # confirm that adjusted affirmative threshold is equivalent to the old one
      # vi.adj.est could be NA if there is only 1 study with C^*i = 1 and affirm = 1 (or same for affirm = 0)
      if ( all( !is.na(d$vi.adj.est) ) ) {
        expect_equal( d$affirm,
                      d$yi.adj.est > d$tcrit.adj.est * sqrt(d$vi.adj.est) )
      }
      
      
    } else {
      # i.e., if Ci = 0 always
      d$yi.adj.est = d$yi
      d$vi.adj.est = d$vi
      d$tcrit.adj.est = d$tcrit
      # so that returned sanity checks will be okay:
      MhatB = NA
      shat2B = NA
    }
    
    
    
    # ~ Dataset Subsets for Various Methods ------------------------------
    
    # dataset of only favored AND published results
    # overwrite the previous one to get the new variables (e.g., yi.adj.est)
    #dp = d %>% filter(Fi == 1 & Di.across == 1)  # throwing weird error now?
    dp = d[ d$Fi == 1 & d$Di.across == 1, ]
    
    # empirical eta in underlying data
    # useful for weird pub bias methods like carter_censor
    num = conditional_prob( d$Di.across[ d$Fi == 1 ] == 1, d$affirm[ d$Fi == 1 ] == 1 )
    denom = conditional_prob( d$Di.across[ d$Fi == 1 ] == 1, d$affirm[ d$Fi == 1 ] == 0 )
    eta_emp = num/denom
    
    # published nonaffirmatives only
    dpn = dp[ dp$affirm == FALSE, ]
    
    # published affirmatives only
    dpa = dp[ dp$affirm == TRUE, ]
    
    if ( verbose == TRUE ) {
      if ( i == 1 ) cat("\n\nHEAD OF DP:\n")
      if ( i == 1 ) print(head(dp))
    }
    
    # ~ Start Values ------------------------------
    
    # initialize rep.res st run_method_safe and other standalone estimation fns
    #  will correctly recognize it as having 0 rows
    rep.res = data.frame()
    
    Mhat.start = p$Mu
    Shat.start = p$S
    
    # in case we're not doing rtma or it fails
    Mhat.MaxLP = NA
    Shat.MaxLP = NA
    
    Mhat.MAP = NA
    Shat.MAP = NA
    
    
    # ~ Existing Methods ------------------------------
    
    # ~~ Naive Meta-Analysis (All PUBLISHED Draws)
    
    if ( "naive" %in% all.methods ) {
      rep.res = run_method_safe(method.label = c("naive"),
                                method.fn = function() {
                                  mod = rma( yi = dp$yi,
                                             vi = dp$vi,
                                             method = "REML",
                                             knha = TRUE )
                                  
                                  report_meta(mod, .mod.type = "rma")
                                },
                                .rep.res = rep.res )
    }
    
    rep.res
    
    
    # ~~ 2PSM (All Published Draws)
    
    if ( "2psm" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("2psm"),
                                method.fn = function() {
                                  mod = weightfunct( effect = dp$yi,
                                                     v = dp$vi,
                                                     steps = c(0.025, 1),
                                                     table = TRUE ) 
                                  
                                  H = mod[[2]]$hessian
                                  ses = sqrt( diag( solve(H) ) )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = mod[[2]]$par[2],
                                                            MLo = mod[[2]]$par[2] - qnorm(.975) * ses[2],
                                                            MHi = mod[[2]]$par[2] + qnorm(.975) * ses[2],
                                                            
                                                            Shat = sqrt( mod[[2]]$par[1] ),
                                                            # truncate lower limit at 0
                                                            SLo = sqrt( max( 0, mod[[2]]$par[1] - qnorm(.975) * ses[1] ) ),
                                                            SHi = sqrt( mod[[2]]$par[1] + qnorm(.975) * ses[1] ),
                                                            
                                                            # estimated OVERALL selection ratio (eta*gamma)
                                                            # index [3] would change for PSM with more cut points
                                                            EtaGammaHat = 1/mod$output_adj$par[3]
                                  ) )
                                  
                                  
                                },
                                .rep.res = rep.res )
      
    }
    
    
    rep.res
    
    
    # then implement the other easy DGP things and make sure they all run locally
    # you got this! oh yeah! :)
    
    
    # ~~ 2PSM (All Published Draws)
    
    if ( "beta-sm" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("beta-sm"),
                                method.fn = function() {
                                  
                                  # must start with naive fit for selmodel
                                  naive = rma( yi = dp$yi,
                                               vi = dp$vi,
                                               method = "REML",
                                               knha = TRUE )
                                  
                                  mod = selmodel(x = naive,
                                                 type = "beta",
                                                 alternative = "two.sided")
                                  
                                  # not returning any info about selection parameters because they
                                  #  don't have same interpretation as eta
                                  report_meta(mod, .mod.type = "rma")
                                  
                                  # IMPORTANT: 
                                  # beta-sm is prone to the warning: "error when trying to invert Hessian"
                                  #  in which case there can be a point estimate with no SEs or tau
                                  #  in this case, selmodel will still return something,
                                  #  but report_meta will throw the missing value error from tau_CI
                                  #  this behavior is what we want, but does mean that the overall.error
                                  #  for beta-sm is not very informative
  
                                },
                                .rep.res = rep.res )
      
    }
    
    
    rep.res
    
    
    # ~ New Methods ------------------------------
    
    # ~~ ****** MBMA ------------------------------
    
    
    # using the identifiable, reweighting-based sample estimate of muB ("MhatB")
    if ( "mbma-MhatB" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("mbma-MhatB"),
                                method.fn = function() {
                                  
                                  # from inside PublicationBias::corrected_meta;
                                  #  only change is that we want affirm indicator to be that of the *confounded* estimates, not the adjusted ones
                                  # weight for model
                                  weights = rep( 1, length(dp$yi.adj.est) )
                                  
                                  
                                  # set weights based on SAS type
                                  # weight based on the affirm indicator of the *confounded* estimates
                                  weights[ dp$affirm == FALSE ] = p$eta  # default
                                  if ( p$SAS.type == "carter_censor" ) {
                                    # this SAS mechanism doesn't use p$eta, so need to get empirical estimate
                                    #**mention in paper
                                    weights[ dp$affirm == FALSE ] = eta_emp
                                  }
                                  
                                  # initialize a dumb (unclustered and uncorrected) version of tau^2
                                  # which is only used for constructing weights
                                  meta.re = rma.uni( yi = dp$yi.adj.est,
                                                     vi = dp$vi)
                                  t2hat.naive = meta.re$tau2  # could subtract off the sig2B here, but would also need to account for some studies' being unconfounded
                                  
                                  
                                  # fit weighted robust model
                                  meta.robu = robu( yi.adj.est ~ 1,
                                                    studynum = 1:nrow(dp),
                                                    data = dp,
                                                    userweights = weights / (vi + t2hat.naive),
                                                    var.eff.size = vi,
                                                    small = TRUE )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = as.numeric(meta.robu$b.r),
                                                            MLo = meta.robu$reg_table$CI.L,
                                                            MHi = meta.robu$reg_table$CI.U,
                                                            
                                                            Shat = NA,
                                                            SLo = NA,
                                                            SHi = NA,
                                                            
                                                            EtaEmpCarter = eta_emp ) ) 
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done mbma-MhatB if applicable")
      
    }
    
    rep.res
    
    
    
    # NEW: same as mbma-MhatB, but uses eta*gamma as selection ratio (hard coded for now)
    if ( "mbma-MhatB-gamma" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("mbma-MhatB-gamma"),
                                method.fn = function() {
                                  
                                  # from inside PublicationBias::corrected_meta;
                                  #  only change is that we want affirm indicator to be that of the *confounded* estimates, not the adjusted ones
                                  # weight for model
                                  weights = rep( 1, length(dp$yi.adj.est) )
                                  # weight based on the affirm indicator of the *confounded* estimates
                                  # **note this uses the empirical gamma from underlying data to accommodate hacking methods in which 
                                  #  we don't specify the population gamma
                                  EtaGammaAssumed = p$eta * dp$gamma[1]
                                  weights[ dp$affirm == FALSE ] = EtaGammaAssumed
                                  
                                  # initialize a dumb (unclustered and uncorrected) version of tau^2
                                  # which is only used for constructing weights
                                  meta.re = rma.uni( yi = dp$yi.adj.est,
                                                     vi = dp$vi)
                                  t2hat.naive = meta.re$tau2  # could subtract off the sig2B here, but would also need to account for some studies' being unconfounded
                                  
                                  
                                  # fit weighted robust model
                                  meta.robu = robu( yi.adj.est ~ 1,
                                                    studynum = 1:nrow(dp),
                                                    data = dp,
                                                    userweights = weights / (vi + t2hat.naive),
                                                    var.eff.size = vi,
                                                    small = TRUE )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = as.numeric(meta.robu$b.r),
                                                            MLo = meta.robu$reg_table$CI.L,
                                                            MHi = meta.robu$reg_table$CI.U,
                                                            
                                                            Shat = NA,
                                                            SLo = NA,
                                                            SHi = NA,
                                                            
                                                            EtaGammaAssumed = EtaGammaAssumed) ) 
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done mbma-MhatB-gamma if applicable")
      
    }
    
    rep.res
    
    
    # Benchmark: using MhatB but with TRUE tau^2
    if ( "mbma-MhatB-true-t2" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("mbma-Mhat-true-t2"),
                                method.fn = function() {
                                  
                                  # from inside PublicationBias::corrected_meta;
                                  #  only change is that we want affirm indicator to be that of the *confounded* estimates, not the adjusted ones
                                  # weight for model
                                  weights = rep( 1, length(dp$yi.adj.est) )
                                  # weight based on the affirm indicator of the *confounded* estimates
                                  weights[ dp$affirm == FALSE ] = p$eta 
                                  
                                  # fit weighted robust model
                                  meta.robu = robu( yi.adj.est ~ 1,
                                                    studynum = 1:nrow(dp),
                                                    data = dp,
                                                    # here uses true t2:
                                                    userweights = weights / (vi + p$t2a + p$t2w),
                                                    var.eff.size = vi,
                                                    small = TRUE )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = as.numeric(meta.robu$b.r),
                                                            MLo = meta.robu$reg_table$CI.L,
                                                            MHi = meta.robu$reg_table$CI.U,
                                                            
                                                            Shat = NA,
                                                            SLo = NA,
                                                            SHi = NA ) ) 
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done mbma-Mhat-true-t2 if applicable")
      
    }
    
    
    # using the true muB
    if ( "mbma-muB" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("mbma-muB"),
                                method.fn = function() {
                                  
                                  # from inside PublicationBias::corrected_meta;
                                  #  only change is that we want affirm indicator to be that of the *confounded* estimates, not the adjusted ones
                                  # weight for model
                                  weights = rep( 1, length(dp$yi.adj.true) )
                                  # weight based on the affirm indicator of the *confounded* estimates
                                  weights[ dp$affirm == FALSE ] = p$eta
                                  
                                  # initialize a dumb (unclustered and uncorrected) version of tau^2
                                  # which is only used for constructing weights
                                  meta.re = rma.uni( yi = dp$yi.adj.true,
                                                     vi = dp$vi)
                                  t2hat.naive = meta.re$tau2  # could subtract off the sig2B here, but would also need to account for some studies' being unconfounded
                                  
                                  
                                  # fit weighted robust model
                                  meta.robu = robu( yi.adj.true ~ 1,
                                                    studynum = 1:nrow(dp),
                                                    data = dp,
                                                    userweights = weights / (vi + t2hat.naive),
                                                    var.eff.size = vi,
                                                    small = TRUE )
                                  
                                  # follow the same return structure as report_meta
                                  list( stats = data.frame( Mhat = as.numeric(meta.robu$b.r),
                                                            MLo = meta.robu$reg_table$CI.L,
                                                            MHi = meta.robu$reg_table$CI.U,
                                                            
                                                            Shat = NA,
                                                            SLo = NA,
                                                            SHi = NA ) ) 
                                },
                                .rep.res = rep.res )
      
      cat("\n doParallel flag: Done mbma-muB if applicable")
      
    }
    
    
    
    
    #srr()
    
    # ~~ ********* RTMA WITH CONFOUNDING ADJUSTMENT ------------------------------
    
    # RTMA with true bias parameters
    if ( "rtma-adj-muB" %in% all.methods ) {
      # # temp for refreshing stan code
      # path = "/home/groups/manishad/MBMA"
      #setwd(path)
      # source("helper_MBMA.R")
      # source("init_stan_model_MBMA.R")
      
      # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("rtma-adj-muB-pmean",
                                                 "rtma-adj-muB-pmed",
                                                 "rtma-adj-muB-max-lp-iterate"),
                                # note that we're now passing the confounding-adjusted estimates, variances,
                                #  and critical values
                                method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi.adj.true,
                                                                                   .sei = sqrt(dpn$vi.adj.true),
                                                                                   .tcrit = dpn$tcrit.adj.true,
                                                                                   .Mu.start = Mhat.start,
                                                                                   # can't handle start value of 0:
                                                                                   .Tt.start = max(0.01, Shat.start),
                                                                                   .stan.adapt_delta = p$stan.adapt_delta,
                                                                                   .stan.maxtreedepth = p$stan.maxtreedepth), .rep.res = rep.res )
      
      
      Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "rtma-adj-muB-max-lp-iterate" ]
      Shat.MaxLP = rep.res$Shat[ rep.res$method == "rtma-adj-muB-max-lp-iterate" ]
      
      cat("\n doParallel flag: Done rtma-adj-muB if applicable")
    }
    
    #srr()
    
    # ~~ Change Starting Values -----
    if ( !is.na(Mhat.MaxLP) ) Mhat.start = Mhat.MaxLP
    if ( !is.na(Shat.MaxLP) ) Shat.start = Shat.MaxLP 
    
    
    
    # RTMA with identifiable (estimated) bias parameters
    
    if ( "rtma-adj-MhatB" %in% all.methods ) {
      # # temp for refreshing stan code
      # path = "/home/groups/manishad/MBMA"
      #setwd(path)
      # source("helper_MBMA.R")
      # source("init_stan_model_MBMA.R")
      
      # this one has two labels in method arg because a single call to estimate_jeffreys_mcmc
      #  returns 2 lines of output, one for posterior mean and one for posterior median
      # order of labels in method arg needs to match return structure of estimate_jeffreys_mcmc
      rep.res = run_method_safe(method.label = c("rtma-adj-MhatB-pmean",
                                                 "rtma-adj-MhatB-pmed",
                                                 "rtma-adj-MhatB-max-lp-iterate"),
                                # note that we're now passing the confounding-adjusted estimates, variances,
                                #  and critical values
                                method.fn = function() estimate_jeffreys_mcmc_RTMA(.yi = dpn$yi.adj.est,
                                                                                   .sei = sqrt(dpn$vi.adj.est),
                                                                                   .tcrit = dpn$tcrit.adj.est,
                                                                                   .Mu.start = Mhat.start,
                                                                                   # can't handle start value of 0:
                                                                                   .Tt.start = max(0.01, Shat.start),
                                                                                   .stan.adapt_delta = p$stan.adapt_delta,
                                                                                   .stan.maxtreedepth = p$stan.maxtreedepth), .rep.res = rep.res )
      
      
      Mhat.MaxLP = rep.res$Mhat[ rep.res$method == "rtma-adj-MhatB-max-lp-iterate" ]
      Shat.MaxLP = rep.res$Shat[ rep.res$method == "rtma-adj-MhatB-max-lp-iterate" ]
      
      cat("\n doParallel flag: Done rtma-adj-MhatB if applicable")
    }
    
    #srr()
    
    # ~~ Change Starting Values -----
    if ( !is.na(Mhat.MaxLP) ) Mhat.start = Mhat.MaxLP
    if ( !is.na(Shat.MaxLP) ) Shat.start = Shat.MaxLP 
    
    
    
    
    # ~~ ****** MAP (SD param) ------------------------------
    
    # THIS IS STILL USING THE TRUE ADJUSTMENT
    if ( "jeffreys-adj-sd" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("jeffreys-adj-sd"),
                                method.fn = function() estimate_jeffreys_RTMA(yi = dpn$yi.adj.true,
                                                                              sei = sqrt(dpn$vi.adj.true), 
                                                                              par2is = "Tt",
                                                                              tcrit = dpn$tcrit.adj.true, 
                                                                              Mu.start = Mhat.start,
                                                                              par2.start = Shat.start,
                                                                              usePrior = TRUE,
                                                                              get.CIs = p$get.CIs,
                                                                              CI.method = "wald",
                                                                              run.optimx = p$run.optimx
                                ),
                                .rep.res = rep.res )
      
      Mhat.MAP = rep.res$Mhat[ rep.res$method == "jeffreys-adj-sd" ]
      Shat.MAP = rep.res$Shat[ rep.res$method == "jeffreys-adj-sd" ]
    }
    
    # ~~ ******** MAON WITH CONFOUNDING ADJUSTMENT --------------------------------------
    
    # using the true muB
    if ( "maon-adj-muB" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("maon-adj-muB"),
                                method.fn = function() {
                                  mod = robu( yi.adj.true ~ 1, 
                                              data = dpn, 
                                              studynum = 1:nrow(dpn),
                                              var.eff.size = vi,  # using original variance
                                              small = TRUE )
                                  
                                  report_meta(mod, .mod.type = "robu")
                                },
                                .rep.res = rep.res )
      
    }
    
    # using the sample estimate of muB
    if ( "maon-adj-MhatB" %in% all.methods ) {
      
      rep.res = run_method_safe(method.label = c("maon-adj-MhatB"),
                                method.fn = function() {
                                  mod = robu( yi.adj.est ~ 1, 
                                              data = dpn, 
                                              studynum = 1:nrow(dpn),
                                              var.eff.size = vi,  # using original variance
                                              small = TRUE )
                                  
                                  report_meta(mod, .mod.type = "robu")
                                },
                                .rep.res = rep.res )
      
    }
    
    if ( verbose == TRUE ) srr()
    
    # ~ Add Scen Params and Sanity Checks --------------------------------------
    
    # add in scenario parameters
    # do NOT use rbind here; bind_cols accommodates possibility that some methods' rep.res
    #  have more columns than others
    rep.res = p %>% bind_cols( rep.res )
    
    # add more info
    rep.res = rep.res %>% add_column( rep.name = i, .before = 1 )
    rep.res = rep.res %>% add_column( scen.name = scen, .before = 1 )
    rep.res = rep.res %>% add_column( job.name = jobname, .before = 1 )
    
    
    cat("\ndoParallel flag: Before adding sanity checks to rep.res")
    
    
    
    # add info about simulated datasets
    rep.res = rep.res %>% add_column(   sancheck.dp.k = nrow(dp),
                                        sancheck.dp.k.affirm = sum(dp$affirm == TRUE),
                                        sancheck.dp.k.nonaffirm = sum(dp$affirm == FALSE),
                                        
                                        sancheck.dp.k.conf = sum( dp$Ci == 1 ),
                                        sancheck.dp.k.affirm.conf = ifelse( sum(dp$affirm) > 0, sum( dp$Ci[ dp$affirm == 1 ] ), NA ),
                                        sancheck.dp.k.nonaffirm.conf = ifelse( sum(dp$affirm) < 1, sum( dp$Ci[ dp$affirm == 0 ] ), NA ),
                                        
                                        # sanity checks for MhatB
                                        # E[Bi^* | Ci^* = 1], the target for MhatB:
                                        sancheck.MhatB = MhatB,
                                        # note: Fi = 1 here is FAVORING indicator, so this is still the mean Bi among underlying (pre-SAS) estimates
                                        # this is an underlying SAMPLE estimate of the truth; should approximately agree with MhatB and muB
                                        sancheck.EBsti = ifelse( mean(d$Ci == 1 & d$Fi == 1) > 0,
                                                                 mean( d$Bi[ d$Ci == 1 & d$Fi == 1] ),
                                                                 NA ),
                                        
                                        sancheck.shat2B = shat2B )
    
    
    # MORE SANITY CHECKS THAT MAY FAIL IF, E.G., CIi=0 ALWAYS:
    # # add info about simulated datasets
    # # "ustudies"/"udraws" refers to underlying studies/draws prior to hacking or publication bias
    # sancheck.prob.published.is.confounded = mean( dp$Ci == 1 )
    # sancheck.prob.published.affirm.is.confounded = mean( dp$Ci[ dp$affirm == 1 ] == 1 )
    # sancheck.prob.published.nonaffirm.is.confounded = mean( dp$Ci[ dp$affirm == 0 ] == 0 )
    # 
    # 
    # ( sancheck.prob.ustudies.published =  mean( d.first$study %in% unique(dp$study) ) )
    # expect_equal( sancheck.prob.ustudies.published, nrow(dp)/nrow(d.first) )
    # # this one should always be 100% unless there's also publication bias:
    # ( sancheck.prob.unhacked.ustudies.published =  mean( d.first$study[ d.first$hack == "no" ] %in% unique( dp$study[ dp$hack == "no" ] ) ) )
    # # under affirm hacking, will be <100%:
    # ( sancheck.prob.hacked.ustudies.published =  mean( d.first$study[ d.first$hack != "no" ] %in% unique( dp$study[ dp$hack != "no" ] ) ) )
    # 
    # # might NOT be 100% if you're generating multiple draws per unhacked studies but favoring, e.g., a random one:
    # ( sancheck.prob.unhacked.udraws.published =  mean( d$study.draw[ d$hack == "no" ] %in% unique( dp$study.draw[ dp$hack == "no" ] ) ) )
    # ( sancheck.prob.hacked.udraws.published =  mean( d$study.draw[ d$hack != "no" ] %in% unique( dp$study.draw[ dp$hack != "no" ] ) ) )
    # 
    # 
    # 
    # #*this one is especially important: under worst-case hacking, it's analogous to prop.retained in
    # #  TNE since it's the proportion of the underlying distribution that's nonaffirmative
    # ( sancheck.prob.unhacked.udraws.nonaffirm =  mean( d$affirm[ d$hack == "no" ] == FALSE ) )
    # # a benchmark for average power:
    # ( sancheck.prob.unhacked.udraws.affirm =  mean( d$affirm[ d$hack == "no" ] ) )
    # ( sancheck.prob.hacked.udraws.nonaffirm =  mean( d$affirm[ d$hack != "no" ] == FALSE ) )
    # ( sancheck.prob.hacked.udraws.affirm =  mean( d$affirm[ d$hack != "no" ] ) )
    # 
    # # probability that a published, nonaffirmative draw is from a hacked study
    # # under worst-case hacking, should be 0
    # ( sancheck.prob.published.nonaffirm.is.hacked = mean( dp$hack[ dp$affirm == 0 ] != "no" ) )
    # # this will be >0
    # ( sancheck.prob.published.affirm.is.hacked = mean( dp$hack[ dp$affirm == 1 ] != "no" ) )
    # 
    # rep.res = rep.res %>% add_column(   sancheck.dp.k = nrow(dp),
    #                                     sancheck.dp.k.affirm = sum(dp$affirm == TRUE),
    #                                     sancheck.dp.k.nonaffirm = sum(dp$affirm == FALSE),
    # 
    #                                     sancheck.prob.published.is.confounded = sancheck.prob.published.is.confounded,
    #                                     sancheck.prob.published.affirm.is.confounded = sancheck.prob.published.affirm.is.confounded,
    #                                     sancheck.prob.published.nonaffirm.is.confounded = sancheck.prob.published.nonaffirm.is.confounded,
    # 
    #                                     sancheck.dp.k.affirm.unhacked = sum(dp$affirm == TRUE & dp$hack == "no"),
    #                                     sancheck.dp.k.affirm.hacked = sum(dp$affirm == TRUE & dp$hack != "no"),
    #                                     sancheck.dp.k.nonaffirm.unhacked = sum(dp$affirm == FALSE & dp$hack == "no"),
    #                                     sancheck.dp.k.nonaffirm.hacked = sum(dp$affirm == FALSE & dp$hack != "no"),
    # 
    #                                     # means draws per HACKED, published study
    #                                     sancheck.dp.meanN.hacked = mean( dp$N[dp$hack != "no"] ),
    #                                     sancheck.dp.q90N.hacked = quantile( dp$N[dp$hack != "no"], 0.90 ),
    # 
    #                                     # average yi's of published draws from each study type
    #                                     sancheck.mean.yi.unhacked.pub.study = mean( dp$yi[ dp$hack == "no"] ),
    #                                     sancheck.mean.yi.hacked.pub.study = mean( dp$yi[ dp$hack != "no"] ),
    # 
    # 
    #                                     sancheck.mean.mui.unhacked.pub.nonaffirm = mean( dp$mui[ dp$hack == "no" & dp$affirm == FALSE ] ),
    #                                     sancheck.mean.yi.unhacked.pub.nonaffirm = mean( dp$yi[ dp$hack == "no" & dp$affirm == FALSE ] ),
    #                                     sancheck.mean.yi.unhacked.pub.affirm = mean( dp$yi[ dp$hack == "no" & dp$affirm == TRUE ] ),
    # 
    #                                     sancheck.mean.yi.hacked.pub.nonaffirm = mean( dp$yi[ dp$hack != "no" & dp$affirm == FALSE ] ),
    #                                     sancheck.mean.yi.hacked.pub.affirm = mean( dp$yi[ dp$hack != "no" & dp$affirm == TRUE ] ),
    # 
    #                                     # average Zi's
    #                                     sancheck.mean.Zi.unhacked.pub.study = mean( dp$Zi[ dp$hack == "no"] ),
    #                                     sancheck.mean.Zi.hacked.pub.study = mean( dp$Zi[ dp$hack != "no"] ),
    # 
    #                                     sancheck.mean.Zi.unhacked.pub.nonaffirm = mean( dp$Zi[ dp$hack == "no" & dp$affirm == FALSE ] ),
    #                                     sancheck.mean.Zi.unhacked.pub.affirm = mean( dp$Zi[ dp$hack == "no" & dp$affirm == TRUE ] ),
    # 
    #                                     sancheck.mean.Zi.hacked.pub.nonaffirm = mean( dp$Zi[ dp$hack != "no" & dp$affirm == FALSE ] ),
    #                                     sancheck.mean.Zi.hacked.pub.affirm = mean( dp$Zi[ dp$hack != "no" & dp$affirm == TRUE ] ),
    # 
    # 
    #                                     sancheck.prob.ustudies.published = sancheck.prob.ustudies.published,
    #                                     sancheck.prob.unhacked.ustudies.published = sancheck.prob.unhacked.ustudies.published,
    #                                     sancheck.prob.hacked.ustudies.published = sancheck.prob.hacked.ustudies.published,
    # 
    #                                     sancheck.prob.unhacked.udraws.published = sancheck.prob.unhacked.udraws.published,
    #                                     sancheck.prob.hacked.udraws.published = sancheck.prob.hacked.udraws.published,
    # 
    #                                     sancheck.prob.unhacked.udraws.nonaffirm = sancheck.prob.unhacked.udraws.nonaffirm,
    #                                     sancheck.prob.unhacked.udraws.affirm = sancheck.prob.unhacked.udraws.affirm,
    #                                     sancheck.prob.hacked.udraws.nonaffirm = sancheck.prob.hacked.udraws.nonaffirm,
    #                                     sancheck.prob.hacked.udraws.affirm = sancheck.prob.hacked.udraws.affirm,
    # 
    #                                     sancheck.prob.published.nonaffirm.is.hacked = sancheck.prob.published.nonaffirm.is.hacked,
    # 
    #                                     # sanity checks for MhatB
    #                                     # E[Bi^* | Ci^* = 1], the target for MhatB:
    #                                     sancheck.MhatB = MhatB,
    #                                     #@note: Di = 1 here is FAVORING indicator, so this is still the mean Bi among underlying (pre-SAS) estimates
    #                                     # this is an underlying SAMPLE estimate of the truth; should approximately agree with MhatB and muB
    #                                     sancheck.EBsti = mean( d$Bi[ d$Ci == 1 & d$Di == 1] ),
    # 
    #                                     sancheck.shat2B = shat2B )
    
    rep.res
    
  }  ### end foreach loop
  
} )[3]  # end system.time





# dim(rs)
# # quick look
# rs %>% mutate(MhatWidth = MHi - MLo,
#               MhatCover = as.numeric( MHi > Mu & MLo < Mu ) ) %>%
#   dplyr::select(method, Mhat, MhatWidth, MhatCover,EtaGammaHat,EtaGammaAssumed,
#                 sancheck.MhatB, sancheck.EBsti) %>%
# 
#   group_by(method) %>%
#   summarise_if(is.numeric, function(x) round( meanNA(x), 2 ) )
# 
# any(is.na(rs$sancheck.MhatB))



# ~~ End of ForEach Loop ----------------
# estimated time for 1 simulation rep
# use NAs for additional methods so that the SUM of the rep times will be the
#  total computational time
nMethods = length( unique(rs$method) )

print(nMethods)
print(unique(rs$method))
table(rs$method)

print(nrow(rs))

rs$doParallel.seconds = doParallel.seconds


rs$rep.seconds = doParallel.seconds/sim.reps
rs$rep.seconds[ rs$method != unique(rs$method)[1] ] = NA

#rs$rep.seconds = rep( c( doParallel.seconds / sim.reps,
#                         rep( NA, nMethods - 1 ) ), sim.reps )

expect_equal( as.numeric( sum(rs$rep.seconds, na.rm = TRUE) ),
              as.numeric(doParallel.seconds) )



# ~ QUICK RESULTS SUMMARY ---------------------------------------------------

if ( run.local == TRUE ) {
  #if ( TRUE ) {  
  # quick look locally
  # rs %>% mutate_if(is.numeric, function(x) round(x,2) )
  
  agg = rs %>% group_by(method) %>%
    summarise( PropMhatNA = mean(is.na(Mhat)),
               PropCI.NA = mean(is.na(MLo)),
               
               Mhat = meanNA(Mhat),
               MhatMSE = meanNA( (Mhat - Mu)^2 ),
               MhatBias = meanNA(Mhat - Mu),
               MhatEmpSE = sd( Mhat, na.rm = TRUE ),
               #ShatMn = meanNA(Shat),
               
               MhatCover = meanNA(MLo < Mu & MHi > Mu),
               MhatWidth = meanNA( MHi - MLo ),
               MLo = meanNA(MLo),
               MHi = meanNA(MHi),
               
               EtaGammaAssumed = meanNA(EtaGammaAssumed),
               EtaGammaHat = meanNA(EtaGammaHat),
               EtaEmpCarter = meanNA(EtaEmpCarter) )
  
  # round
  agg = as.data.frame( agg %>% mutate_if( is.numeric,
                                          function(x) round(x,2) ) )
  
  agg
  
  
  # # scenario diagnostics for scenario
  # keepers = namesWith("sancheck.", rs)
  # agg.checks = rs %>% summarise_at( keepers,
  #                                   function(x) round( mean(x), 2) )
  #t(agg.checks)
  
}




# ~ WRITE LONG RESULTS ------------------------------
if ( run.local == FALSE ) {
  setwd("/home/groups/manishad/MBMA/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
}