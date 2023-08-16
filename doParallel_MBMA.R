
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
run.local = FALSE

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
           "stringr",
           "weightr")

if ( run.local == TRUE | interactive.cluster.run == TRUE ) toLoad = c(toLoad, "here")


# SET UP FOR CLUSTER OR LOCAL RUN ------------------------------

# ~~ Cluster Run ----------------------------------------
if (run.local == FALSE) {
  
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
  
  # load command line arguments (must be after loading packages)
  args = commandArgs(trailingOnly = TRUE)
  cat("\n\n args received from sbatch file:", args)
  jobname = args[1]
  #scen = args[2]
  

  # alt: avoid stringr::str_split becuase it gives weird error on cluster
  scens.to.run = unlist(strsplit(paste0(args[2], ","), ",")) # this will be a string like "1,2,3,4"
  
  
  cat("\n\n Parsed scenarios:", scens.to.run)
  
  
  
  # helper code
  path = "/home/groups/manishad/MBMA"
  setwd(path)
  source("helper_MBMA.R")
  source("analyze_sims_helper_MBMA.R")  # for make_agg_data
  
  # ~~ Cluster Run ----------------------------------------
  
  if ( interactive.cluster.run == FALSE ) {
    # get scen parameters (made by genSbatch.R)
    setwd(path)
    scen.params = read.csv( "scen_params.csv" )
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
  
  # simulation reps per scenario
  if ( interactive.cluster.run == FALSE ) sim.reps = 1000  # expect about 1 hr per 1000 sim.reps if running 8 scens per doParallel
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
  
  # exactly as on cluster

  ### 2022-7-23 - debugging set ###
  scen.params = tidyr::expand_grid(
    
    rep.methods = "naive ; mbma-MhatB",
    
    
    # args from sim_meta_2
    Nmax = 5,  # later code will set this to 1 if prob.hacked = 0
    true.dist = c("expo", "norm"),
    true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)",  # original setting close to empirical distribution
                      "0.02 + rexp(n = 1, rate = 1)"),  # larger SEs overall
    Mu = c(0, 0.5),
    t2a = c(0, 0.25^2, 0.5^2),
    t2w = c(0),
    m = 50,
    
    # SWS args
    # remember: method affirm will only work if prob.hacked < 1 else will never have nonaffirms
    hack = c("affirm2", "favor-lowest-p", "favor-gamma-ratio"),  
    rho = c(0),
    prob.hacked = c(1, 0),
    
    # SAS args
    k.pub.nonaffirm = c(5, 10, 15, 30, 50),
    eta = c(1, 5, 10),
    gamma = 2,  # only used for method favor-gamma-ratio
    SAS.type = c("2psm", "carter"), # "2psm" (original) or "carter"
    
    
    
    # confounding parameters
    # NOT using log scale here b/c underlying data are continuous
    muB = c(0.1, 0.25, 0.5),
    sig2B = 0.5,
    prob.conf = c(0.5),
    
    # Stan control args - only relevant if running RTMA - remove these args?
    stan.maxtreedepth = 25,
    stan.adapt_delta = 0.995,
    
    get.CIs = TRUE,
    run.optimx = FALSE )
  
  # if there are multiple hacking types, remove redundant combos
  # i.e., only need 1 hack type with p.hacked = 0
  first.hack.type = unique(scen.params$hack)[1]
  scen.params = scen.params %>% filter( prob.hacked > 0 | (prob.hacked == 0 & hack == first.hack.type) )
  
  # no need to make extra draws is prob.hacked = 0
  scen.params$Nmax[ scen.params$prob.hacked == 0 ] = 1
  
  scen.params %>% group_by(prob.hacked, hack) %>%
    summarise( mean(Nmax) )
  
  # add scen numbers
  start.at = 1
  scen.params = scen.params %>% add_column( scen = start.at : ( nrow(scen.params) + (start.at - 1) ),
                                            .before = 1 )
  
  scen.params$scen = 1:nrow(scen.params)
  
  
  sim.reps = 10  # reps to run in this iterate
  
  # set the number of local cores
  registerDoParallel(cores=8)

  scens.to.run=7144

  # make sure we made enough scens to actually run them
  expect_equal( length(scens.to.run) <= length(scen.params$scen), TRUE )
  
  # just to avoid errors in doParallel script below
  jobname = "job_1"
  i = 1
}



# RUN SIMULATION ------------------------------

if ( exists("rs") ) rm(rs)

# ~ ********** Beginning of ForEach Loop -----------------------------

# system.time is in seconds
doParallel.seconds = system.time({
  
  for ( scen in scens.to.run ) {
    
    cat("\n\n~~~~~~~~~~~~~~~~ BEGIN SCEN", scen, "~~~~~~~~~~~~~~~~")
    
    new_rs = foreach( i = 1:sim.reps, .combine = bind_rows ) %dopar% {
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
    
      
      d$Zi = d$yi / sqrt(d$vi)
      
      
      # ~ For MBMA: Gold-Standard Confounding Adjustment (TRUE muB, sigB) ------------------------------
      
      
      if ( any(d$Ci == 1) ) {
        # this section of code is fine even if d$Ci = 1 always
        d$yi.adj.true = d$yi - d$Ci * p$muB
        d$vi.adj.true = d$vi + d$Ci * p$sig2B
        
        d$tcrit.adj.true = d$tcrit
        d$tcrit.adj.true[ d$Ci == 1 ] = ( d$tcrit[ d$Ci == 1 ] * sqrt(d$vi[ d$Ci == 1 ]) - p$muB ) / sqrt(d$vi.adj.true[ d$Ci == 1 ])
        
        
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
      
      # dataset of only favored AND published resultsdoParallel
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
  
        
        
        # ~ ******* Adjusted yi, vi, tcrit in the dataset -------------------
        
        # adjusted yi's using the estimated MhatB
        # using the underlying dataset rather than dp so that we can do sanity checks later
        d$yi.adj.est = d$yi - d$Ci * MhatB
        
        d$vi.adj.est = d$vi + d$Ci * shat2B
        
        d$tcrit.adj.est = d$tcrit
        d$tcrit.adj.est[ d$Ci == 1 ] = ( d$tcrit[ d$Ci == 1 ] *
                                           sqrt(d$vi[ d$Ci == 1 ]) - MhatB ) / sqrt(d$vi.adj.est[ d$Ci == 1 ])
        

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
      
      ### @temp debugging - figure out which of san checks throws "arg of length zero"
      cat( nrow(dp) )
      cat( sum(dp$affirm == TRUE) )
      cat( sum(dp$affirm == FALSE) )
      
      cat( sum( dp$Ci == 1 ) )
      cat( ifelse( sum(dp$affirm) > 0, sum( dp$Ci[ dp$affirm == 1 ] ), NA ) )
      cat( ifelse( sum(dp$affirm) < 1, sum( dp$Ci[ dp$affirm == 0 ] ), NA ) )
      
      # sanity checks for MhatB
      # E[Bi^* | Ci^* = 1], the target for MhatB:
      cat( MhatB )
      # note: Fi = 1 here is FAVORING indicator, so this is still the mean Bi among underlying (pre-SAS) estimates
      # this is an underlying SAMPLE estimate of the truth; should approximately agree with MhatB and muB
      cat( ifelse( mean(d$Ci == 1 & d$Fi == 1) > 0,
                               mean( d$Bi[ d$Ci == 1 & d$Fi == 1] ),
                               NA ) )
      
      cat( shat2B )
      ### @end temp debugging
      
      
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
      
      
     
      rep.res
      
    }  ### end foreach loop for ONE scen in scens.to.run
    
    
    # aggregate results across multiple scens
    
    
    # # flags in case script fails
    cat( paste("\n\ndoParallel flag. nrow(new_rs):" ) ); nrow(new_rs)
    cat( paste("\n\ndoParallel flag. nrow(rs):" ) ); if (exists("rs")) nrow(rs)
    cat( paste("\n\ndoParallel flag. scens.to.run[1]:" ) ); scens.to.run[1]
    
    if ( scen == scens.to.run[1] ) {
      cat( paste("\n\ndoParallel flag: IF"))
      # MUST be superassignment or else rs won't exist at the end
      #  not sure why since this is in for-loop rather than for-each
      rs <<- new_rs
    } else {
      cat( paste("\n\ndoParallel flag: ELSE"))
      rs <<- bind_rows(rs, new_rs)
    }
    
    cat( paste("\n\ndoParallel flag13. nrow(rs):" ) ); if (exists("rs")) nrow(rs)
    
  } # end "for ( scen in scens.to.run )"
  
} )[3]  # end system.time



cat( paste("\n\ndoParallel flag. Done entire for-loop." ) )



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


# temp debugging
cat( paste("\ndoParallel flag. head(rs):" ) )
print(head(rs))

rs$rep.seconds = doParallel.seconds/(sim.reps * length(scens.to.run))
rs$rep.seconds[ rs$method != unique(rs$method)[1] ] = NA


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
  

}




# ~ WRITE LONG AND SHORT RESULTS ------------------------------
if ( run.local == FALSE ) {
  cat("\n\n doParallel flag: about to write results \n")
  setwd("/home/groups/manishad/MBMA/long_results")
  fwrite( rs, paste( "long_results", jobname, ".csv", sep="_" ) )
  cat("\n\n doParallel flag: done writing results \n")
  
  
  
  # pre-aggregate 
  agg_job = make_agg_data(rs)
  setwd("/home/groups/manishad/MBMA/short_results")
  fwrite( agg_job, paste( "short_results", jobname, ".csv", sep="_" ) )
  
}


