
# NOTES ---------------------------------------------------------------

# keeping this script in general Code dir because it's a living work in progress  

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study


# ANALYSIS FNS ---------------------------------------------------------------


# In order to catch errors from individual estimation methods safely and informatively,
#  in general the estimation method fns are structured st they can be run within the
#  fn run_method_safe, which automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error.

# ~~ Wrapper Fn to Safely Run a Method -------

# See note at the beginning of this script
#  this fn automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error

# Important: this fn works if method.fn() returns multiple rows
# BUT in that case, it assumes that the CIs are shared for all rows of that method

# expects global vars: all.errors, rep.res
# directly edits res via superassignment
run_method_safe = function( method.label,
                            method.fn,
                            .rep.res ) {
  
  cat( paste("\n run_method_safe flag 1: about to try running method", method.label) )
  
  
  tryCatch({
    
    method.output = method.fn()
    new.rows = method.output$stats
    
    if ( !exists("new.rows") ) {
      cat("\n\n**** Object new.rows didn't exist for method", method.label)
      cat("\nHere is method.output:\n")
      print(method.output)
    }
    
    cat( paste("\n run_method_safe flag 2: done calling method.fn() for", method.label) )
    
    error = NA
    
  }, error = function(err) {
    # needs to be superassignment because inside the "error" fn
    error <<- err$message
    
    # only need one variable in the blank dataframe since bind_rows below
    #  will fill in the rest
    new.rows <<- data.frame( method = method.label )
    
  })
  
  new.rows = new.rows %>% add_column( method = method.label, .before = 1 )
  new.rows$overall.error = error
  
  # optimx.dataframe is itself a df, so needs to be handled differently
  # if ( !is.null(optimx.dataframe) ) new.row = bind_cols(new.row, optimx.dataframe)
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.rows else .rep.res = bind_rows(.rep.res, new.rows)
  return(.rep.res) 
  
}


# ANALYSIS FNS: APPLIED EXAMPLES -----------------------------------------------


# ~ Fit Diagnostics -----------------------------------------------


# 2022-3-12
# fit diagnostics
# get CDF of (non-iid) marginal Z-scores (Zi.tilde)
#  given a fitted Shat
# .affirm: VECTOR with same length as x for affirm status
#  including the affirms is useful for 2PSM
yi_cdf = function(yi,
                  sei,
                  Mhat,
                  Shat) {
  
  
  affirm = (yi/sei) > qnorm(0.975)
  
  dat = data.frame( yi = yi,
                    sei = sei,
                    affirm = affirm,
                    cdfi = NA)
  
  if ( any(dat$affirm == FALSE) ) {
    dat$cdfi[ dat$affirm == FALSE ] = ptruncnorm(q = dat$yi[ dat$affirm == FALSE ],
                                                 a = -Inf,
                                                 b = qnorm(.975) * dat$sei[ dat$affirm == FALSE ],
                                                 mean = Mhat,
                                                 sd = sqrt(Shat^2 + sei[ dat$affirm == FALSE ]^2) )
  }
  
  if ( any(dat$affirm == TRUE) ) {
    dat$cdfi[ dat$affirm == TRUE ] = ptruncnorm(q = dat$yi[ dat$affirm == TRUE ],
                                                a = qnorm(.975) * dat$sei[ dat$affirm == TRUE ],
                                                b = Inf,
                                                mean = Mhat,
                                                sd = sqrt(Shat^2 + sei[ dat$affirm == TRUE ]^2))
  }
  
  return(dat)
  
}


# yi: published nonaffirmative estimates
# sei: their SEs
# Mhat: estimated Mu from RTMA
# Shat: estimated tau from RTMA
yi_qqplot = function(yi,
                     sei,
                     Mhat,
                     Shat){
  
  # get theoretical CDFs for each yi, given its affirm status
  cdf.dat = yi_cdf(yi = yi,
                   sei = sei,
                   Mhat = Mhat,
                   Shat = Shat)
  
  ecdf_fn = ecdf(yi)
  cdf.dat$ecdfi = ecdf_fn(yi)
  
  ggplot( data = cdf.dat,
          aes( x = cdfi,
               y = ecdfi) ) +
    geom_abline( slope = 1, 
                 intercept = 0,
                 color = "red") +
    geom_point( size = 2,
                alpha = 0.5 ) +
    xlab("Fitted CDF of point estimates") +
    ylab("Empirical CDF of point estimates") +
    theme_classic()
  
  
}

# DATA SIMULATION ---------------------------------------------------------------

# Simulate a meta-analysis in which some proportion of underlying studies (prior to publication)
#  are hacked, following various hacking mechanisms. Each study makes multiple draws (hypothesis tests)
#  until some stopping criterion based on Nmax and the hacking mechanism.

# - Nmax: max number of draws (hypothesis tests) that each hacked study can make before giving up
# - Mu: overall mean for meta-analysis
# - t2a: across-study heterogeneity (NOT total heterogeneity)
# Study parameters, assumed same for all studies:
#  - m: sample size
#  - t2w: within-study heterogeneity across draws
#  - true.sei.expr: quoted expression to evaluate to simulate a single study's standard error
#  - rho: autocorrelation of draws from a given study
#  - hack: mechanism of p-hacking (see sim_one_study_set for details)
# - k.pub.nonaffirm: number of published nonaffirmative studies desired in meta-analysis
#    (will simulate as many studies as needed to achieve that number)
# - prob.hacked: probability that an underlying study is hacked
# - SAS.type: "2psm" is what I had originally; "carter" is carter_censor in helper code
sim_meta_2 = function(Mu,
                      true.dist = "norm", # dist of population effects for this study set
                      t2a,  
                      k.pub.nonaffirm,  # number of published nonaffirmatives
                      return.only.published = FALSE,
                      
                      # study parameters, assumed same for all studies:
                      Nmax,
                      m,  # sample size for this study
                      t2w,  # within-study heterogeneity
                      true.sei.expr,  # TRUE SE string to evaluate
                      
                      # confounding arguments
                      muB = 0,
                      sig2B = 0,
                      prob.conf = 0,  # probability that an underlying study is confounded
                      
                      # SWS arguments
                      rho = 0,  # autocorrelation of muin's
                      hack, # mechanism of hacking for studies that DO hack (so not "no")
                      prob.hacked,
                      gamma, # only used for favor-gamma-ratio
                      
                      # SAS arguments
                      SAS.type = "2psm",  
                      eta = 1  # across-study SAS
                      
                      
) {
  
  
  # collect arguments to pass to sim_one_study_set
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k.pub.nonaffirm",
                                      "prob.hacked",
                                      "true.sei.expr",
                                      "eta",
                                      "prob.conf",
                                      "SAS.type")]
  
  
  if ( hack == "no" ) stop("hack should not be 'no' for this fn")
  
  k.pub.nonaffirm.achieved = 0
  i = 1
  
  # simulate studies until we reach the desired number of favored AND published nonaffirmatives
  while( k.pub.nonaffirm.achieved < k.pub.nonaffirm ) {
    
    Ci = rbinom(n = 1, size = 1, prob = prob.conf)
    is.hacked = rbinom(n = 1, size = 1, prob = prob.hacked)
    true.se = eval( parse( text = true.sei.expr) )
    
    if ( is.hacked == 0 ) {
      # to generate unhacked studies, need to change argument "hack"
      .argsUnhacked = .args
      .argsUnhacked[ names(.args) == "hack" ] = "no"
      .argsUnhacked$se = true.se
      .argsUnhacked$Ci = Ci      
      
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsUnhacked )
      
    } else if ( is.hacked == 1 ) {
      
      # for unhacked studies, no need to change argument "hack"
      .argsHacked = .args
      .argsHacked$se = true.se
      .argsHacked$Ci = Ci
      
      # might be multiple rows if return.only.published = FALSE
      newRows = do.call( sim_one_study_set, .argsHacked )
      
    }
    
    # add study ID
    newRows = newRows %>% add_column( .before = 1,
                                      study = i )
    
    # add study-draw ID
    study.draw = paste(newRows$study, 1:nrow(newRows), sep = "_")
    newRows = newRows %>% add_column( .after = 1,
                                      study.draw )
    
    
    #### Apply SAS to the study's favored draw only
    # IMPORTANT: Fi represents SWS selection indicator; Di.across represents OVERALL SWS & SAS indicator
    #   since it is set to 0 whenever Fi = 0
    newRows$Di.across.prob = rep(0, nrow(newRows))  # default to zero to avoid applying SAS to non-favored draws
    
    if ( SAS.type == "2psm" ) {
      newRows$Di.across.prob[ newRows$Fi == 1 & newRows$affirm == TRUE ] = 1
      newRows$Di.across.prob[ newRows$Fi == 1 & newRows$affirm == FALSE ] = 1/eta
      
    } else if ( SAS.type == "carter" ) {
      if ( sum(newRows$Fi) > 1 ) stop("Your hack type generates multiple favored draws per study, but SAS.type carter not implemented for this case")
      newRows$Di.across.prob[ newRows$Fi == 1 ] = carter_censor(pObs = newRows$pval[ newRows$Fi == 1 ],
                                                                direction = 1,  # always set direction to 1 regardless of sign for 2-tailed selection
                                                                posSign_NS_baseRate = 0.3,
                                                                negSign_NS_baseRate = 0.05,
                                                                counterSig_rate = 0.50)
    }
    newRows$Di.across = rbinom( n=nrow(newRows), size = 1, prob = newRows$Di.across.prob )
    #### done creating SAS & SWS indicator
    
    
    if ( i == 1 ) .dat = newRows else .dat = rbind( .dat, newRows )
    
    i = i + 1
    
    # the ".dat$Fi == 1" is included for clarity but is unnecessary given that Di.across == 0 whenever
    #   Fi == 0
    k.pub.nonaffirm.achieved = sum( .dat$affirm == FALSE & .dat$Fi == 1 & .dat$Di.across == 1 ) 
    
  }  # end "while( k.pub.nonaffirm.achieved < k.pub.nonaffirm )"
  
  
  # add more info to dataset
  .dat$k.underlying = length(unique(.dat$study))
  .dat$k.nonaffirm.underlying = length( unique( .dat$study[ .dat$affirm == FALSE ] ) )
  
  # record overall within-study eta (called gamma in SAPH)
  # can't do this within studies because either the numerator or denominator will be zero
  #  (b/c only 1 favored draw per study)
  # all hack types need to return Fi indicator
  # BEWARE: need to calculate overall gamma BEFORE restricting to studies, as we do here
  .dat$gamma = conditional_prob(.dat$Fi == 1, .dat$affirm == 1) / conditional_prob(.dat$Fi == 1, .dat$affirm == 0 )
  
  
  
  # return only draws that are both favored AND published
  if ( return.only.published == TRUE ) {
    .dat = .dat[ .dat$Fi == 1 & .dat$Di.across == 1, ]
  }
  
  return(.dat)
}



# # confounding only
# d = sim_meta_2(  # test only
#   Nmax = 1,
#   Mu = 1,
#   t2a = 0.1,
#   m = 50,
#   t2w = .5,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   k.pub.nonaffirm = 30,
#   prob.hacked = 0,
#   eta = 1,
#   muB = 1,
#   sig2B = 0.5,
#   prob.conf = 0.5)
# 
# d %>% group_by(Ci) %>%
#   summarise(mean(Bi),
#             mean(mui),
#             mean(yi))


# # confounding + SAS
# d = sim_meta_2(  # test only
#   Nmax = 1,
#   Mu = 0.5,
#   t2a = sqrt(0.2),
#   m = 50,
#   t2w = .20,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   k.pub.nonaffirm = 100,
#   prob.hacked = 0,
#   eta = 5,
#   muB = log(2),
#   sig2B = 0.5,
#   prob.conf = 0.5)
# 
# # before SAS
# d %>%
#   group_by(Ci, affirm) %>%
#   summarise(n(),
#             mean(Bi),
#             mean(mui),
#             mean(yi))
# 
# # after SAS
# d %>%
#   filter(Di == 1 & Di.across == 1) %>%
#   group_by(Ci, affirm) %>%
#   summarise(n(),
#             mean(Bi),
#             mean(mui),
#             mean(yi))


# # SAS only
# d = sim_meta_2(  # test only
#   Nmax = 1,
#   Mu = 1,
#   t2a = 0.1,
#   m = 50,
#   t2w = .5,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   k.pub.nonaffirm = 30,
#   prob.hacked = 0,
#   eta = 5)
# 
# # should be 1 and 1/eta
# d %>% group_by(affirm) %>%
#   summarise(mean(Di.across))


# 
# # with SWS
# d = sim_meta_2(  
#   Nmax = 20,
#   Mu = 1,
#   t2a = 0.1,
#   m = 50,
#   t2w = .5,
#   true.sei.expr = "runif( n = 1, min = 0.5, max = 2 )",
#   hack = "affirm",
#   return.only.published = FALSE,
#   rho=0,
#   k.pub.nonaffirm = 30,
#   prob.hacked = 0.4 )
# 
# d$k.nonaffirm.underlying[1]
# d$k.underlying[1]
# table(d$Di, d$affirm)




# ~ Simulate a single study ----------------- 

# Simulate study from potentially heterogeneous meta-analysis distribution;
#  within-study draws have their own heterogeneous within-study distribution

### Hacking types ###

# - "no": Makes exactly Nmax results and treats the last one as the reported one,
#     so the final result could be affirmative or nonaffirmative

# - "affirm" (worst-case hacking): Makes draws until the first affirmative is obtained
#    but if you reach Nmax, do NOT report any result at all. (Hack type "signif" is the 
#    same but favors all significant results.)

# - "affirm2": (NOT worst-case hacking): Similar to "affirm", but always reports the last draw,
#    even if it was nonaffirm (no file drawer)

# - "favor-best-affirm-wch" (worst-case hacking): Always makes Nmax draws. If you get any affirmatives,
#    publish the one with the lowest p-value. If you don't get any affirmatives, don't publish anything.

# - "favor-lowest-p": Always make Nmax draws. Favor the one with lowest p-value. (Collect data on P(Fin* = 1 | Ain*=1) / P(Fin* = 1 | Ain*=0), the within-study contribution to eta.)

# - "favor-gamma-ratio" (staackable hacking): Always make Nmax draws. Decide whether to favor each draw based on within-study selection ratio, gamma. As such, studies can have multiple favored draws. (2PSM will be correctly specificed, I think also subject to constraints on heterogeneity - see SAPH).

# If Nmax is small, rhoEmp (empirical autocorrelation of muin's) will be smaller
#  than rho. That's okay because it reflects small-sample bias in autocorrelation
# estimate itself, not a problem with the simulation code
# For more about the small-sample bias: # https://www.jstor.org/stable/2332719?seq=1#metadata_info_tab_contents

sim_one_study_set = function(Nmax,  # max draws to try
                             true.dist = "norm", # dist of population effects for this study set
                             Mu,  # overall mean for meta-analysis
                             t2a,  # across-study heterogeneity (NOT total heterogeneity)
                             m,  # sample size for this study
                             t2w,  # within-study heterogeneity
                             se,  # TRUE SE for this study
                             return.only.published = FALSE,
                             hack, # should this study set be hacked? ("no", "affirm","affirm2", "signif")
                             
                             muB = 0,
                             sig2B = 0,
                             Ci = 0,  # should this study set be confounded?
                             
                             # for correlated draws; see make_one_draw
                             rho = 0,
                             
                             # within-study selection ratio; only used when hack=favor-gamma-ratio
                             # important: other hacking methods will IGNORE gamma because can't be directly specified
                             gamma = 1
) {  
  
  
  # # test only
  # Nmax = 1
  # Mu = -0.15
  # t2a = 0.1
  # m = 50
  # t2w = .5
  # se = 1
  # hack = "no"
  # muB = 0.2
  # sig2B = 0
  # Ci = 1  # should this study set be confounded?
  # rho=0
  
  # ~~ Mean (potentially confounded) for this study set ----
  
  # bias factor due to confounding
  if ( Ci == 0 ) Bi = 0
  if ( Ci == 1 ) Bi = rnorm( mean = muB, sd = sqrt(sig2B), n = 1 )
  
  # doesn't have t2w because that applies to results within this study set
  if ( true.dist == "norm" ){
    mui = Mu + Bi + rnorm(mean = 0,
                          sd = sqrt(t2a),
                          n = 1)
  }
  
  if ( true.dist == "expo" ){
    # set the rate so the heterogeneity is correct
    mui = rexp( n = 1, rate = sqrt(1/t2a) )
    # now the mean is sqrt(t2a) rather than Mu + Bi
    # shift to have the correct mean (in expectation)
    mui = mui + ( (Mu + Bi) - sqrt(t2a))
    
  }
  
  
  # TRUE SD (not estimated)
  sd.y = se * sqrt(m)
  
  # collect all args from outer fn, including default ones
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  .args$mui = mui
  .args$sd.y = sd.y
  
  
  stop = FALSE  # indicator for whether to stop drawing results
  N = 0  # counts draws actually made
  
  # ~~ Draw until study reaches its stopping criterion ----
  # we use this loop whether there's hacking or not
  while( stop == FALSE & N < Nmax ) {
    
    if ( rho == 0 ) {
      # make uncorrelated draw
      newRow = do.call( make_one_draw, .args )
    } else {
      # make correlated draw
      if ( N == 0 ) .args$last.muin = NA  # on first draw, so there's no previous one
      if ( N > 0 ) .args$last.muin = d$muin[ nrow(d) ]
      newRow = do.call( make_one_draw, .args ) 
    }
    
    
    # number of draws made so far
    N = N + 1
    
    # add new draw to dataset
    if ( N == 1 ) d = newRow
    if ( N > 1 ) d = rbind( d, newRow )
    
    # check if it's time to stop drawing results
    if (hack == "signif") {
      stop = (newRow$pval < 0.05)
    } else if ( hack %in% c("affirm", "affirm2") ) {
      stop = (newRow$pval < 0.05 & newRow$yi > 0)
    } else if ( hack %in% c("no", "favor-best-affirm-wch", "favor-lowest-p", "favor-gamma-ratio") ) {
      # if this study set is unhacked, then stopping criterion
      #  is just whether we've reached Nmax draws
      # and for favor-best-affirm-wch or favor-lowest-p, we always do Nmax draws so 
      #  we can pick the smallest p-value
      stop = (N == Nmax)
    } else {
      stop("No stopping criterion implemented for your chosen hack mechanism")
    }
    
  }  # end while-loop until N = Nmax or we succeed
  
  # record info in dataset
  d$N = N
  d$hack = hack
  d$Ci = Ci
  d$Bi = Bi
  
  # ~~ Empirical correlation of muin's ----
  #  but note this will be biased for rho in small samples (i.e., Nmax small)
  if ( nrow(d) > 1 ) {
    # get lag-1 autocorrelation
    d$rhoEmp = suppressWarnings( cor( d$muin[ 2:length(d$muin) ],
                                      d$muin[ 1: ( length(d$muin) - 1 ) ] ) )
    
    # mostly for debugging; could remove later
    d$covEmp = suppressWarnings( cov( d$muin[ 2:length(d$muin) ],
                                      d$muin[ 1: ( length(d$muin) - 1 ) ] ) )
    
  } else {
    d$rhoEmp = NA
    d$covEmp = NA
  }
  
  # convenience indicators for significance and affirmative status
  d$signif = d$pval < 0.05
  d$affirm = (d$pval < 0.05 & d$yi > 0)
  
  
  # ~~ Decide which draw to favor & publish ----
  # for these hack types, Fi=1 for only the last draw IF we got an affirm result
  #  but if we didn't, then it will always be 0
  if ( hack == "signif" ) d$Fi = (d$signif == TRUE)
  if ( hack == "affirm" ) d$Fi = (d$affirm == TRUE)
  
  # CHANGED Di to Fi
  # if no hacking or affirmative hacking without file drawer,
  #   assume only LAST draw is published,
  #   which could be affirm or nonaffirm
  if ( hack %in% c("no", "affirm2") ) {
    d$Fi = 0
    d$Fi[ length(d$Fi) ] = 1
  }
  
  # CHANGED Di to Fi
  # for favor-best-affirm-wch, consider all AFFIRMATIVE draws and favor the one
  #  with the lowest p-value
  if ( hack %in% c("favor-best-affirm-wch") ) {
    d$Fi = 0
    # if there was at least 1 affirm, publish it 
    if ( any(d$affirm == TRUE) ) {
      best.affirm.pval = min( d$pval[d$affirm == TRUE] )
      d$Fi[ d$pval == best.affirm.pval & d$affirm == TRUE ] = 1
    }
    # ...otherwise don't publish any draw
    # sanity check:
    #View(d%>%select(Di,affirm,pval,yi))
  }
  
  
  # NEW - using Fi as within-study indicator
  # for favor-lowest-p, consider all draws and favor the one with the lowest two-tailed p-value regardless of estimate sign
  if ( hack %in% c("favor-lowest-p") ) {
    d$Fi = 0
    best.pval = min(d$pval)
    d$Fi[ d$pval == best.pval ] = 1
  }
  
  if ( hack %in% c("favor-gamma-ratio") ) {
    d = d %>% rowwise() %>%
      mutate( Fi.prob = ifelse(affirm == TRUE, 1, 1/gamma),
              Fi = rbinom(n = 1, size = 1, prob = Fi.prob) )
    
  }
  
  if ( return.only.published == TRUE ){
    stop("return.only.published==TRUE not implemented")
    #  d = d[ d$Di == 1, ]
  }
  
  return(d)
  
}


# two boolean vectors
# P(vec1 == 1 | vec2 == 1)
conditional_prob = function(vec1, vec2){
  
  if ( mean(vec2 == 1) > 0 ){
    return( mean( vec1 == 1 & vec2 == 1 ) / mean(vec2 == 1) )
  } else {
    return(NA)
  }
}


### example
# 
# d = sim_one_study_set(Nmax = 20,
#                       Mu = 1,
#                       t2a = 1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "favor-lowest-p",
#                       return.only.published = FALSE)
# nrow(d)
# d

### example
# 
# d = sim_one_study_set(Nmax = 20,
#                       Mu = 1,
#                       t2a = 1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "favor-best-affirm-wch",
#                       return.only.published = FALSE)
# nrow(d)
# d



# ### example
# d = sim_one_study_set(Nmax = 5,
#                       Mu = 0.1,
#                       t2a = 0.1,
#                       m = 50,
#                       t2w = .5,
#                       se = 1,
#                       hack = "affirm2",
#                       return.only.published = FALSE)
# nrow(d)
# d




# ~ Draw one unbiased result within one study ------------------
# muin should be NA if either we want uncorrelated draws OR it's the first of a series of correlated draws
make_one_draw = function(mui,  # mean for this study set
                         t2w,
                         sd.y,  # TRUE SD
                         m,  # sample size
                         
                         # for making correlated draws
                         rho = 0,  # autocorrelation of muin's (not yi's)
                         last.muin = NA,  
                         ...) {
  
  
  # true mean for draw n (based on within-study heterogeneity)
  # either we want uncorrelated draws OR it's the first of a series of correlated draws
  # note that the within-study draws are normal regardless of true.dist
  #  (that arg only controls the study SET's population effect)
  if ( rho == 0 | is.na(last.muin) ) {
    muin = rnorm(mean = mui,
                 sd = sqrt(t2w),
                 n = 1)
  }
  
  # make correlated draw
  if ( rho != 0 & !is.na(last.muin) ) {
    # draw from BVN conditional, given muin from last draw
    # conditional moments given here:
    #  https://www.amherst.edu/system/files/media/1150/bivarnorm.PDF
    muin = rnorm(mean = mui + rho * (last.muin - mui),
                 sd = sqrt( t2w * (1 - rho^2) ),
                 n = 1)
  }
  
  
  # draw subject-level data from this study's population effect
  y = rnorm( mean = muin,
             sd = sd.y,
             n = m)
  
  # run a one-sample t-test
  # two-tailed p-value
  test = t.test(y,
                alternative = "two.sided")
  
  pval = test$p.value
  tstat = test$statistic
  vi = test$stderr^2  # ESTIMATED variance
  
  
  # if (hack == "signif") success = (pval < 0.05)
  # if (hack == "affirm") success = (pval < 0.05 & mean(y) > 0)
  
  return( data.frame(pval = pval,
                     tstat = tstat,
                     tcrit = qt(0.975, df = m-1),
                     mui = mui,
                     muin = muin,
                     yi = mean(y),
                     vi = vi,
                     viTrue = sd.y^2 / m,  # true variance; will equal p$se^2
                     m = m ) )
  #success = success,
  #N = Nmax ) )
}


# ~ Carter pub bias mechanism ---------------------------------------------------------------

# Hong & Reed's implementation: https://osf.io/wnyhg

# clamp x to values between 0 and 1
clamp <- function(x) {if (x > 1) x=1; if (x < 0) x = 0; return(x)}


#' @param p observed p value
#' @param p_range Range of observed p-values where the easing takes place
#' @param from_prob Probability of publication (starting position)
#' @param to_prob Probability of publication (end position)
easeOutExpo <- function (p, p_range, from_prob, to_prob) {
  p_start <- p_range[1]
  p_range_length <- p_range[2] - p_range[1]
  (to_prob-from_prob) * (-2^(-10 * (p-p_range[1])/p_range_length) + 1) + from_prob;
}

easeInExpo <- function (p, p_range, from_prob, to_prob) {
  p_start <- p_range[1]
  p_range_length <- p_range[2] - p_range[1]
  (to_prob-from_prob) * 2^(10 * (((p-p_range[1])/p_range_length) - 1)) + from_prob;
}


#' @param pObs two-tailed p-value
#' @param posSign_NS_baseRate What's the probability that a p > .10 in the right direction enters the literature?
#' @param negSign_NS_baseRate What's the probability that a p > .01 in the wrong direction enters the literature? (Left anchor at p = .01)
#' @param counterSig_rate What's the probability that a p < .001 in the wrong direction enters the literature?
#' @param direction +1: Expected direction, -1: wrong direction

carter_censor <- function(pObs, direction, posSign_NS_baseRate = 0.3, negSign_NS_baseRate = 0.05, counterSig_rate = 0.50){
  
  # ---------------------------------------------------------------------
  # Correct direction of effect
  
  if (direction > 0 & pObs < .05){       #right direction, sig
    pubProb = 1
  } else if(direction > 0 & pObs >= .05 & pObs < .1){ #right direction, trending
    pubProb = easeOutExpo(p=pObs, p_range=c(.05, .1), from_prob=1, to_prob=posSign_NS_baseRate)
  } else if (direction > 0 & pObs >= .1){	# right direction; non-significant (p > .1)
    pubProb =posSign_NS_baseRate
    
    # ---------------------------------------------------------------------
    # Wrong direction of effect	
  } else if(direction <= 0 & pObs < .001){	# wrong direction, highly sig.
    pubProb= counterSig_rate
  } else if(direction <= 0 & pObs >= .001 & pObs < .01){ # wrong direction, standard sig.
    pubProb= easeOutExpo(p=pObs, p_range=c(.001, .01), from_prob=counterSig_rate, to_prob=negSign_NS_baseRate)
  } else if(direction <= 0 & pObs >= .01){	# wrong direction, non-sig.
    pubProb=negSign_NS_baseRate
  }
  return(pubProb)
}

# carter_censor(pObs = 0.1,
#        direction = 1,
#        posSign_NS_baseRate = 0.3,
#        negSign_NS_baseRate = 0.05,
#        counterSig_rate = 0.50)



# DATA WRANGLING ---------------------------------------------------------------

# corrObject: something returned by correct_dataset_phack
# looks for (or makes) global object, "res"
add_method_result_row = function(repRes = NA,
                                 corrObject,
                                 methName) {
  
  # newRow = bind_cols( corrObject$metaCorr,
  #                 corrObject$sanityChecks )
  #TEMP: DON'T KEEP THE SANITY CHECKS BECAUSE CORRECT_META_PHACK2 doesn't have it
  newRow = corrObject$metaCorr
  
  newRow = newRow %>% add_column(.before = 1,
                                 methName = methName )
  
  
  # "if" condition is hacky way to deal with repRes = NA case
  if ( is.null( nrow(repRes) ) ) repRes = newRow else repRes = bind_rows(repRes, newRow)
  return(repRes)
}



# quickly look at results when running doParallel locally
srr = function() {
  
  if( "optimx.Mhat.winner" %in% names(rep.res) ) {
    cat("\n")
    print( rep.res %>% select(method, Mhat, MLo, MHi,
                              Shat,
                              optim.converged,
                              optimx.Mhat.winner,
                              optimx.Nconvergers,
                              optimx.Pagree.of.convergers.Mhat.winner) %>%
             mutate_if(is.numeric, function(x) round(x,2)) ) %>%
      cat("\n")
  } else {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) %>%
             mutate(MhatWidth = MHi - MLo))
    cat("\n")
  }
}



# nicely report a metafor or robumeta object with optional suffix to denote which model
report_meta = function(.mod,
                       .mod.type = "rma",  # "rma" or "robu"
                       .suffix = "") {
  
  if ( !is.null(.mod) ) {

    if ( .mod.type == "rma" ) {
      tau.CI = tau_CI(.mod)
      .res = data.frame( .mod$b,
                         .mod$ci.lb,
                         .mod$ci.ub,
                         
                         sqrt(.mod$tau2),
                         tau.CI[1],
                         tau.CI[2] )
    }
    
    
    if ( .mod.type == "robu" ) {
      
      .res = data.frame( .mod$b.r,
                         .mod$reg_table$CI.L,
                         .mod$reg_table$CI.U,
                         
                         sqrt(.mod$mod_info$tau.sq),
                         NA,
                         NA )
    }
    
  } else {
    .res = data.frame( rep(NA, 6) )
  }
  
  
  names(.res) = paste( c("Mhat", "MLo", "MHi", "Shat", "SLo", "SHi"), .suffix, sep = "" )
  row.names(.res) = NULL
  

  return( list(stats = .res) )
}

# SMALL GENERIC HELPERS ---------------------

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}


# check CI coverage
covers = function( truth, lo, hi ) {
  return( (lo <= truth) & (hi >= truth) )
}

# get names of dataframe containing a string
namesWith = function(pattern, dat){
  names(dat)[ grepl(pattern = pattern, x = names(dat) ) ]
}


# quick length(unique)
nuni = function(x) {
  length(unique(x))
}

# (re-)install package AND its dependencies
# useful for stupid rstan issues in which rstan itself it UTD but not its dependencies
# https://stackoverflow.com/questions/21010705/update-a-specific-r-package-and-its-dependencies
instPkgPlusDeps <- function(pkg, install = FALSE,
                            which = c("Depends", "Imports", "LinkingTo"),
                            inc.pkg = TRUE) {
  stopifnot(require("tools")) ## load tools
  ap <- available.packages() ## takes a minute on first use
  ## get dependencies for pkg recursively through all dependencies
  deps <- package_dependencies(pkg, db = ap, which = which, recursive = TRUE)
  ## the next line can generate warnings; I think these are harmless
  ## returns the Priority field. `NA` indicates not Base or Recommended
  pri <- sapply(deps[[1]], packageDescription, fields = "Priority")
  ## filter out Base & Recommended pkgs - we want the `NA` entries
  deps <- deps[[1]][is.na(pri)]
  ## install pkg too?
  if (inc.pkg) {
    deps = c(pkg, deps)
  }
  ## are we installing?
  if (install) {
    install.packages(deps)
  }
  deps ## return dependencies
}

# example
# instPkgPlusDeps("fields")


# CLUSTER FNS ---------------------------------------------------------------

# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
sbatch_skeleton <- function() {
  return(
    "#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=JOBNAME
#################  
#a file for job output, you can check job progress
#SBATCH --output=OUTFILE
#################
# a file for errors from the job
#SBATCH --error=ERRORFILE
#################
#time you think you need; default is one hour
#SBATCH --time=JOBTIME
#################
#quality of service; think of it as job priority
#SBATCH --qos=QUALITY
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners,qsu
#################
#number of nodes you are requesting
#SBATCH --nodes=NODENUMBER
#################
#memory per node; default is 4000 MB
#SBATCH --mem=MEMPERNODE
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=MAILTYPE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=USER_EMAIL
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=TASKS_PER_NODE
#################
#SBATCH --cpus-per-task=CPUS_PER_TASK
#now run normal batch commands

ml load v8
ml load R/4.1.2
R -f PATH_TO_R_SCRIPT ARGS_TO_R_SCRIPT")
}



generateSbatch <- function(sbatch_params,
                           runfile_path = NA,
                           run_now = F) {
  
  #sbatch_params is a data frame with the following columns
  #jobname: string, specifies name associated with job in SLURM queue
  #outfile: string, specifies the name of the output file generated by job
  #errorfile: string, specifies the name of the error file generated by job
  #jobtime: string in hh:mm:ss format, max (maybe soft) is 48:00:00 
  #specifies the amoung of time job resources should be allocated
  #jobs still running after this amount of time will be aborted
  #quality: kind of like priority, normal works
  #node_number, integer: the number of nodes (computers w/16 cpus each) to allocate 
  #mem_per_node, integer: RAM, in MB, to allocate to each node
  #mailtype, string: ALL, BEGIN, END, FAIL: what types of events should you be notified about via email
  #user_email string: email address: email address to send notifications
  #tasks_per_node: integer, number of tasks, you should probably use 1
  #cpus_per_task: integer, 1-16, number of cpus to use, corresponds to number of available cores per task
  #path_to_r_script: path to r script on sherlock
  #args_to_r_script: arguments to pass to r script on command line
  #write_path: where to write the sbatch file
  #server_sbatch_path: where sbatch files will be stored on sherlock
  #runfile_path is a string containing a path at which to write an R script that can be used to run
  #the batch files generated by this function. 
  #if NA, no runfile will be written
  #run_now is a boolean specifying whether batch files should be run as they are generated
  
  sbatches <- list()
  if (!is.na(runfile_path)) {
    outfile_lines <- c(paste0("# Generated on ",  Sys.time()))
  }
  for (sbatch in 1:nrow(sbatch_params) ) {
    gen_batch <- sbatch_skeleton()
    #set job name
    if (is.null(sbatch_params$jobname[sbatch])) { 
      gen_batch <- gsub("JOBNAME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBNAME", sbatch_params$jobname[sbatch], gen_batch) 
    }
    #set outfile name
    if (is.null(sbatch_params$outfile[sbatch])) { 
      gen_batch <- gsub("OUTFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("OUTFILE", sbatch_params$outfile[sbatch], gen_batch) 
    }
    #set errorfile name
    if (is.null(sbatch_params$errorfile[sbatch])) { 
      gen_batch <- gsub("ERRORFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ERRORFILE", sbatch_params$errorfile[sbatch], gen_batch) 
    }
    #set jobtime
    if (is.null(sbatch_params$jobtime[sbatch])) { 
      gen_batch <- gsub("JOBTIME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBTIME", sbatch_params$jobtime[sbatch], gen_batch) 
    }
    #set quality
    if (is.null(sbatch_params$quality[sbatch])) { 
      gen_batch <- gsub("QUALITY", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("QUALITY", sbatch_params$quality[sbatch], gen_batch) 
    }
    #set number of nodes
    if (is.null(sbatch_params$node_number[sbatch])) { 
      gen_batch <- gsub("NODENUMBER", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("NODENUMBER", sbatch_params$node_number[sbatch], gen_batch) 
    }
    #set memory per node
    if (is.null(sbatch_params$mem_per_node[sbatch])) { 
      gen_batch <- gsub("MEMPERNODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MEMPERNODE", sbatch_params$mem_per_node[sbatch], gen_batch) 
    }
    #set requested mail message types
    if (is.null(sbatch_params$mailtype[sbatch])) { 
      gen_batch <- gsub("MAILTYPE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MAILTYPE", sbatch_params$mailtype[sbatch], gen_batch) 
    }
    #set email at which to receive messages
    if (is.null(sbatch_params$user_email[sbatch])) { 
      gen_batch <- gsub("USER_EMAIL", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("USER_EMAIL", sbatch_params$user_email[sbatch], gen_batch) 
    }
    #set tasks per node
    if (is.null(sbatch_params$tasks_per_node[sbatch])) { 
      gen_batch <- gsub("TASKS_PER_NODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("TASKS_PER_NODE", sbatch_params$tasks_per_node[sbatch], gen_batch) 
    }
    #set cpus per task
    if (is.null(sbatch_params$cpus_per_task[sbatch])) { 
      gen_batch <- gsub("CPUS_PER_TASK", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("CPUS_PER_TASK", sbatch_params$cpus_per_task[sbatch], gen_batch) 
    }
    #set path to r script
    if (is.null(sbatch_params$path_to_r_script[sbatch])) { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", sbatch_params$path_to_r_script[sbatch], gen_batch) 
    }
    #set args to r script
    if (is.null(sbatch_params$args_to_r_script[sbatch])) { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", sbatch_params$args_to_r_script[sbatch], gen_batch) 
    }
    
    #write batch file
    if (is.null(sbatch_params$write_path[sbatch])) { 
      cat(gen_batch, file = paste0("~/sbatch_generated_at_", gsub(" |:|-", "_", Sys.time()) ), append = F)
    } else { 
      cat(gen_batch, file = sbatch_params$write_path[sbatch], append = F)
    }
    
    if (!is.na(sbatch_params$server_sbatch_path[sbatch])) {
      outfile_lines <- c(outfile_lines, paste0("system(\"sbatch ", sbatch_params$server_sbatch_path[sbatch], "\")"))
    } 
    sbatches[[sbatch]] <- gen_batch
  }
  if (!is.na(runfile_path)) {
    cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
  }
  if(run_now) { system(paste0("R -f ", runfile_path)) } 
  
  return(sbatches)
}


# looks at results files to identify sbatches that didn't write a file
# .max.sbatch.num: If not passed, defaults to largest number in actually run jobs.

sbatch_not_run = function(.results.singles.path,
                          .results.write.path,
                          .name.prefix,
                          .max.sbatch.num = NA ) {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # extract job numbers
  sbatch.nums = as.numeric( unlist( lapply( strsplit( keepers, split = "_"), FUN = function(x) x[5] ) ) )
  
  # check for missed jobs before the max one
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums)
  all.nums = 1 : .max.sbatch.num
  missed.nums = all.nums[ !all.nums %in% sbatch.nums ]
  
  # give info
  print( paste("The max job number is: ", max(sbatch.nums) ) )
  print( paste( "Number of jobs that weren't run: ",
                ifelse( length(missed.nums) > 0, length(missed.nums), "none" ) ) )
  
  if( length(missed.nums) > 0 ) {
    setwd(.results.write.path)
    write.csv(missed.nums, "missed_job_nums.csv")
  }
  
  return(missed.nums)
  
}


# for running multiple scens in each doParallel 
group_scens <- function(x, n) {
  # Calculate the number of groups
  num_groups <- ceiling(length(x) / n)
  
  # Split the vector into groups of length n
  groups <- split(x, rep(1:num_groups, each=n, length.out=length(x)))
  
  # Combine the groups into strings
  strings <- sapply(groups, function(group) {
    paste0(group, collapse=",")
  })
  
  return( as.vector(strings) )
}

# x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# group_scens(x, 5)


# FN: STITCH RESULTS FILES -------------------------------------

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # .results.singles.path = "/home/groups/manishad/MRM/sim_results/long"
  # .results.stitched.write.path = "/home/groups/manishad/MRM/sim_results/overall_stitched"
  # .name.prefix = "long_results"
  # .stitch.file.name="stitched.csv"
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  names = names( read.csv(keepers[1] )[-1] )
  
  # read in and rbind the keepers
  tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )
  s <- do.call(rbind, tables)
  
  names(s) = names( read.csv(keepers[1], header= TRUE) )
  
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
  return(s)
}


# quickly look at results from job #1

res1 = function() {
  setwd("/home/groups/manishad/SAPH/long_results")
  rep.res = fread("long_results_job_1_.csv")
  srr()
  
  cat("\nErrors by method:" )
  print( rep.res %>% group_by(method) %>%
           summarise(prop.error = mean( overall.error != "" ) ) )
  
  #table(rep.res$overall.error)
  
  cat("\n\nDim:", dim(rep.res))
  cat("\n\nReps completed:", nrow(rep.res)/nuni(rep.res$method))
}


