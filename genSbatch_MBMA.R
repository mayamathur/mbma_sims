
# SET SIMULATION PARAMETERS MATRIX -----------------------------------------

# FOR CLUSTER USE
path = "/home/groups/manishad/MBMA"
setwd(path)
source("helper_MBMA.R")

allPackages = c("here",
                "magrittr",
                "dplyr",
                "data.table",
                "tidyverse",
                "tidyr",
                "metafor",
                "robumeta",
                "testthat",
                "truncdist",
                "gmm",
                "stringr",
                "tmvtnorm",
                "doParallel",
                "foreach")


( packagesNeeded = allPackages[ !( allPackages %in% installed.packages()[,"Package"] ) ] )
if( length(packagesNeeded) > 0 ) install.packages(packagesNeeded)

# load all packages
lapply( allPackages,
        require,
        character.only = TRUE)

#**you need to see all "TRUE" printed by this in order for the package to actually be loaded


### 2022-7-23 - debugging set ###
# RAN FINE WITH SIM.REPS = 100 (didn't try higher)
# scen.params = tidyr::expand_grid(
#   
#   rep.methods = "naive ; mbma-MhatB ; mbma-MhatB-gamma ; maon-adj-MhatB ; 2psm ; beta-sm",
#   
#   
#   # args from sim_meta_2
#   Nmax = 5,  # later code will set this to 1 if prob.hacked = 0
#   true.dist = c("expo", "norm"),
#   true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)",  # original setting close to empirical distribution
#                     "0.02 + rexp(n = 1, rate = 1)"),  # larger SEs overall
#   Mu = c(0, 0.5),
#   t2a = c(0, 0.25^2, 0.5^2),
#   t2w = c(0),
#   m = 50,
#   
#   # SWS args
#   # remember: method affirm will only work if prob.hacked < 1 else will never have nonaffirms
#   hack = c("affirm2", "favor-lowest-p", "favor-gamma-ratio"),  
#   rho = c(0),
#   prob.hacked = c(1, 0),
#   
#   # SAS args
#   k.pub.nonaffirm = c(5, 10, 15, 30, 50),
#   eta = c(1, 5, 10),
#   gamma = 2,  # only used for method favor-gamma-ratio
#   SAS.type = c("2psm", "carter"), # "2psm" (original) or "carter"
#   
#   
#   
#   # confounding parameters
#   # NOT using log scale here b/c underlying data are continuous
#   muB = c(0.1, 0.25, 0.5),
#   sig2B = 0.5,
#   prob.conf = c(0.5),
#   
#   # Stan control args - only relevant if running RTMA - remove these args?
#   stan.maxtreedepth = 25,
#   stan.adapt_delta = 0.995,
#   
#   get.CIs = TRUE,
#   run.optimx = FALSE )


### 2023-4-9 - full set ###
scen.params = tidyr::expand_grid(

  rep.methods = "naive ; mbma-MhatB ; mbma-MhatB-gamma ; 2psm ; beta-sm",


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
  hack = c("affirm2", "favor-lowest-p"),
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
  prob.conf = c(0.5, 1),

  # Stan control args - only relevant if running RTMA - remove these args?
  stan.maxtreedepth = 25,
  stan.adapt_delta = 0.995,

  get.CIs = TRUE,
  run.optimx = FALSE )


# ### 2022-7-23 - as in RSM_0 ###
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
#   # be careful about fractional values; might make overall eta hard to calculate 
#   prob.hacked = c(0),
#   
#   eta = c(1, 5, 10),
#   gamma = 2,
#   
#   true.sei.expr = c("0.02 + rexp(n = 1, rate = 3)"),
#   
#   # confounding parameters
#   # NOT using log scale here b/c underlying data are continuous
#   muB = c(0.25, 0.5),
#   sig2B = 0.5,
#   prob.conf = c(0.5), 
#   
#   # Stan control args - only relevant if running RTMA - remove these args?
#   stan.maxtreedepth = 25,
#   stan.adapt_delta = 0.995,
#   
#   get.CIs = TRUE,
#   run.optimx = FALSE )



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


( n.scen = nrow(scen.params) )
# look at it
head( as.data.frame(scen.params) )

# write the csv file of params (to Sherlock)
setwd(path)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
source("helper_MBMA.R")

# number of sbatches to generate (i.e., iterations within each scenario)
# n.reps.per.scen = 1000
#n.reps.in.doParallel = 1000
#( n.files = ( n.reps.per.scen / n.reps.in.doParallel ) * n.scen )

# sim.reps is now set only within doParallel
n.scens.per.doParallel = 15  # doesn't need to be set inside doParallel
( n.files = n.scen / n.scens.per.doParallel )


path = "/home/groups/manishad/MBMA"

#scen.name = rep( scen.params$scen, each = ( n.files / n.scen ) )
# group scenarios so each sbatch can receive multiple scen names (as a single string)
scen.name = group_scens(x = scen.params$scen,
                           n = n.scens.per.doParallel)

expect_equal( length(scen.name), n.files )

jobname = paste("job", 1:n.files, sep="_")
outfile = paste("/home/groups/manishad/MBMA/rmfiles/rm_", 1:n.files, ".out", sep="")
errorfile = paste("/home/groups/manishad/MBMA/rmfiles/rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")

sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            #  2023-4-9 sims: used 10:00:00 for initial run, which worked for ~75% of sbatch files
                            #  the others (mostly with large-k scens) timed out, so I then increased jobtime to 1 day
                            #  and re-ran the missed sbatches
                            #  info about job time limits: https://www.sherlock.stanford.edu/docs/advanced-topics/job-management/#job-submission-limits
                            # jobtime = "10:00:00", 
                            jobtime = "1-00:00:00", 
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel_MBMA.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

#@temp: only keep one of them
#sbatch_params = sbatch_params[1,]
generateSbatch(sbatch_params, runfile_path)

n.files

# run just the first one
# sbatch -p qsu,owners,normal /home/groups/manishad/MBMA/sbatch_files/1.sbatch

# 2023-04-09: 864
path = "/home/groups/manishad/MBMA"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 1:864) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/MBMA/sbatch_files/", i, ".sbatch", sep="") )
}



######## If Running Only Some Jobs To Fill Gaps ########

# run in Sherlock ml load R
path = "/home/groups/manishad/MBMA"
setwd(path)
source("helper_MBMA.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/MBMA/short_results",
                              "/home/groups/manishad/MBMA/short_results",
                              .name.prefix = "short_results",
                              .max.sbatch.num = 864 )



setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/MBMA/sbatch_files/", i, ".sbatch", sep="") )
}