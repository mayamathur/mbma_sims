
# EDIT: THIS RUNS OUT OF MEMORY. NOW USING THE NEW COPY.

# run this interactively in ml load R or via:
#   sbatch -p qsu,owners,normal /home/groups/manishad/MBMA/job_stitch.sbatch
# scontrol show job 33834701
# look at its out file:
# cd /home/groups/manishad/MBMA
# cd /home/users/mmathur
# less rm_stitch.out

# for non-huge simulations, can often run this script interactively in a higher-memory
#  Sherlock session:
# ml load R/4.1.2
# srun --mem=32G --time=3:00:00 --pty bash
# R


# to be run by stitch.sbatch or manually
# To quickly run this script in high-mem interactive session:
# setwd("/home/groups/manishad/MBMA"); source("stitch_on_sherlock_MBMA.R")

# # load command line arguments
# args = commandArgs(trailingOnly = TRUE)
# start.num = as.numeric( args[1] )  # starting results number to stitch
# stop.num = as.numeric( args[2] )  # stopping results number to stitch



path = "/home/groups/manishad/MBMA"
setwd(path)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")

# PRELIMINARIES ----------------------------------------------

library(data.table)
library(dplyr)
library(testthat)


# ##### NEW - stitch pre-aggregated things
# 
# .results.singles.path = "/home/groups/manishad/MBMA/short_results"
# .results.stitched.write.path = "/home/groups/manishad/MBMA/stitched_results"
# .name.prefix = "short_results"
# .stitch.file.name="agg.csv"
# 
# # get list of all files in folder
# all.files = list.files(.results.singles.path, full.names=TRUE)
# 
# # we only want the ones whose name includes .name.prefix
# keepers = all.files[ grep( .name.prefix, all.files ) ]
# length(keepers)
# 
# # grab variable names from first file
# names = names( read.csv(keepers[1] ) )
# 
# # read in and rbind the keepers
# tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )
# 
# # sanity check: do all files have the same names?
# # if not, could be because some jobs were killed early so didn't get doParallelTime
# #  variable added at the end
# #  can be verified by looking at out-file for a job without name "doParallelTime"
# allNames = lapply( tables, names )
# # # find out which jobs had wrong number of names
# # lapply( allNames, function(x) all.equal(x, names ) )
# # allNames[[1]][ !allNames[[1]] %in% allNames[[111]] ]
# 
# # bind_rows works even if datasets have different names
# #  will fill in NAs
# agg <- do.call(bind_rows, tables)
# 
# names(agg) = names( read.csv(keepers[1], header= TRUE) )
# 
# if( is.na(agg[1,1]) ) agg = agg[-1,]  # delete annoying NA row
# # write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
# 
# cat("\n\n nrow(agg) =", nrow(agg))
# cat("\n nuni(agg$scen.name) =", nuni(agg$scen.name) )








##### OLD - stitch individual rows
# MAKE STITCHED DATA ----------------------------------------------

# s = stitch_files(.results.singles.path = "/home/groups/manishad/MBMA/sim_results/long",
#                  .results.stitched.write.path = "/home/groups/manishad/MBMA/sim_results/overall_stitched",
#                  .name.prefix = "long_results",
#                  .stitch.file.name="stitched.csv")

.results.singles.path = "/home/groups/manishad/MBMA/long_results"
.results.stitched.write.path = "/home/groups/manishad/MBMA/stitched_results"
.name.prefix = "long_results"
.stitch.file.name="stitched.csv"


# get list of all files in folder
all.files = list.files(.results.singles.path, full.names=TRUE)

# we only want the ones whose name includes .name.prefix
keepers = all.files[ grep( .name.prefix, all.files ) ]
length(keepers)

# grab variable names from first file
names = names( read.csv(keepers[1] ) )

# read in and rbind the keepers
tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )

# sanity check: do all files have the same names?
# if not, could be because some jobs were killed early so didn't get doParallelTime
#  variable added at the end
#  can be verified by looking at out-file for a job without name "doParallelTime"
allNames = lapply( tables, names )
# # find out which jobs had wrong number of names
# lapply( allNames, function(x) all.equal(x, names ) )
# allNames[[1]][ !allNames[[1]] %in% allNames[[111]] ]

# bind_rows works even if datasets have different names
#  will fill in NAs
s <- do.call(bind_rows, tables)

names(s) = names( read.csv(keepers[1], header= TRUE) )

if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
# write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )

cat("\n\n nrow(s) =", nrow(s))
cat("\n nuni(s$scen.name) =", nuni(s$scen.name) )

# ~ Check for Bad Column Names ---------------------------

# not sure why this is needed - has NA columns at end
names(s)
any(is.na(names(s)))

if ( any(is.na(names(s))) ) {
  NA.names = which( is.na(names(s) ) )
  s = s[ , -NA.names ]
  
}

s = s %>% filter(!is.na(scen.name))



# ~ Quick Sanity Check  ---------------------------

# s %>% group_by(scen.name, method) %>%
#   summarise(n(),
#             mean(Mu),
#             meanNA(Mhat),
#             mean(is.na(Mhat)))

# ~ Write stitched.csv ---------------------------

setwd(.results.stitched.write.path)
fwrite(s, .stitch.file.name)

# also make a zipped version
string = paste("zip -m stitched.zip", .stitch.file.name)
system(string)




# MAKE AGG DATA ----------------------------------------------

path = "/home/groups/manishad/MBMA"
setwd(path)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")


agg = make_agg_data(s)

setwd(.results.stitched.write.path)
fwrite(agg, "agg.csv")


# table(agg$method.pretty)

cat("\n\n nrow(agg) =", nrow(agg))
cat("\n nuni(agg$scen.name) =", nuni(agg$scen.name) )

# # OPTIONAL: shorter version of agg that's nicer to look at
# agg.short = agg %>% select(prob.hacked, prob.conf, method, Mhat, MhatBias, MhatCover, MhatWidth,
#                            #MhatEstFail, MhatCIFail,
#                            sancheck.dp.k.nonaffirm.unhacked, sancheck.dp.k.nonaffirm.hacked) %>%
#   filter( !( method %in% c("rtma-adj-pmed", "rtma-adj-pmean") ) ) %>%
#   mutate_if(is.numeric, function(x) round(x,2))
#   
# setwd(.results.stitched.write.path)
# fwrite(agg.short, "agg_short.csv")




# # look again at failures
# t = agg %>% group_by(k.pub.nonaffirm, method) %>%
#   summarise( mean(MhatEstFail),
#              mean(MhatCIFail)
#              )
# as.data.frame(t)
# 
# 
# # errors of 2PSM when it fails
# table( s$overall.error[ s$method == "2psm" & is.na(s$Mhat) ] )



##### Move to Local #####

# # stitched and agg -> local directory
# scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/MBMA/stitched_results/* /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/2021/Sensitivity\ analysis\ for\ p-hacking\ \(MBMA\)/Linked\ to\ OSF\ \(MBMA\)/Sherlock\ simulation\ results/Pilot\ simulations

# LOOK FOR MISSED JOBS ----------------------------------------------


if (FALSE) {
  path = "/home/groups/manishad/MBMA"
  setwd(path)
  source("helper_MBMA.R")
  source("analyze_sims_helper_MBMA.R")
  
  # look for missed jobs
  missed.nums = sbatch_not_run( "/home/groups/manishad/MBMA/long_results",
                                "/home/groups/manishad/MBMA/long_results",
                                .name.prefix = "long",
                                .max.sbatch.num = 864)
  
  setwd( paste(path, "/sbatch_files", sep="") )
  for (i in missed.nums) {
    system( paste("sbatch -p qsu,owners,normal /home/groups/manishad/MBMA/sbatch_files/", i, ".sbatch", sep="") )
  }
}


