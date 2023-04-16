

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




# PRELIMINARIES ----------------------------------------------

library(data.table)
library(dplyr)
library(testthat)
library(doParallel)

# set the number of cores
num_cores <- 16
registerDoParallel(num_cores)

path = "/home/groups/manishad/MBMA"
setwd(path)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")




# MAKE AGG DATA - FROM PRE-AGGREGATED SHORT RESULTS ----------------------------------------------


# can be run in interactive session - :)!
.results.singles.path = "/home/groups/manishad/MBMA/short_results"
.results.stitched.write.path = "/home/groups/manishad/MBMA/stitched_results"
.name.prefix = "short_results"
.stitch.file.name="agg.csv"


# get list of all files in folder
all.files = list.files(.results.singles.path, full.names=TRUE)

# we only want the ones whose name includes .name.prefix
keepers = all.files[ grep( .name.prefix, all.files ) ]
length(keepers)

# grab variable names from first file
names = names( read.csv(keepers[1] ) )

#bm
# read in and rbind the keepers in parallel
rbind_tables <- function(files) {
  tables <- lapply(files, fread)
  rbindlist(tables, fill = TRUE)
}

split_files <- split(keepers, 1:length(keepers) %% num_cores)

tables <- foreach(i = 1:length(split_files), .packages = c("data.table")) %dopar% {
  rbind_tables(split_files[[i]])
}

# combine the results
agg <- rbindlist(tables, fill = TRUE)



# write it
setwd(.results.stitched.write.path)
fwrite(agg, .stitch.file.name)



# ZIP THE SHORT RESULTS ----------------------------------------------


# zip the results (can be run from any directory on Sherlock)
# zip -r /home/groups/manishad/MBMA/stitched_results/short_results.zip /home/groups/manishad/MBMA/short_results



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


