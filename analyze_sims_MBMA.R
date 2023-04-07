
# Note: This script shares a stats_for_paper.csv with the applied examples, so if you wr(), 
#  you'll also need to rerun the applied examples.


# PRELIMINARIES ----------------------------------------------------

#rm(list=ls())

# data-wrangling packages
library(here)
library(plotly)  # must be BEFORE dplyr or else plotly::select will take over
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(tidyverse)
library(fastDummies)
library(xlsx)
# meta-analysis packages
library(metafor)
library(robumeta)
# other
library(xtable)
library(testthat)
library(Deriv)
library(mosaic)
library(hpa)
library(pracma)
library(truncnorm)
library(tmvtnorm)
library(RColorBrewer)
library(sjmisc)

# prevent masking
select = dplyr::select


# ~~ User-specified global vars -------------------------
# no sci notation
options(scipen=999)

# control which results should be redone and/or overwritten
# e.g., sim_plot_multiple_outcomes
#@ not all fns respect this setting
overwrite.res = TRUE


# ~~ Set directories -------------------------
code.dir = here()

data.dir = str_replace( string = here(),
                        pattern = "Code",
                        replacement = "Results/temp" )


# temp while sims are ongoing
results.dir = data.dir
# results.dir = str_replace( string = here(),
#                            pattern = "Code",
#                            replacement = "Results/*2022-7-23 More scens for manuscript; add mbma-MhatB-true-t2" )


#@temp: commented out to avoid overwriting
# overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/figures/sims"
# 
# # for stats_for_paper.csv
# # same dir as for applied examples so that they'll write to single file
# overleaf.dir.nums = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/R_objects"


# # alternative for running new simulations
# data.dir = str_replace( string = here(),
#                         pattern = "Code",
#                         replacement = "Results2" )
# results.dir = data.dir


setwd(code.dir)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")


# ~~ Temp only: check on sims in real time -------------------------

setwd(data.dir)
aggo = fread("agg.csv")

# check when the dataset was last modified to make sure we're working with correct version
file.info("aggo.csv")$mtime

# # for working with stitched directly
# s = fread( "stitched.csv")
# # check when the dataset was last modified to make sure we're working with correct version
# file.info("stitched.csv")$mtime
# 
# if ( "method.1" %in% names(s) ) s = s %>% select(-method.1)
# 
# s %>% group_by(scen.name, method) %>%
#   summarise(n(),
#             mean(Mu),
#             meanNA(Mhat),
#             mean(is.na(Mhat)))
# agg = make_agg_data(s)

agg = wrangle_agg_local(aggo)

table(agg$method.pretty)
table(agg$evil.selection)

# which key scen params have run so far?
library(tableone)
param.vars.short = c(
               "hack",
               "k.pub.nonaffirm",
               "prob.hacked",
               
               "eta",
               "SAS.type",
               
               "true.dist",
               "true.sei.expr",
               
               "muB",
               "t2a"
               )

CreateCatTable(vars = param.vars.short,
               data = agg)



# for after running sims
# # ~~ Get agg data -------------------------
# 
# # if only analyzing a single set of sims (no merging):
# setwd(data.dir)
# agg = fread( "agg.csv")
# # check when the dataset was last modified to make sure we're working with correct version
# file.info("agg.csv")$mtime
# 
# 
# dim(agg)  # will exceed number of scens because of multiple methods
# expect_equal( 90, nuni(agg$scen.name) )
# 
# # prettify variable names
# agg = wrangle_agg_local(agg)
# table(agg$method.pretty)
# 
# # look at number of actual sim reps
# table(agg$sim.reps.actual)
# 
# 
# # one-off for paper
# update_result_csv( name = "n scens",
#                    value = nuni(agg$scen.name),
#                    .results.dir = results.dir,
#                    .overleaf.dir = overleaf.dir.nums )


# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# *** BEST AND WORST PERFORMANCE ACROSS SCENS -------------------------

t = agg %>%
  filter( method %in% c("naive", "mbma-MhatB", "2psm", "beta-sm") ) %>%
  group_by(method) %>%
  summarise(BiasMin = min(MhatBias),
            BiasMd = median(MhatBias),
            BiasMax = max(MhatBias),
            
            AbsBiasMin = min(MhatAbsBias),
            AbsBiasMd = median(MhatAbsBias),
            AbsBiasMax = max(MhatAbsBias),
            
            MhatEstFail = median(MhatEstFail),
            
            #*express coverage as percent
            CoverMin = 100*min(MhatCover),
            CoverMd = 100*median(MhatCover),
            
            WidthMd = median(MhatWidth),
            WidthMin = min(MhatWidth),
            WidthMax = max(MhatWidth)) 

t

# go through each outcome column in table (e.g., BiasMin) and write
#  results for each method to csv
for ( .col in names(t)[ 2 : ncol(t) ] ) {
  
  update_result_csv( name = paste( t$method, .col, sep = " " ),
                     value = t[[.col]],
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir.nums )
}


# ******** RANKED PERFORMANCE TABLES -------------------------

# all scenarios
( t1.mn = make_winner_table(.agg = agg,
                       summarise.fun.name = "mean" ) )
( t1.worst = make_winner_table(.agg = agg,
                          summarise.fun.name = "worst10th" ) )


# scens with SWS
temp = agg %>% filter(prob.hacked == 1); dim(temp)
( t1.mn = make_winner_table(.agg = temp,
                            summarise.fun.name = "mean" ) )
( t1.worst = make_winner_table(.agg = temp,
                               summarise.fun.name = "worst10th" ) )


# skewed effects
temp = agg %>% filter(true.dist == "expo"); dim(temp)
( t1.mn = make_winner_table(.agg = temp,
                            summarise.fun.name = "mean" ) )
( t1.worst = make_winner_table(.agg = temp,
                               summarise.fun.name = "worst10th" ) )

# later (not run yet): evil.selection == 1 or == 0
temp = agg %>% filter(evil.selection == 1); dim(temp)
( t1.mn = make_winner_table(.agg = temp,
                            summarise.fun.name = "mean" ) )
( t1.worst = make_winner_table(.agg = temp,
                               summarise.fun.name = "worst10th" ) )

# ******** PLOTS (BIG AND NOT PRETTIFIED) -------------------------

Ynames = rev(MhatYNames)

# alternatively, run just a subset:
Ynames = c("MhatWidth", "MhatCover", "MhatBias", "MhatRMSE")

# to help decide which vars to include in plot:
param.vars.manip2


# in case you want to filter scens:
# full set for reference:
# c("naive", "gold-std", "maon", "2psm", "pcurve", "jeffreys-mcmc-pmean", 
#   "jeffreys-mcmc-pmed", "jeffreys-mcmc-max-lp-iterate", "jeffreys-sd", 
#   "jeffreys-var", "mle-sd", "csm-mle-sd", "mle-var", "2psm-csm-dataset", 
#   "prereg-naive", "ltn-mle-sd")
( all.methods = unique(agg$method) )
# toDrop = c("rtma-adj-MhatB-pmed",
#            "rtma-adj-muB-pmed",
#            "rtma-adj-MhatB-pmean",
#            'rtma-adj-muB-pmean')
toDrop = NULL
method.keepers = all.methods[ !all.methods %in% toDrop ]


# make facetted plotly

for ( .t2a in unique(agg$t2a) ) {
  
  # # test only
  # .hack = "affirm"
  # .t2a = 0.25
  
  cat( paste("\n\n -------- STARTING t2a=", .t2a) )
  
  aggp = agg %>% filter(method %in% method.keepers &
                          t2a == .t2a)
  # to label the plots
  prefix = paste( "2022-7-23 sims; ",
                  "t2a=", .t2a,
                  sep = "")
  
  
  # temporarily set wd
  # results.dir.temp = paste(results.dir,
  #                          "/Big unprettified plots/",
  #                          .Mu,
  #                          "/hack=",
  #                          .hack,
  #                          sep = "")
  
  results.dir.temp = paste(results.dir,
                           "/Big unprettified plots",
                           sep = "")
  
  
  # set facetting variables for plots
  aggp$tempFacetVar1 = paste( "muB=", aggp$muB,
                              sep = "")
  
  aggp$tempFacetVar2 = as.factor( paste( "eta=", aggp$eta,
                                         sep = "") )
  levels(aggp$tempFacetVar2) = c("eta=1", "eta=5", "eta=10")
  # force factor ordering
  table(aggp$tempFacetVar2)
  
  
  for ( Yname in Ynames) {
    
    # to run "manually"
    #Yname = "MhatBias"
    #Yname = "MhatCover"
    
    y.breaks = NULL
    if ( Yname == "MhatBias") y.breaks = round( seq(-0.6, 0.7, 0.1), 2)
    if ( Yname == "MhatWidth") y.breaks = seq(0, 2, 0.2)
    
    p  = quick_5var_agg_plot(.Xname = "k.pub.nonaffirm",
                             .Yname = Yname,
                             .colorVarName = "method",
                             .facetVar1Name = "tempFacetVar1",
                             .facetVar2Name = "tempFacetVar2",
                             .dat = aggp,
                             .ggtitle = prefix,
                             .y.breaks = y.breaks,
                             .writePlot = FALSE,
                             .results.dir = results.dir.temp)
    
    
    pl = ggplotly(p)
    
    # in filename, mark the most important plots with asterisk
    if ( Yname %in% c("MhatBias", "MhatCover", "MhatWidth") ){
      new.prefix = paste("*", prefix, sep = "")
    } else {
      new.prefix = prefix
    }
    
    # how to save a plotly as html
    # https://www.biostars.org/p/458325/
    setwd(results.dir.temp)
    string = paste(new.prefix, Yname, "plotly.html", sep="_")
    htmlwidgets::saveWidget(pl, string)
    
  }
  
}




# ******** PLOTS (SIMPLE AND PRETTY FOR MAIN TEXT) -------------------------

# for each hack type, arrange plots so each facet row is an outcome
( all.methods = unique(agg$method.pretty) )
toDrop = c("MBMA (true t2)", "MAN adjusted")  # methods to exclude from plots
( method.keepers = all.methods[ !is.na(all.methods) &
                                  !(all.methods %in% toDrop) ] )


# outcomes to show in main text figures
YnamesMain = c("MhatBias", "MhatCover", "MhatWidth")



pretty.results.dir = paste(results.dir, "/Prettified plots", sep = "")

t2a.levels = unique(agg$t2a)

# this dataset will be one full-page figure in main text or Supp depending on hack type
# by default, these write only to Overleaf dir

for ( i in 1:length(t2a.levels) ) {
  taui = sqrt(t2a.levels[i])
  sim_plot_multiple_outcomes( .t2a = t2a.levels[i],
                              .ggtitle = bquote( tau ~ " = " ~ .(taui) ),
                              .local.results.dir = pretty.results.dir )
}

 

