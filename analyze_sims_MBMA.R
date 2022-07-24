
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
#@ not all fns respect this setting
overwrite.res = FALSE


# ~~ Set directories -------------------------
code.dir = here()

data.dir = str_replace( string = here(),
                        pattern = "Code",
                        replacement = "Results" )



results.dir = str_replace( string = here(),
                           pattern = "Code",
                           replacement = "Results" )


overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/figures/sims"

# # alternative for running new simulations
# data.dir = str_replace( string = here(),
#                         pattern = "Code \\(git\\)",
#                         replacement = "Simulation results" )
# 
# results.dir = data.dir


setwd(code.dir)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")


# ~~ Get agg data -------------------------

# if only analyzing a single set of sims (no merging):
setwd(data.dir)
agg = fread( "agg.csv")
# check when the dataset was last modified to make sure we're working with correct version
file.info("agg.csv")$mtime


dim(agg)  # will exceed number of scens because of multiple methods
expect_equal( 90, nuni(agg$scen.name) )

# prettify variable names
agg = wrangle_agg_local(agg)

# look at number of actual sim reps
table(agg$sim.reps.actual)


# ~~ List variable names -------------------------

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()


# BEST AND WORST PERFORMANCE IN CORRECTLY-SPECIFIED SCENS -------------------------

# scens where RTMA is correctly specified, but all other methods aren't
as.data.frame( agg %>% 
                 group_by(method) %>%
                 summarise( min(MhatBias),
                            median(MhatBias),
                            max(MhatBias),
                            min(MhatCover),
                            median(MhatCover),
                            max(MhatWidth)) )


# ******** PLOTS (BIG AND NOT PRETTIFIED) -------------------------

Ynames = rev(MhatYNames)

# alternatively, run just a subset:
Ynames = c("MhatWidth", "MhatCover", "MhatBias")

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

#bm: Remaining 6 ones keep failing -- is it the missing MhatB again?
#  Then re-stitch.
#  then edit sim_plot_multiple_outcomes below :)

# for each hack type, arrange plots so each facet row is an outcome
( all.methods = unique(agg$method.pretty) )
( method.keepers = all.methods[ !is.na(all.methods) &
                                  all.methods != "Gold standard"] )


# outcomes to show in main text figures
YnamesMain = c("MhatBias", "MhatCover", "MhatWidth")

# outcomes to show in supplement figures
YnamesSupp = c("MhatBias", "MhatCover", "MhatWidth")

# this dataset will be one full-page figure in main text or Supp depending on hack type
# by default, these write only to Overleaf dir
pl1 = sim_plot_multiple_outcomes(.hack = "favor-best-affirm-wch",
                                 .ggtitle = bquote( "SWS favors best affirmative; stringent SAS;" ~ mu ~ "= 0.5" ),
                                 .local.results.dir = results.dir )


pl2 = sim_plot_multiple_outcomes(.hack = "affirm",
                                 .ggtitle = bquote( "SWS favors first affirmative; stringent SAS; " ~ mu ~ "= 0.5" ),
                                 .local.results.dir = results.dir)



pl3 = sim_plot_multiple_outcomes(.hack = "affirm2",
                                 .ggtitle = bquote( "SWS favors first affirmative; no SAS; " ~ mu ~ "= 0.5" ),
                                 .local.results.dir = results.dir)





# 2022-4-4: EFFECT OF SCEN PARAMS ON DATASETS -------------------------

param.vars.manip2 = drop_vec_elements(param.vars.manip, "method")

t = agg %>% group_by_at( param.vars.manip2 ) %>%
  # keep only scens in main text
  filter(Mu == 0.5 & true.sei.expr == "0.02 + rexp(n = 1, rate = 3)") %>%
  select( all_of( contains("sancheck") ) ) %>%
  mutate_if( is.numeric, function(x) round(x, 2) )

setwd(results.dir)
setwd("Sanity checks")
write.xlsx( as.data.frame(t), "table_sanchecks.xlsx")

