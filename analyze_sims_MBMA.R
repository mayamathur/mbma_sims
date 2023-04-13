
# Note: This script shares a stats_for_paper.csv with the applied examples, so if you wr(), 
#  you'll also need to rerun the applied examples.


# PRELIMINARIES ----------------------------------------------------

#  rm(list=ls())

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
library(tableone)
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
                        replacement = "Results/2023-04-09 full sims for manuscript (RSM_1)" )


# temp while sims are ongoing
results.dir = data.dir
# results.dir = str_replace( string = here(),
#                            pattern = "Code",
#                            replacement = "Results/*2022-7-23 More scens for manuscript; add mbma-MhatB-true-t2" )



overleaf.dir.figs = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/figures/sims"

# for stats_for_paper.csv
# same dir as for applied examples so that they'll write to single file
overleaf.dir.nums = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/R_objects"


setwd(code.dir)
source("helper_MBMA.R")
source("analyze_sims_helper_MBMA.R")


# ~~ Get agg data -------------------------

setwd(data.dir)
aggo = fread("agg.csv")

# check when the dataset was last modified to make sure we're working with correct version
file.info( paste(data.dir, "agg.csv", sep="/") )$mtime

nrow(aggo)  # will exceed number of scens because of multiple methods
expect_equal( 12980, nuni(aggo$scen.name) )

# proportion done
nuni(aggo$scen.name) / 12980


# prettify variable names
agg = wrangle_agg_local(aggo)


# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()

#@temp
table(agg$method)
agg = agg %>% filter(!(method %in% c("mbma-MhatB-gamma", "maon-adj-MhatB")))

# make new var
agg$MhatEstConverge = 1 - agg$MhatEstFail

# for analyzing sims as they run:
# which key scen params have run so far?
CreateCatTable(vars = param.vars.manip,  # param.vars.manip is from init_var_names
               data = agg)





# ~~ List variable names -------------------------




# ONE-OFF STATS FOR PAPER  -------------------------

# scenario counts
update_result_csv( name = "n scens",
                   value = nuni(agg$scen.name),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )


update_result_csv( name = "n scens evilselect0",
                   value = nuni(agg$scen.name[ agg$evil.selection == 0]),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "n scens evilselect1",
                   value = nuni(agg$scen.name[ agg$evil.selection == 1]),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "n scens smallk",
                   value = nuni(agg$scen.name[ agg$k.pub.nonaffirm == 5]),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = "n scens skewed",
                   value = nuni(agg$scen.name[ agg$true.dist == "expo"]),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

# convergence
update_result_csv( name = "beta-sm convergence",
                   value = round( 100*median(agg$MhatEstConverge[agg$method == "beta-sm"]), 0 ),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

# sample sizes
update_result_csv( name = paste( "sancheck.dp.k", c("Q1", "median", "Q3") ),
                   value = round( summary(agg$sancheck.dp.k)[2:4], 0 ),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )

update_result_csv( name = paste( "sancheck.dp.k.affirm", c("Q1", "median", "Q3") ),
                   value = round( summary(agg$sancheck.dp.k.affirm)[2:4], 0 ),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )



# BEST AND WORST PERFORMANCE ACROSS SCENS -------------------------

#@ maybe don't put in paper since redundant with winner tables below?

t = agg %>%
  filter( method %in% c("naive", "mbma-MhatB", "2psm", "beta-sm") ) %>%
  group_by(method) %>%
  summarise(BiasMin = min(MhatBias),
            BiasMd = median(MhatBias),
            BiasMax = max(MhatBias),
            
            AbsBiasMin = min(MhatAbsBias),
            AbsBiasMd = median(MhatAbsBias),
            AbsBiasMax = max(MhatAbsBias),
            
            RMSEMin = min(MhatAbsBias),
            RMSEMd = median(MhatAbsBias),
            RMSEMax = max(MhatAbsBias),
            
            #*express coverage as percent
            CoverMin = 100*min(MhatCover),
            CoverMd = 100*median(MhatCover),
            
            WidthMd = median(MhatWidth),
            WidthMin = min(MhatWidth),
            WidthMax = max(MhatWidth),
            
            MhatEstFail = median(MhatEstFail) ) 

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
# here, reason SM-step appears unbiased is that it's positively biased under evil.selection=0
#   but negatively biased under evil.selection=1
make_both_winner_tables(.agg = agg)


# # diagnosis: scens with more confounding
# make_both_winner_tables(.agg = agg %>% filter(muB == 0.5) )


# scenarios where our method IS correctly spec
make_both_winner_tables(.agg = agg %>% filter(evil.selection == 0) )

# scenarios where our method is NOT correctly spec
make_both_winner_tables(.agg = agg %>% filter(evil.selection == 1) )

# small k.pub.nonaffirm
make_both_winner_tables(.agg = agg %>% filter(k.pub.nonaffirm == 5) )

# skewed effects
make_both_winner_tables(.agg = agg %>% filter(true.dist=="expo") )


# other possibilities suggested by reviewers:
# - large within-study SE
# - muB = 0.1
# - Mu = 0 (null)




# CARTER PUB BIAS PLOT -------------------------


dp = expand.grid(pval = seq(0.0001, 0.15, by = 0.001) )

# parameters for carter_censor taken from sim_one_study_set
dp = dp %>% rowwise() %>%
  mutate( pub.prob = carter_censor( pObs = pval,
                                    direction = 1, 
                                    posSign_NS_baseRate = 0.3,
                                    negSign_NS_baseRate = 0.05,
                                    counterSig_rate = 0.50 ) )


# key cutpoints
x.breaks = seq(0, 0.15, 0.025)



p = ggplot( data = dp, aes(
  x = pval, 
  y = pub.prob) ) +
  geom_line() +
  
  # use only some values
  #scale_x_log10( breaks = c(500, 1000) ) +
  
  xlab("Two-tailed p-value") +
  scale_x_continuous( breaks = x.breaks ) +
  
  
  ylab("Publication probability") +
  scale_y_continuous( breaks = seq(0, 1, 0.1) ) +
  
  theme(legend.position = "bottom") +
  
  # base_size controls all text sizes; default is 11
  # https://ggplot2.tidyverse.org/reference/ggtheme.html
  theme_bw(base_size = 16) +
  theme( text = element_text(face = "bold"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank() ) 

p

if ( overwrite.res == TRUE ) {
  my_ggsave( name = "carter_pubbias_plot.pdf",
             .plot = p,
             .width = 8,
             .height = 6,
             .results.dir = results.dir,
             .overleaf.dir = overleaf.dir.figs )
} else {
  message("\n\nNot writing the plot to local dir or Overleaf because overwrite.res = FALSE")
}









###### TO REMOVE?

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



