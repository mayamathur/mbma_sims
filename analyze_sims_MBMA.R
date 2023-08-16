
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
library(broom)

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

# proportion of full-factorial scens that were able to run
nuni(aggo$scen.name) / 12980


# prettify variable names
agg = wrangle_agg_local(aggo)


# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names()

# remove experimental methods
table(agg$method)
agg = agg %>% filter(!(method %in% c("mbma-MhatB-gamma", "maon-adj-MhatB")))


# for analyzing sims as they run:
# which key scen params have run so far?
CreateCatTable(vars = param.vars.manip,  # param.vars.manip is from init_var_names
               data = agg)




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

update_result_csv( name = "prop scens evilselect1",
                   value = round( 100*mean(agg$evil.selection == 1) ),
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

update_result_csv( name = "n scens large SEs",
                   value = nuni(agg$scen.name[ agg$true.sei.expr == "0.02 + rexp(n = 1, rate = 1)"]),
                   .results.dir = results.dir,
                   .overleaf.dir = overleaf.dir.nums )


update_result_csv( name = "n scens small muB",
                   value = nuni(agg$scen.name[ agg$muB == 0.1]),
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




# ******** WINNER TABLES -------------------------

# ~ Winner tables shown in paper ------------------------------

# can toggle output of fn below by changing the default arg of 
#  make_winner_table between display = "dataframe" (easy viewing)
#  and display = "xtable" (Overleaf)

# 1: all scenarios
make_both_winner_tables(.agg = agg)
# here, reason SM-step appears unbiased is that it's positively biased under evil.selection=0
#   but negatively biased under evil.selection=1


# 2: scenarios where our method IS correctly spec
make_both_winner_tables(.agg = agg %>% filter(evil.selection == 0) )

# 3: scenarios where our method is NOT correctly spec
make_both_winner_tables(.agg = agg %>% filter(evil.selection == 1) )

# 4: small k.pub.nonaffirm
make_both_winner_tables(.agg = agg %>% filter(k.pub.nonaffirm == 5) )

# 5: small mean bias
make_both_winner_tables(.agg = agg %>% filter(muB == 0.1) )
# c.f.: other values
# make_both_winner_tables(.agg = agg %>% filter(muB > 0.1) )


# other possibilities suggested by reviewers:
# - large within-study SE
# - muB = 0.1
# - Mu = 0 (null)


# ~ Winner tables with results discussed only in prose  --------------------------------------------

# 6: skewed effects
make_both_winner_tables(.agg = agg %>% filter(true.dist=="expo") )


# For Reviewer 3's questions
# 7: larger within-study SEs
make_both_winner_tables(.agg = agg %>% filter(true.sei.expr == "0.02 + rexp(n = 1, rate = 1)") )
# c.f.: smaller SEs
# make_both_winner_tables(.agg = agg %>% filter(true.sei.expr == "0.02 + rexp(n = 1, rate = 3)") )



# ~ Winner tables for own curiosity --------------------------------------------

# scens without vs. with p-hacking
make_both_winner_tables(.agg = agg %>% filter(prob.hacked == 0) )
make_both_winner_tables(.agg = agg %>% filter(prob.hacked == 1) )


# scens with SAS.type=carter
make_both_winner_tables(.agg = agg %>% filter(SAS.type == "carter") )




# REGRESSIONS -------------------------

# under what conditions are our method and SM-step NEGATIVELY biased?

.method = "mbma-MhatB"
( string = paste( "(MhatBias < 0) ~", paste(param.vars.manip, collapse = "+"), sep="" ) )
tidy( lm( eval(parse(text=string)),
          data = agg %>% filter(method == .method) ) )


.method = "2psm"
tidy( lm( eval(parse(text=string)),
          data = agg %>% filter(method == .method) ) )



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



