rm(list = ls())
library(gtools)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(foreach)
library(doParallel)

#setwd("/communal_monte_carlo_R/communal_monte_carlo_R")
source("code/wasserstein_distance_functions.R")
source("code/barycenter.r")
source("code/communal_monte_carlo_functions.R")


#MVN mixture Gibbs Sampler
shard_num      = 4
global_steps   = 1
n              = 10000
quit_after_n   = 10000
K              = 20
d              = 2
scale          = 100
np             = 1000
experiment_num = 1
gstep_com      = FALSE



priors_list    = get_default_priors(K, d, scale, np, shard_num)
posterior_list = do_communal_mc_MVN_mix(global_steps, shard_num, K, d, n, priors_list, quit_after_n, experiment_num, gstep_com)

save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, quit_after_n)


plot(read.csv(paste('data/K=',    toString(K),
                    '/d=',        toString(d),
                    '_n=',        toString(n),
                    '_file_num_', toString(1),
                    '.csv',sep = "")))
for(i in 1:shard_num){plot_means(posterior_list[[i]], i+1)}



