rm(list = ls())
library(gtools)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(foreach)
library(doParallel)

setwd("/Users/o/Google Drive/school/Williamson Research/communal_monte_carlo_R")
source("code/wasserstein_distance_functions.R")
source("code/barycenter.R")
source("code/communal_monte_carlo_functions.R")


#MVN mixture Gibbs Sampler
shard_num    = 16
global_steps = 6
n            = 10000
K            = 20
d            = 2
scale        = 100
np           = 1000



priors_list    = get_default_priors(K, d, scale, np, shard_num)
posterior_list = do_communal_mc_MVN_mix(global_steps, shard_num, K, d, n, priors_list)

save_particles(posterior_list, shard_num, global_steps, K, d, n, np)


plot(read.csv(paste('data/K=',    toString(K),
                    '/d=',        toString(d),
                    '_n=',        toString(n),
                    '_file_num_', toString(1),
                    '.csv',sep = "")))
for(i in 1:shard_num){plot_means(posterior_list[[i]], i+1)}



