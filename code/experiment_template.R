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

options(max.print=999999)

#MVN mixture Gibbs Sampler
experiment_num = 1
shard_num      = 4
np             = 10
gstep_com      = TRUE

n              = 10000
K              = 20
d              = 2
scale          = 100

max_global_step = 10
max_file_num    = 10
test_dat = read.csv("data/K=20/d=2_n=10000_file_num_96.csv")
log_lik_matrix = matrix(NA, nrow = max_global_step, ncol = 3)
colnames(log_lik_matrix) = c("global_steps", "mean", "sd")
for(g in 1:100){#max_global_step){
  
    global_steps = 2#g
    log_lik = c()
    for(fn in 1:max_file_num){
        print(paste("global step = ",g," file number = ",fn))
        file_num       = fn
      
        data_file_name = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
        dat            = read.csv(data_file_name)[1:1000, ]
        priors_list    = get_default_priors(K, d, scale, np, shard_num)
        posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                            global_steps, 
                                                            shard_num,  
                                                            priors_list, 
                                                            experiment_num, 
                                                            gstep_com)
        #save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num)
        
        #fraction = ((g-1)+fn)/(max_global_step*max_file_num)
        #log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))

        
        
    }
    log_lik_matrix[g,"global_steps"] = g
    log_lik_matrix[g,"mean"]         = mean(log_lik)
    log_lik_matrix[g,"sd"]           = sd(log_lik)/sqrt(length(log_lik))
}
plot(log_lik_matrix[,"global_steps"], log_lik_matrix[,"mean"])



