rm(list = ls())
library(gtools)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(foreach)
library(doParallel)
library(Matrix )
library(lpSolve)

setwd("/work/03330/ochavez/stampede2/communal_monte_carlo/communal_monte_carlo_R")
source("code/wasserstein_distance_functions.R")
source("code/barycenter.r")
source("code/communal_monte_carlo_functions.R")

options(max.print=999999)


#MVN mixture Gibbs Sampler
shard_num      = 10    #processor number
np             = 200  #particle number
N              = 1e+06 #number of datapoints per file up to n below
gstep_com      = TRUE 


n              = 500000
K              = 10
d              = 2
scale          = 100

max_file_num   = 2

experiment_num = 0
global_steps   = 1





exp_results_fout_name = paste0("results/log_lik_experiment_num=",experiment_num,".RData")
test_dat = read.csv(paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_10.csv"))[500,]

#embarasingly parrallel
gstep_com = FALSE
count     = 1
log_lik = c()
for(fn in 1:max_file_num){
  
  file_num       = fn
  
  data_file_name = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
  dat            = read.csv(data_file_name)[1:N, ]
  priors_list    = get_default_priors(K, d, scale, np, shard_num)
  posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                      global_steps, 
                                                      shard_num,  
                                                      priors_list, 
                                                      experiment_num, 
                                                      gstep_com)
  save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num, type="emb_par")
  
  #fraction = count/(max_file_num)
  #log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))
  #count = count + 1
  
}
#all_loglik = cbind(log_lik, all_loglik)

#base comparison - one processor

np             = np#*shard_num
count     = 1
log_lik = c()
for(fn in 1:max_file_num){
  
  file_num       = fn
  
  data_file_name = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
  dat            = read.csv(data_file_name)[1:floor(N), ]
  priors_list    = get_default_priors(K, d, scale, np, shard_num)
  posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                      global_steps = 1, 
                                                      shard_num = 1,  
                                                      priors_list, 
                                                      experiment_num, 
                                                      gstep_com = FALSE)
  save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num,type="S=1")
  
  #fraction = count/(max_file_num)
  #log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))
  #count = count + 1
  
}
