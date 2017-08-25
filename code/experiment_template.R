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
shard_num      = 6
np             = 100
gstep_com      = TRUE 
N              = 1000

n              = 1000
K              = 10
d              = 2
scale          = 100

max_global_step = 10
max_file_num    = 10
files_list = c(1:10)
global_steps_list = c(1:6)
test_dat = read.csv(paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_10.csv"))[500,]
log_lik_matrix = matrix(NA, nrow = max_global_step, ncol = 3)
colnames(log_lik_matrix) = c("global_steps", "mean", "sd")
count = 1
for(g in 1:max_global_step){
  
    global_steps = g
    log_lik = c()
    
    for(fn in 1:max_file_num){
        print(paste("global step = ",g," file number = ",fn))
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
        save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num)
        
        fraction = count/(max_global_step*max_file_num)
        log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))
        count = count + 1
    }
    
    log_lik_matrix[g,"global_steps"] = g
    log_lik_matrix[g,"mean"]         = median(log_lik)
    log_lik_matrix[g,"sd"]           = sd(log_lik)/sqrt(length(log_lik))
    if(g == 1){
      all_loglik = as.data.frame(log_lik)
    }else{
      all_loglik = cbind(all_loglik, log_lik)
    }
    
}

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
  #save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num)
  
  fraction = count/(max_file_num)
  log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))
  count = count + 1
  
}
all_loglik = cbind(log_lik, all_loglik)

#base comparison - one processor

np             = np*shard_num
count     = 1
log_lik = c()
for(fn in 1:max_file_num){
  
  file_num       = fn
  
  data_file_name = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
  dat            = read.csv(data_file_name)[1:N, ]
  priors_list    = get_default_priors(K, d, scale, np, shard_num)
  posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                      global_steps = 1, 
                                                      shard_num = 1,  
                                                      priors_list, 
                                                      experiment_num, 
                                                      gstep_com = FALSE)
  #save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num)
  
  fraction = count/(max_file_num)
  log_lik = c(log_lik, experiment_log_lik(test_dat, posterior_list, fraction))
  count = count + 1
  
}


all_loglik = cbind(log_lik, all_loglik)
plot(log_lik_matrix[,"global_steps"], log_lik_matrix[,"mean"])
boxplot(all_loglik, ylim = c(-10,0))
bs_means = bootstrap_columns(all_loglik, nrep = 500, statistic = "mean")
    

