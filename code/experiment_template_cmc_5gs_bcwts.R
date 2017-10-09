rm(list = ls())
library(gtools)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(foreach)
library(doParallel)
library(Matrix )
library(lpSolve)

#setwd("/work/03330/ochavez/stampede2/communal_monte_carlo/communal_monte_carlo_R")
source("code/wasserstein_distance_functions.R")
source("code/barycenter.r")
source("code/communal_monte_carlo_functions.R")

options(max.print=999999)

#MVN mixture Gibbs Sampler
experiment_num = 6
shard_num      = 10    #processor number
np             = 200  #particle number
N              = 1e+06 #number of datapoints per file up to n below
gstep_com      = TRUE 


n              = 1e+06
K              = 10
d              = 2
scale          = 100

max_file_num    = 10

global_steps_list = c(5)
lgs               = length(global_steps_list)

for(g in 1:lgs){
  
    #global_steps = g
    log_lik = c()
    
    for(fn in 1:max_file_num){
        print(paste("global step = ",global_steps_list[g]," file number = ",fn))
        file_num       = fn
      
        data_file_name = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
        dat            = read.csv(data_file_name)[1:N, ]
        priors_list    = get_default_priors(K, d, scale, np, shard_num)
        posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                            global_steps_list[g], 
                                                            shard_num,  
                                                            priors_list, 
                                                            experiment_num, 
                                                            gstep_com)
        save_particles(posterior_list, shard_num, global_steps_list[g], K, d, n, np, experiment_num, file_num)
    }
}
