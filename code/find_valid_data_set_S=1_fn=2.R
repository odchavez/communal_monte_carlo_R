rm(list = ls())
library(gtools)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(foreach)
library(doParallel)
library(Matrix )
library(lpSolve)

tacc_wd = "/work/03330/ochavez/stampede2/communal_monte_carlo/communal_monte_carlo_R"
#if(file.exists(tacc_wd)){
    setwd(tacc_wd)    
#}

source("code/wasserstein_distance_functions.R")
source("code/barycenter.r")
source("code/communal_monte_carlo_functions.R")

options(max.print=999999)


#MVN mixture Gibbs Sampler
file_num       = 2
N              = 20000 #number of datapoints per file up to n below
experiment_num = 0

shard_num      = 1     #processor number
np             = 5000  #particle number
 

n              = 20000
K              = 10
d              = 2
scale          = 10



global_steps   = 1
gstep_com      = FALSE

#one processor



data_file_name = paste0("data/K=",K,"_n=",n ,"/d=",d,"_n=",n,"_file_num_",file_num,".csv")
dat            = read.csv(data_file_name)[1:floor(N), ]
priors_list    = get_default_priors(K, d, scale, np, shard_num)

posterior_list = do_communal_mc_MVN_mix_single_file(dat, 
                                                    global_steps, 
                                                    shard_num = 1,  
                                                    priors_list, 
                                                    experiment_num, 
                                                    gstep_com = FALSE,
                                                    file_num)

save_particles(posterior_list, shard_num, global_steps, K, d, n, np, experiment_num, file_num,type="S=1")
  

