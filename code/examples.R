library(gtools)
library(MASS)
library(MCMCpack)

setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")

sample_particles = function(log_wts, candidate_particles)

lik_weight = function(data_line, candidate_particles)

sample_Sig_k = function(data_line, priors, z_i){#nu, Sig_prior, N, Y, mu, prior_mu, z_vec, K){
  
  n = 1
    
  x_bar = data_line
  p_mu = priors$mean[ z_i,]
  d = length(p_mu)
  nu = priors$z_count[z_i]
  T     = (matrix(priors$sig[z_i,], ncol = d))*nu

  M     = (nu * n)/(nu + n) * 
    as.matrix(p_mu-x_bar, ncol = 1) %*% t(as.matrix(p_mu - x_bar))	

  T_n    = T + M
  post_COV = T_n/(nu * (nu - d + 1))
  post_pred_cov = post_COV*(nu + 1)
  Sig = solve(riwish(v = n+nu, S = T_n))
  output = as.vector(Sig)
  return(output)
}


sample_mu_k = function(data_line, priors, z_i){
  
  prior_obs_size = priors$z_count[z_i] + priors$mean_size_init[z_i]
  prior_mu       = priors$mean[z_i,]
  n              = 1
  post_mean = (prior_obs_size * prior_mu + n*data_line)/(prior_obs_size + 1)
  print(post_mean)
  d = length(prior_mu)
  #print(class(z_i))
  #print(priors$sig[z_i,])
  #print(matrix(priors$sig[z_i,], ncol = d))
  output = mvrnorm(n = 1, post_mean, matrix(priors$sig[z_i,], ncol = d))
  return(post_mean)
}

sample_z_i = function(priors){
  K = priors$K
  prob_z_n = priors$prob
  z_i = sample(1:K, 1, prob = prob_z_n)
  return(z_i)
}

sample_dirichlet = function(priors){
  n = sum(priors$z_count)
  alpha_vec = (priors$z_count + priors$alpha_k)#/(n - 1 + sum(priors$alpha_k))
  prob_z_n = rdirichlet(1, alpha_vec)
  print(prob_z_n)
  return(prob_z_n)
}

particle_filter_MVN_iter = function(data_line, priors){
  
  #weights
  prob_z_n = sample_dirichlet(priors)
  priors$prob = prob_z_n
  #counts and class assignment
  z_i      = sample_z_i(priors)
  priors$z_count[z_i] = priors$z_count[z_i] + 1
  #mean
  mu_k = sample_mu_k(data_line, priors, z_i)
  priors$mean[z_i,] = mu_k
  
  Sig_k = sample_Sig_k(data_line, priors, z_i)
  priors$sig[z_i,] = Sig_k
  
  return(priors)
}


particle_filter_MVN = function(file_name, priors_list, np){
  
  stop = FALSE
  f = file(file_name, "r")
  while(!stop) {
    
    next_line = readLines(f, n = 1)
    if(length(next_line) != 0){
      data_line = as.numeric(strsplit(next_line, ",")[[1]])
      print("doing an iteration on")
      print(paste("data_line = ", toString(data_line)))
      ## Insert some if statement logic here
      candidate_particles = list()
      for(p in 1:np){
          candidate_particles[[p]] = particle_filter_MVN_iter(data_line, priors_list[[p]])
      }
      #log_wts = log_lik_weight(data_line, candidate_particles)
      #posterior = sample_particles(log_wts, candidate_particles)
      priors_list = candidate_particles#posterior
    }else{
      stop = TRUE
      close(f)
    }
  }
  return(priors)
}

#MVN mixture Gibbs Sampler
data_size = 3000
K         = 5
d         = 2
shard_num = 1
scale     = 5
np        = 2000
data_file = paste('data/MVN_train_data_K=',
                  toString(K),'_data_size=',
                  toString(data_size),'.csv',
                  sep = "")

#priors
prob      = matrix(rep(1/K, K), ncol = K, nrow = 1)
mean      = matrix(rep(0, d*K), ncol = d, nrow = K)
sig       = t(matrix(as.vector(diag(d))*scale, ncol = K, nrow = d*d))
#params
z_count         = matrix(rep(0, K), ncol = K)
alpha_k = matrix(rep(1, K), ncol = K)

priors = list(prob = prob, 
              z_count = z_count, 
              alpha_k = alpha_k, 
              mean = mean,
              mean_size_init = alpha_k,
              sig = sig, 
              K = K)

priors_list = list(priors)[rep(1,np)]

final_params = particle_filter_MVN(data_file, priors_list, np)

plot(read.csv(data_file))  
