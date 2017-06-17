library(gtools)
library(MASS)
library(MCMCpack)

setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")

plot_means = function(final_params){
  P = length(final_params)
  for(i in 1:P){
    points(final_params[[i]]$mean)
  }
}

sample_particles = function(log_lik_wts, candidate_particles){
  np = length(log_lik_wts)
  select = sample(np, np, prob = exp(log_lik_wts), replace = TRUE)
  particles_out = list()
  for(p in 1:np){
    cp = select[p]
    particles_out[[p]] = candidate_particles[[cp]]
    #points(particles_out[[p]]$mean, col = 'red')
  }
  return(particles_out)
}

#lik_weight = function(data_line, candidate_particles)
  
log_lik_mvn = function(mu, Sig, dat){
  
  Sig_Inv = solve(Sig)
  k = length(dat)
  X = as.numeric(dat)
  
  A = (-1/2) * k *log(2*pi) 
  B = (-1/2) * log(det(Sig)) 
  C = -(1/2) * (X - mu)%*% Sig_Inv %*% (X - mu)
  output = A + B + C
  return(output)
}

log_lik_mvn_mix = function(param_sample, dat){
  
  w = param_sample$prob
  M = param_sample$mean
  S = param_sample$sig
  k = ncol(w)
  comps = rep(NA, k)
  for(i in 1:k){
    mu  = M[i,]
    Sig =  matrix(S[i,], ncol = 2)
    comps[i] = log(w[i]) + log_lik_mvn(mu, Sig, dat)
  }
  #print(comps)
  output = log(sum(exp(comps)))
  return(output)
}


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
  nu             = priors$nu[z_i]
  post_mean = (prior_obs_size * prior_mu + n*data_line)/(prior_obs_size + n)
  #print(post_mean)
  d = length(prior_mu)
  #print(class(z_i))
  #print(priors$sig[z_i,])
  #print(matrix(priors$sig[z_i,], ncol = d))
  post_cov = matrix(priors$sig[z_i,], ncol = d)/(nu * (nu - d + 1))
  output = mvrnorm(n = 1, post_mean, post_cov)
  
  return(output)
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
  #print(prob_z_n)
  return(prob_z_n)
}

particle_filter_MVN_iter = function(data_line, priors){
  
  #weights
  prob_z_n = sample_dirichlet(priors)
  priors$prob = prob_z_n
  #counts and class assignment
  z_i      = sample_z_i(priors)
  priors$z_count[z_i] = priors$z_count[z_i] + 1
  priors$nu[z_i] = priors$nu[z_i] + 1
  #mean
  mu_k = sample_mu_k(data_line, priors, z_i)
  priors$mean[z_i,] = mu_k
  
  Sig_k = sample_Sig_k(data_line, priors, z_i)
  priors$sig[z_i,] = Sig_k
  
  return(priors)
}


particle_filter_MVN = function(file_name, priors_list, np, data_size){
  
  stop = FALSE
  f = file(file_name, "r")
  count = 0
  while(!stop) {
    
    next_line = readLines(f, n = 1)
    if(length(next_line) != 0){
      count = count+1
      print(paste(100*count/data_size, "% complete"))
      data_line = as.numeric(strsplit(next_line, ",")[[1]])
      #print("doing an iteration on")
      #print(paste("data_line = ", toString(data_line)))
      ## Insert some if statement logic here
      candidate_particles = list()
      log_lik_wts = rep(0, np)
      for(p in 1:np){
          candidate_particles[[p]] = particle_filter_MVN_iter(data_line, priors_list[[p]])
          log_lik_wts[p]           = log_lik_mvn_mix(candidate_particles[[p]], data_line)
      }
      #print(log_lik_wts)
      posterior = sample_particles(log_lik_wts, candidate_particles)
      priors_list = posterior
    }else{
      stop = TRUE
      close(f)
    }
  }
  return(priors_list)
}

#MVN mixture Gibbs Sampler
data_size = 3000
K         = 5
d         = 2
shard_num = 1
scale     = 5
np        = 100
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
nu      = matrix(rep(2, K), ncol = K)
priors = list(prob = prob, 
              z_count = z_count, 
              alpha_k = alpha_k, 
              mean = mean,
              mean_size_init = alpha_k,
              sig = sig, 
              K = K,
              nu = nu)

priors_list = list(priors)[rep(1,np)]
plot(read.csv(data_file))
final_params = particle_filter_MVN(data_file, priors_list, np, data_size)
plot_means(final_params)
  
