plot_means = function(final_params, color){
  P = length(final_params)
  for(i in 1:P){
    points(final_params[[i]]$mean, col = color)
  }
}

sample_particles = function(log_lik_wts, candidate_particles){
  np = length(log_lik_wts)
  print("in sample_particles")
  print( exp(log_lik_wts))
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
  nu = priors$nu[z_i]
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

sample_z_i = function(priors, data_line){
  #print("enter sample_z_i")
  K = priors$K
  prob_z_n = rep(0,K)
  w = priors$prob
  d = ncol(priors$mean)
  for(k in 1:K){
    mu = priors$mean[k,]
    Sig = matrix(priors$sig[k,], ncol = d)
    #print(mu)
    #print(Sig)
    prob_z_n[k] = w[k]*dmvnorm(data_line, mu, Sig, log=FALSE)
  }
  
  #prob_z_n = priors$prob
  print("in sample_z_i")
  print(prob_z_n)
  z_i = sample(1:K, 1, prob = prob_z_n)
  #print("exit sample_z_i")
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
  #print("entered(particle_filter_MVN_iter)")
  #weights
  prob_z_n = sample_dirichlet(priors)
  priors$prob = prob_z_n
  #print("2 ran")
  #counts and class assignment
  z_i      = sample_z_i(priors, data_line)
  priors$z_count[z_i] = priors$z_count[z_i] + 1
  priors$nu[z_i] = priors$nu[z_i] + 1
  #print("3 ran")
  #mean
  mu_k = sample_mu_k(data_line, priors, z_i)
  priors$mean[z_i,] = mu_k
  #print("4 ran")
  
  Sig_k = sample_Sig_k(data_line, priors, z_i)
  priors$sig[z_i,] = Sig_k
  #print("5 ran")
  return(priors)
}


particle_filter_MVN = function(file_name, priors_list, np, data_size, quit_after_n, gs, sn, experiment_num, global_steps){
  
  stop = FALSE
  f = file(file_name, "r")
  count = 0
  while(!stop){
    
    next_line = readLines(f, n = 1)
    #print(paste("next_line = ", next_line))
    if(length(next_line) != 0 & count < quit_after_n){
      
      
      
      #print(paste(100*count/data_size, "% complete"))
      
      #sink("logfiles/log.txt", append=TRUE)
      #cat()
      
      cat(paste("log_gs=",gs, "out of _",global_steps, " ", 
                100*count/min(quit_after_n,data_size), 
                "% complete", "\n"), 
          file=paste0("logfiles/experiment_num=",
                      experiment_num,"_shard=",sn,".txt"), 
          append=TRUE)
      
      data_line = as.numeric(strsplit(next_line, ",")[[1]])
      if(is.na(data_line)){next}
      count = count+1
      print(count)
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
      #if(count%%500 == 0){
      #  plot_means(priors_list)  
      #}
      
    }else{
      stop = TRUE
      close(f)
    }
  }
  return(priors_list)
}
get_default_priors = function(K, d, scale, np, shard_num){
  #priors
  prob      = matrix(rep(1/K, K), ncol = K, nrow = 1)
  mean      = t(matrix(rep(0, K*d) , nrow = d))#matrix(rep(0, d*K), ncol = d, nrow = K)
  sig       = t(matrix(as.vector(diag(d))*scale, ncol = K, nrow = d*d))
  #params
  z_count = matrix(rep(0, K), ncol = K)
  alpha_k = matrix(rep(1, K), ncol = K)
  nu      = matrix(rep(2, K), ncol = K)
  priors  = list(prob = prob, 
                 z_count = z_count, 
                 alpha_k = alpha_k, 
                 mean = mean,
                 mean_size_init = alpha_k,
                 sig = sig, 
                 K = K,
                 nu = nu)
  
  priors_list = list(priors)[rep(1,np)]
  output = list(priors_list)[rep(1,shard_num)]
  return(output)
}

get_new_priors = function(final_params, shard_num, np){
  print("in get_new_priors")
  master_list   = c()
  output_priors = list()
  for(i in 1:shard_num){
    master_list = c(master_list, final_params[[i]])
  }
  distMat = get_distMat(final_params)
  print("distMat")
  print(distMat)
  colMeasure = get_colMeasure(final_params)
  print("colMeasure")
  print(colMeasure)
  con_rel_sol = Barycenter_measure(colMeasure, distMat, maxIter = 20, lambda = 0.1)
  print(con_rel_sol)
  for(i in 1:shard_num){

    index = sample(np*shard_num, np, replace = TRUE, prob = con_rel_sol)
    temp = list()
    for(ind in 1:np){
      temp[[ind]] = master_list[[index[ind]]] 
    }
    output_priors[[i]] = temp
    
  }
  return(output_priors)
}

save_particles = function(particles, 
                          shard_num, 
                          global_steps, 
                          K, d, n, np, 
                          experiment_num, 
                          quit_after_n){
  
  file_name = paste('particles/',     
                    'K=',   toString(K),
                    'd=',   toString(d),
                    '_n=',  toString(min(n, quit_after_n)*shard_num*global_steps),
                    '_pn=', toString(np),
                    '_sn=', toString(shard_num),
                    '_gs=', toString(global_steps),
                    '_experiment_number=',      toString(experiment_num) ,
                    '.RData',sep = "")
  
  save(particles, file = file_name)
}

omars_file_name = function(K, d, n, file_num){
  output = paste('data/K=',    toString(K),
                 '/d=',        toString(d),
                 '_n=',        toString(n),
                 '_file_num_', toString(file_num),
                 '.csv',sep = "")
  return(output)
}

#open_log_file = function(){
#  file_name = paste('logfiles/',     toString(K),
#                    'd=',             toString(d),
#                    '_n=',            toString(n*shard_num*global_steps),
#                    '_pn=',     toString(np),
#                    '_sn=',    toString(shard_num),
#                    '_gs=', toString(global_steps),
#                    '_time=',      toString(as.numeric(Sys.time())) ,
#                    '.RData',sep = "")
#  writeLines(c(""), paste("logfiles/log.txt"))
#  
#}

do_communal_mc_MVN_mix = function(global_steps, shard_num, K, d, n, priors_list, quit_after_n, experiment_num, gstep_com){
  for(gs in 1:global_steps){
    
    
    
    cl <- makeCluster(shard_num)
    registerDoParallel(cl)
    final_params = foreach(sn = 1:shard_num, 
                           .packages = c("MCMCpack", "mvtnorm")) %dopar%{
      source("code/communal_monte_carlo_functions.R")
                             
      writeLines(c(""), paste0("logfiles/experiment_num=",experiment_num,"_shard=",sn,".txt"))
                             
      #sn = 1                       
      file_num  = sn + (gs - 1)*shard_num
      data_file = omars_file_name(K, d, n, file_num)
      particle_filter_MVN(data_file, priors_list[[sn]], length(priors_list[[sn]]), n, quit_after_n, gs, sn,experiment_num, global_steps)
      
      
      
      #particle_filter_MVN(data_file, priors_list, np, n)
    }
    stopCluster(cl)
    
    if(gstep_com == FALSE){
      priors_list = final_params
    }else{
      print(paste("performing global step", gs))
      priors_list = get_new_priors(final_params, shard_num, np)
      print(paste("global step", gs, "complete"))
    }
    
  }
  return(priors_list)
}

experiment_log_lik = function(test_dat, particles, EXP){
  print(paste0("analyzing experiment ",EXP))
  count = 1
  N = nrow(test_dat)
  S = length(particles)
  P = length(particles[[1]])
  experiment_log_lik_vec = rep(NA, N*S*P)
  for(n in 1:N){
    if(100*N/n %% 10 == 0){
      print(paste0("analyzing experiment ",EXP,".....", 100*n/N,"% complete"))
    }
    for(s in 1:S){
      for(p in 1:P){
        experiment_log_lik_vec[count] = log_lik_mvn_mix(particles[[s]][[p]], test_dat[n,])
        count = count + 1
      }
    }
  }
  return(experiment_log_lik_vec)
}