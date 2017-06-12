
setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")

sample_dirichlet = function(priors){
  
}

particle_filter_MVN_iter = function(data_line, priors){
  
  sample_dirichlet(priors)
  return()
}


particle_filter_MVN = function(file_name, priors){
  
  stop = FALSE
  f = file(file_name, "r")
  while(!stop) {
    
    next_line = readLines(f, n = 1)
    if(length(next_line) != 0){
      data_line = as.numeric(strsplit(next_line, ",")[[1]])
      print("doing an iteration on")
      print(paste("data_line = ", toString(data_line)))
      ## Insert some if statement logic here
      output = particle_filter_MVN_iter(data_line, priors)
      #priors = output
    }else{
      stop = TRUE
      close(f)
    }
    
    
   
    
  }
  return(output)
}

#MVN mixture Gibbs Sampler
data_size = 3000
K         = 5
d         = 2
shard_num = 1
scale     = 100
np        = 2000
data_file = paste('data/MVN_train_data_K=',
                  toString(K),'_data_size=',
                  toString(data_size),'.csv',
                  sep = "")

#priors
prob = matrix(rep(1/K, K), ncol = K, nrow = 1)
mean = matrix(rep(0, d*K), ncol = d, nrow = K)
sig  = t(matrix(as.vector(diag(d))*scale, ncol = K, nrow = d*d))
z    = rep(-1, np)
priors = list(prob = prob, mean = mean, sig = sig, z = z)
particle_filter_MVN(data_file, priors)
  
