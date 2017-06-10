<<<<<<< HEAD

setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")


particle_filter_MVN_iter = function(data_line, priors){
  return()
}


particle_filter_MVN = function(file_name, priors){
  
  stop = FALSE
  f = file(file_name, "r")
  while(!stop) {
    
    next_line = readLines(f, n = 1)
    if(length(next_line) != 0){
      data_line = as.numeric(strsplit(next_line, ",")[[1]])
      print(paste("data_line = ", data_line))
    }else{
      stop = TRUE
      close(f)
    }
    
    
    ## Insert some if statement logic here
    output = particle_filter_MVN_iter(data_line, priors)
    #priors = output
    
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
  
=======

setwd("/Users/o/Google Drive/school/Williamson Research/Implementations/python code/test_runs")

#MVN mixture Gibbs Sampler
data_size = "100000"
K         = "30"
shard_num = 1
data_file = paste('data/MVN_train_data_K=',toString(K),'_data_size=',toString(data_size),'.csv', sep = "")
dat = read.csv(data_file, header = FALSE,  nrows = 1, skip = 0)
dat


stop = FALSE
f = file(data_file, "r")
while(!stop) {
  next_line = readLines(f, n = 1)
  data_line = as.numeric(strsplit(next_line, ",")[[1]])
  ## Insert some if statement logic here
  if(length(next_line) == 0) {
    stop = TRUE
    close(f)
  }
}
>>>>>>> 11e6a07dec4fe6e96c02590920339c641f03c523
