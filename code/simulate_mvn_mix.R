library(MASS) #to sample MVN
library(MCMCpack) # to sample Inv. Wishart

setwd("/Users/o/Google Drive/school/Williamson Research/communal_monte_carlo_R")
sample_MVN_MIX = function(K, class_counts, MU, COV){
  for(i in 1:K){
    val <- mvrnorm(class_counts[i],mu=MU[i,],Sigma=COV[,,i])
    if(i == 1){allval = val}
    else{allval = rbind(allval,val)}
  }
  output <- allval[sample(n,n),]
  return(output)
}

file_num     = 10
d            = 2
K            = 10
n            = 1000
prob         = rdirichlet(1,rep(1,K))
class_counts = rmultinom(1, n, prob)
COV          = replicate(K, solve(riwish(d, diag(d))))
MU           = mvrnorm( K, mu=rep(0, d),Sigma=diag(d)*sqrt(10))

true_params = list(d = d, 
                   K = K, 
                   n = n, 
                   prob = prob, 
                   class_counts = class_counts, 
                   MU = MU, 
                   COV = COV)
for(fn in 1:file_num){
  file_out = sample_MVN_MIX(K, class_counts, MU, COV)
  write.csv(file_out, file = paste0("data/K=",K,"/d=",d,"_n=",n,"_file_num_",fn,".csv"),
            row.names = FALSE, col.names = NA)  
}
save(true_params, file= paste0("data/K=",K,"/true_param.RData"))  

