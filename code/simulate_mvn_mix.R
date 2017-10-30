library(MASS) #to sample MVN
library(MCMCpack) # to sample Inv. Wishart

#setwd("/Users/o/Google Drive/school/Williamson Research/communal_monte_carlo_R")
sample_MVN_MIX = function(n, K, class_counts, MU, COV){
  for(i in 1:K){
    if(class_counts[i] == 0){next}
    val <- mvrnorm(class_counts[i],mu=MU[i,],Sigma=COV[,,i])
    if(i == 1){allval = val}
    else{allval = rbind(allval,val)}
  }
  output <- allval[sample(n,n),]
  return(output)
}

scale        = 100
file_num     = 11
d            = 2
K            = 20
n            = 1e+06
prob         = rdirichlet(1,rep(1,K))
class_counts = rmultinom(1, n, prob)
COV          = replicate(K, solve(riwish(d, diag(d))))
MU           = mvrnorm( K, mu=rep(0, d),Sigma=diag(d)*sqrt(scale))

true_params = list(d = d, 
                   K = K, 
                   n = n, 
                   prob = prob, 
                   class_counts = class_counts, 
                   MU = MU, 
                   COV = COV)

save(true_params, file= paste0("data/K=",K,"_n=",n,"/true_param.RData"))  

for(fn in 1:file_num){
  file_out = sample_MVN_MIX(n, K, class_counts, MU, COV)
  write.csv(file_out, file = paste0("data/K=",K,"_n=",n,"/d=",d,"_n=",n,"_file_num_",fn,".csv"),
            row.names = FALSE, col.names = NA)  
}

