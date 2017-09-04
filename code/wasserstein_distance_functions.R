organize_sample = function(worker_params_chain, iter_num){
  output = list()
  d = ncol(worker_params_chain$mean_chain[[1]])
  K = ncol(worker_params_chain$prob)
  Mu = matrix(NA, nrow = K, ncol = d)
  output$Pi = worker_params_chain$prob[iter_num,]
  for(k in 1:K){
    Mu[k, ] = worker_params_chain$mean_chain[[k]][iter_num,]
  }
  output$Mu = Mu
  return(output)
}

compute_K_matrix = function(Asamp, Bsamp){
  M_a = Asamp$mean #P_a['Mu']
  M_b = Bsamp$mean #P_b['Mu']
  ##print("M_a = " , M_a)
  ##print("len(M_a) = ", len(M_a))
  ##print("M_b = " , M_b)
  ##print("len(M_b) = ", len(M_b))
  row_num = nrow(M_a)
  col_num = nrow(M_b)
  K_matrix = matrix(0, nrow = row_num, ncol = col_num)
  for(r in 1:row_num){
    for(c in 1:col_num){
      #  #print("M_a[c] - M_b[r]")
      #  #print(M_a[c], " - ", M_b[r])
      K_matrix[r,c] = dist(M_a[c,] - M_b[r,])
    }
    
  }
  
  ##print("returning K_matrix")
  return(K_matrix)
}

make_A_eq = function(N1, N2){
  #output = np.matrix([[0]*N1*N2]*(N1+N2))
  output = matrix(0, ncol = N1*N2, nrow = N1+N2)
  #begin_col = range(0,N1*N2,N1)
  begin_col = seq(1, N1*N2, by = N1)
  #for row in range(N2):	
  #  #	#for col in range(0,N1*N2,N1):
  for(row in 1:N2){
    for( col in seq(1, N1*N2, by = N2)){
      #output[row, begin_col[row]:(begin_col[row]+N1)] = 1
      output[row, begin_col[row]:(begin_col[row]+N1 - 1)] = 1
    }
  }
  
  #next_cols = range(N2,N1+N2)
  next_rows = (N2+1):(N1+N2)
  ##print("next_cols = ", next_cols)
  #for ri in range(len(next_cols)):
  #  #begin_col = range(0,N1*N2,N1)
  #  for bc in begin_col:
  for( ri in 1:length(next_rows)){
    for(bc in begin_col){
      #      row = next_cols[ri]
      #      output[row, bc + ri] = 1    
      row = next_rows[ri]
      output[row, bc + ri - 1] = 1
    }
  }
  
  ##print("returning from make_A_eq")
  #return(output)
  return(output)
}

distance_between_samp_MVN_mix = function(Asamp, Bsamp){
  library(lpSolve)
  #gets the distance betwen
  K_matrix = mu_dist_mat(Asamp$mean, Bsamp$mean, Asamp$sig, Bsamp$sig) 
  #print(K_matrix)
  c = as.vector(t(K_matrix))
  #b = np.array(list(np.array([Q['Pi'], P['Pi']]).flat))
  
  #N1 = len(P['Mu'])
  #N2 = len(Q['Mu'])
  N1 = nrow(Asamp$mean)
  N2 = nrow(Bsamp$mean)
  A = make_A_eq(N1, N2)
  D = diag(ncol(A))
  constraints = rbind(A, D, D)
  
  condition = c(rep("=", nrow(A)), rep("<=", ncol(A)), rep(">=", ncol(A)))
  
  b = c(Asamp$prob[1,],Bsamp$prob[1,],  rep(1, ncol(A)), rep(0, ncol(A)))
  #print(A)
  #X_bounds = [(0, 1) for x in range(N1*N2)]
  #sf = 0.00001
  ##res = linprog(c, A_eq=A*sf, b_eq=b*sf, bounds=X_bounds, options={"disp": True, "tol":1e-6})
  #res = linprog(c, A_eq=A*sf, b_eq=b*sf, bounds=X_bounds, options={"disp": True, "bland": True , "tol": 1e-6})
  LPout = lp ("min", c, constraints, condition, b)
  #output$solution
  M = t(matrix(LPout$solution, nrow = N1))
  #print("M = ")
  #print(M)
  #print("rowSums(M) = ")
  #print(rowSums(M))
  #print("colSums(M) = ")
  #print(colSums(M))
  output = sum(K_matrix*M)
  return(output)
  ##print("return get_dist_between_particles")
  #print("c = ", c)
  #print("res[x] = ", res["x"])
  #print("res = ", res)
  #return(sum(c*res["x"]))
  #return(NULL)
}

#greedy_d = function(A, B, per_num){
#
#  min_val = Inf
#  min_wts = NA
#  #rep_num = 1
#  #temp = rep(NA, per_num)
#  for(i in 1:per_num){
#    Ar = reorder(A)
#    #print(Ar)
#    Br = reorder(B)
#    dist_mat   = mu_dist_mat(Ar$mean, Br$mean)
#    wts_mat = weights_matrix(Ar, Br, dist_mat)
#    temp = sum(wts_mat*dist_mat)
#    #print("greedy_d weights matrix")
#    #print(wts_mat)
#    #output[i] = temp
#    if(temp < min_val){
#      min_wts = wts_mat
#      min_val = temp
#    }
#  }
#  #print("greedy_d weights matrix")
#  #print(min_wts)
#  #print(output)
#  return(min_val)
#  
#}

#greedy_d_all_perm = function(A, B){
#  
#  min_val = Inf
#  min_wts = NA
#  N = length(A$prob)
#  perm_list = permn(1:N)
#  #print(perm_list)
#  per_num = length(perm_list)
#  #rep_num = 1
#  #temp = rep(NA, per_num)
#  for(i in 1:per_num){
#    index = perm_list[[i]]
#    #print(index)
#    Ar = A
#    Ar$prob = matrix(A$prob[,index], nrow = 1)
#    Ar$mean = A$mean[index,]
#    Ar$sig  = A$sig[index,]
#    #print(Ar)
#    #Br = reorder(B)
#    dist_mat   = mu_dist_mat(Ar$mean, B$mean)
#    wts_mat = weights_matrix(Ar, B, dist_mat)
#    temp = sum(wts_mat*dist_mat)
#    #print("greedy_d weights matrix")
#    #print(wts_mat)
#    #output[i] = temp
#    if(temp < min_val){
#      min_dist = dist_mat
#      min_wts  = wts_mat
#      min_val  = temp
#    }
#  }
#  #print("greedy_d weights matrix")
#  print(min_wts)
#  print(min_dist)
#  return(min_val)
#  
#}

#greedy_d_eq_wts = function(A, B, per_num){
#  
#  min_val = Inf
#  min_wts = NA
#  #rep_num = 1
#  #temp = rep(NA, per_num)
#  for(i in 1:per_num){
#    Ar = reorder(A)
#    #print("checking length(B$prob)")
#    #print(length(B$prob))
#    #matrix(A$prob[,index], nrow = 1)
#    Ar$prob = matrix(rep(1/length(A$prob), length(A$prob)), nrow = 1)
#    #print(Ar)
#    Br = reorder(B)
#    Br$prob = matrix(rep(1/length(B$prob), length(B$prob)), nrow = 1)
#    dist_mat   = mu_dist_mat(Ar$mean, Br$mean)
#    wts_mat = weights_matrix(Ar, Br, dist_mat)
#    temp = sum(wts_mat*dist_mat)
#    #print("greedy_d weights matrix")
#    #print(wts_mat)
#    #output[i] = temp
#    if(temp < min_val){
#      min_wts = wts_mat
#      min_val = temp
#    }
#  }
#  #print("greedy_d weights matrix")
#  #print(min_wts)
#  #print(output)
#  return(min_val)
#  
#}

#greedy_d_mass_order = function(A, B, per_num){
#  
#  min_val = Inf
#  min_wts = NA
#  Ar = reorder_decreasing(A)
#  Br = reorder_decreasing(B)
#  dist_mat   = mu_dist_mat(Ar$mean, Br$mean)
#  wts_mat = weights_matrix(Ar, Br, dist_mat)
#  output = sum(wts_mat*dist_mat)
#
#  return(output)
#  
#}

mu_dist_mat = function(mu1, mu2, all_C1, all_C2){
    #if("expm" %in% rownames(installed.packages()) == FALSE) {install.packages("expm")}
    #library(expm)
    N = nrow(mu1)
    d = ncol(mu1)
    output = matrix(0, nrow = N, ncol = N)
    for(r in 1:N){
        for(c in 1:N){
            #print(paste("postition: ", r, ",", c))
            #print(mu1[r,])
            #print(mu2[c,])
            #C1       = matrix(all_C1[r,], nrow = d)
            #print(C1)
            #C2       = matrix(all_C2[c,], nrow = d)
            #print(C2)
            #sqrt_C2  = sqrtm(C2)
            #print(sqrt_C2)
            #temp_mat = sqrtm(sqrt_C2 %*% C1 %*% sqrt_C2) 
            #print(temp_mat)
            output[r,c] = sqrt(sum((mu1[r,] - mu2[c,])^2)) #+ sum(diag(C1 + C2 - 2*temp_mat))#note: sum(diag(M)) == trace(M)
            print("last dist ran")
        }
    }
  #print(output)
  return(output)
}

weights_matrix = function(A, B, dist_mat){
  N          = ncol(A$prob)
  weight_mat =  matrix(0, nrow = N, ncol = N)
  #try many permutations and starting points
  wtsA_before = wtsA = A$prob[1,]
  wtsb_before = wtsB = B$prob[1,]
  valid_row = rep(TRUE, N)
  valid_col = rep(TRUE, N)
  for( r in 1:nrow(dist_mat)){
    
    while(valid_row[r]){
      
      c_index = which.min(dist_mat[r,]) 
      if(wtsA[r] <= wtsB[c_index]){
        
        weight_mat[r,c_index] = wtsA[r]
        wtsB[c_index] = wtsB[c_index] - wtsA[r]
        wtsA[r] = 0
        valid_row[r] = FALSE
        
      }else if(wtsA[r] > wtsB[c_index]){
        
        weight_mat[r,c_index] = wtsB[c_index]
        wtsA[r] = wtsA[r] - wtsB[c_index]
        wtsB[c_index] = 0
        valid_col[c_index] = FALSE
        
      }
      dist_mat[r,c_index] = Inf
      if(all.equal(rowSums(weight_mat), A$prob[1,])==TRUE & 
         all.equal(colSums(weight_mat), B$prob[1,])==TRUE){
        break
      }
    }
  }
  #print(weight_mat)
  return(weight_mat)
}

reorder = function(A){
  N      = length(A$prob)
  index  = sample(N,N)
  A$prob = matrix(A$prob[,index], nrow = 1)
  A$mean = A$mean[index,]
  A$sig  = A$sig[index,]
  return(A)
}

reorder_decreasing = function(A){
  N      = length(A$prob)
  index  = order(A$prob, decreasing = TRUE)
  A$prob = matrix(A$prob[,index], nrow = 1)
  A$mean = A$mean[index,]
  A$sig  = A$sig[index,]
  return(A)
}

R_barpost = function(colMeasure, distMat){
  library(Matrix )
  library(lpSolve)
  #Code from original function is below and was originally writen in Matlab
  #by Srivastava for Wasserstein Barrycenter paper
  
  #nsubset = length(distMat);
  nsubset = length(colMeasure)
  #print(paste("nsubset = ", nsubset))
  #vecColMeasure = cell2mat(cellfun(@(x) x',  colMeasure, 'UniformOutput', false))';
  vecColMeasure = matrix(unlist(colMeasure), ncol = 1)
  #print("vecColMeasure = ")
  #print(vecColMeasure)
  #nsample = size(vecColMeasure, 1);
  nsample = length(vecColMeasure)
  #print("nsample = ")
  #print(nsample)
  
  #vecDistMat = cell2mat(distMat);
  #vecDistMat = vecDistMat(:);
  vecDistMat = matrix(distMat, ncol = 1)
  #print("vecDistMat = ")
  #print(vecDistMat)

  meas_sizes = unlist(lapply(colMeasure, length))
  nr_fmat_cell = sum(meas_sizes)
  fmat_cell_temp = list()
  for(i in 1:nsubset){
    fmat_cell_temp[[i]] = matrix(unlist(rep(list(diag(nr_fmat_cell)),meas_sizes[i])), nrow = nr_fmat_cell)
    #print("fmat_cell_temp[[i]] = ")
    #print(fmat_cell_temp[[i]])
  }
  fmat = as.matrix(.bdiag(fmat_cell_temp))
  #fmatCell = cellfun(@(x) kron(ones(1, size(x, 1)), eye(nsample)), colMeasure, 'UniformOutput', false);
  #print("fmat = ")
  #print(fmat)
  #hmatCell = cellfun(@(x) kron(eye(size(x, 1)), ones(1, nsample)), colMeasure, 'UniformOutput', false);
  hmat_cell_temp = list()
  for(i in 1:nsubset){
    ones    = matrix( rep(1, nsample), nrow = 1)
    to_diag = rep(list(ones), meas_sizes[i])
    hmat_cell_temp[[i]] = as.matrix(.bdiag(to_diag))  
    #print("hmat_cell_temp[[i]] = ")
    #print(hmat_cell_temp[[i]])
  }
  hmat = as.matrix(.bdiag(hmat_cell_temp))
  #print("hmat = ")
  #print(hmat)
  #gmat = kron(ones(nsubset, 1), eye(nsample)); % G
  pre_gmat = rep(list(diag(nsample)),nsubset)
  gmat = do.call(rbind, pre_gmat) 
  #print("gmat = ")
  #print(gmat)
  
  #% Aeq in the matlab linprog fun.
  #aMat = [zeros(1, nsample^2) ones(1, nsample);
  #		fmat                -gmat;
  #		hmat                zeros(nsample, nsample)
  #	];
  zeros_mat = matrix(0, nrow = nsample, ncol = nsample)
  top_row    = matrix(c(rep(0, nsample*nsample), rep(1, nsample)), nrow = 1)
  #print(dim(top_row))
  second_row = cbind(fmat, -gmat)
  #print(dim(second_row))
  third_row  = cbind(hmat, zeros_mat)
  #print(dim(third_row))
  aMat = rbind(top_row, second_row, third_row)
  #print("aMat = ")
  #print(aMat)
  
  #bVec = [1;
  #		zeros(nsubset * nsample, 1);
  #		vecColMeasure
  #	];
  bVec = matrix(c(1, 
                  rep(0, nsubset * nsample), 
                  vecColMeasure, 
                  rep(0, nsample^2 + nsample),
                  rep(1, nsample^2 + nsample)), ncol = 1)
  #print("bVec = ")
  #print(bVec)

  #costVec = [vecDistMat; 
  #zeros(nsample, 1)];
  costVec = matrix(c(vecDistMat, rep(0, nsample)), ncol = 1)
  #print("costVec = ")
  #print(costVec)

  #% upper bds = 1 and lower bds = 0
  #lbd = zeros(nsample^2 + nsample, 1);
  #ubd = ones(nsample^2 + nsample, 1);
  diags = diag(nsample^2 + nsample)
  
  condition = c(rep("=", nrow(aMat)), rep(">=", ncol(diags)), rep("<=", ncol(diags)))
  constraints = rbind(aMat, diags, diags)
  #[optSol, obj, exitflag, output, lambda] = linprog(costVec, [], [], aMat, bVec, lbd, ubd);
  LPout = lp("min", costVec, constraints, condition, bVec)
  sol = LPout$solution
  output = sol[(length(sol) - nsample + 1):length(sol)]
  return(output)
}
get_colMeasure = function(all_chains){
  acL = lapply(all_chains, length)
  output = list()
  for(i in 1:length(acL)){
    N_i = acL[[i]]
    output[[i]] = matrix(1/N_i, nrow = N_i, ncol = 1)  
  }
  return(output)
}

get_distMat = function(all_chains){
  
  ULC = unlist_chians(all_chains)
  acL = length(ULC)

  #cores=detectCores()
  cl <- makeCluster(acL)#cores[1]) 
  registerDoParallel(cl)
  
  pre_distMat = foreach(cn_1 = 1:acL, .packages = "lpSolve", 
                        .export=c("distance_between_samp_MVN_mix", 
                                  "mu_dist_mat", "make_A_eq")) %dopar%{
    
    d_outs = rep(0, acL)
    
    for(cn_2 in cn_1:acL){
      #cn_2_index = which(ind_set == j)
      if(cn_1 < cn_2){
        d_outs[cn_2] = distance_between_samp_MVN_mix(ULC[[cn_1]], ULC[[cn_2]])
      }
    }
    d_outs
  }
  stopCluster(cl)
  unPart_distMat = matrix(unlist(pre_distMat), 
                          nrow = length(pre_distMat), 
                          byrow = TRUE)
  
  output = unPart_distMat + t(unPart_distMat)
  return(output)
}

unlist_chians = function(all_chains){
  NC = length(all_chains)
  output = list()
  count = 1
  for(nc in 1:NC){
    for(index in 1:length(all_chains[[nc]])){
      output[count] = all_chains[[nc]][index]
      count = count + 1
    }
  }
  return(output)
}






















