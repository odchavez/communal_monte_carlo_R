#Original Barycenter Code from Barycenter package
#function (images, maxIter = 10, lambda = 60/median(costMatrix)) 
#{
#    time <- proc.time()
#    dimension <- dim(images[[1]])
#    n <- dimension[1] * dimension[2]
#    coord1 <- seq(0, 1, length.out = dimension[2])
#    coord2 <- rev(seq(0, 1, length.out = dimension[1]))
#    coordinates <- expand.grid(coord1, coord2)
#    costMatrix <- as.matrix(dist(coordinates, diag = TRUE, upper = TRUE))
#    a_tild <- rep(1/n, n)
#    a_hat <- rep(1/n, n)
#    t_0 <- 2
#    t <- t_0
#    for (i in 1:maxIter) {
#        beta <- (t + 1)/2
#        a <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
#        ALPHA <- 0
#        for (j in 1:length(images)) {
#            ALPHA <- Subgradient(a, t(images[[j]]), costMatrix, 
#                lambda) + ALPHA
#        }
#        ALPHA <- (1/length(images)) * ALPHA
#        a_tild <- a_tild * exp(-(t_0) * beta * ALPHA)
#        a_tild <- a_tild/sum(a_tild)
#        a_hat <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
#        t <- t + 1
#    }
#    a <- matrix(a, dimension[1], dimension[2], byrow = TRUE)
#    a.temp <- a[, nrow(a):1]
#    print(image(a.temp))
#    print(proc.time() - time)
#    return(a)
#}
#<environment: namespace:Barycenter>
library(Barycenter)
Barycenter_measure = function(colMeasure, distMat, maxIter = 100, lambda){# = 1/median(distMat)){
  
    time <- proc.time()
    #dimension <- dim(images[[1]])
    #n <- dimension[1] * dimension[2]
    n = length(unlist(colMeasure))
    #coord1 <- seq(0, 1, length.out = dimension[2])
    #coord2 <- rev(seq(0, 1, length.out = dimension[1]))
    #coordinates <- expand.grid(coord1, coord2)
    #costMatrix <- as.matrix(dist(coordinates, diag = TRUE, upper = TRUE))
    
    a_tild <- rep(1/n, n)
    a_hat <- rep(1/n, n)
    t_0 <- 2
    t <- t_0
    for (i in 1:maxIter) {
        beta <- (t + 1)/2
        a <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
        ALPHA <- 0
        for (j in 1:length(colMeasure)) {
            n_j = length(colMeasure[[j]])
            ALPHA <- Subgradient(a, colMeasure[[j]], distMat[,(n_j*(j-1)+1):(n_j*j)], 
                lambda) + ALPHA
        }
        ALPHA <- (1/length(colMeasure)) * ALPHA
        a_tild <- a_tild * exp(-(t_0) * beta * ALPHA)
        a_tild <- a_tild/sum(a_tild)
        a_hat <- (1 - 1/beta) * a_hat + (1/beta) * a_tild
        t <- t + 1
    }
    #a <- matrix(a, dimension[1], dimension[2], byrow = TRUE)
    #a.temp <- a[, nrow(a):1]
    #print(image(a.temp))
    print(proc.time() - time)
    return(a)
}
