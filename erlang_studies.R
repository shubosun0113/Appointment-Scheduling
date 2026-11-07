source("Main.R")
source("helper-functions.R")
source("kdpp-sample.R")
library(extraDistr)

# A t interval for the mean of the independent outer simulation results.
# For erlang_real, each outer result already averages 100 service-time draws.
mean_ci95 <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x))) {
    stop("A confidence interval needs at least two finite simulation results")
  }
  centre <- mean(x)
  half_width <- qt(0.975, df = length(x) - 1L) * sd(x) / sqrt(length(x))
  data.frame(n = length(x), mean = centre,
             lower = centre - half_width, upper = centre + half_width)
}

final_cost_ci95 <- function(obj) {
  lapply(obj[c("N", "OT", "WT", "IT", "V")], mean_ci95)
}

local_search_plus_dpp_erlang <- function(T,slot_length,mu,m,q0,q_breaks=NULL,ci,cw,co,cw_walkin=0.075,
                                         omega=0,dist,nonpunc_prob = c(1/3,1/3,1/3),K=3,start_type="random",B0=100,
                                         Round=FALSE,model='erlang',cl=NULL,SetSeed=FALSE){
  result5 <- list()
  for (start in 1:5) {
    result <- list()
    #Read parameters from files
    # PAR <- SET[[case]][[1]]
    ###Set up beta-binomial distribution
    
    w <- c(
      ci, cw, co, cw_walkin
    )
    
    cov <- 1/sqrt(K) # Erlang coefficient of variation
    eta <- K/mu
    
    aa <- ParBetaBin(mu,cov,m)
    alpha <- aa$alpha
    beta <- aa$beta
    
    
    
    smp <- 0:m
    prob <- dbbinom(smp, m, alpha , beta , log = FALSE)
    
    #
    grid <- seq(0, T-slot_length, by=slot_length)
    sch.set <- grid
    ngrid <- length(grid)
    
    if(length(q0)==1){
      q.set <- rep(q0, ngrid)
    } else if(length(q0)==ngrid && is.null(q_breaks)){
      q.set <- q0
    } else { # designed for 2 or 3 time varying time probability
      if (is.null(q_breaks)) {
        stop("q_breaks must be provided for period-specific probabilities")
      }
      
      if (length(q_breaks) != length(q0) - 1L) {
        stop("length(q_breaks) must equal length(q0) - 1")
      }
      
      if (any(q_breaks <= 0) ||
          any(q_breaks >= T) ||
          any(diff(q_breaks) <= 0)) {
        stop("q_breaks must be increasing and strictly between 0 and T")
      }
      
      q.ind <- findInterval(grid, q_breaks)+1
      q.set <- q0[q.ind]
    } 
    
    q.average <- mean(q.set)
    N0 <- max(1, round(T / ((1 - q.average) * mu)))
    if(start_type=="even"){
      ## Evenly distribute the patients over the session
      ideal_times <- (0:(N0 - 1)) * T / N0
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="random"){
      set.seed(start*1000)
      sch.opt0 <- sample(sch.set,N0)
    } else if(start_type=="front"){
      ## Evenly distribute the patients over the 1st-half session
      ideal_times <- sch.set[1:round(length(sch.set)/2)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="back"){
      ## Evenly distribute the patients over the 2nd-half session
      ideal_times <- sch.set[round(length(sch.set)/2+1):length(sch.set)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="dpp"){
      if(start!=0){
        #The initial L matrix
        L <- diag(ngrid)
        
        ##Find the initial schedule
        Lnew <- L
        L.decomposed <- decompose.matrix(Lnew)
        kdpp.sample <- k.dpp.sample(L.decomposed, k=floor(min(N0,length(sch.set))-2))
        sch.opt0 <- sch.set[kdpp.sample]
      }
    }
    
    sch.opt0 <- sort(sch.opt0)
    
    ## Convert the schedule to the number of appointments in each grid slot
    slot_index <- match(sch.opt0, grid)
    ## Patient-specific no-show probabilities
    q.opt0 <- q.set[slot_index]
    
    idx_sch_ini <- idx_sch_opt0 <- tabulate(
      slot_index,
      nbins = ngrid
    )
    
    ## No-show probability associated with each possible grid time
    # q.set <- rep(q0, ngrid)
    
    ###Set up the no show probability and nonpunctuality
    # omega=0
    # dist="three_point" 
    
    ###First compute the cost for the given 
    obj.opt0 <- obj.ini <-  PFmeasure_parallel(sch.opt0,dist=dist,q=q.opt0,omega,T=T, Round=Round,model=model,nsim=10000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
    
    #The initial L matrix
    
    sch.ini <- rep(sch.set,idx_sch_ini)
    q.ini <- rep(q.set,idx_sch_ini)
    
    loss.min <- mean(obj.ini$V)
    print(mean(loss.min))
    
    idx_sch_new <- idx_sch_ini
    
    ##Update the L matrix
    sch.new <- rep(sch.set,idx_sch_new) 
    q.new <- rep(q.set,idx_sch_new)  
    
    obj.path <-  WorkloadPath_parallel(sch.new,dist=dist,q=q.new,sch.set=sch.set,omega,T=T, Round=Round,model=model,nsim=10000,p=prob,y=smp,cl=cl,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
    wpath <- obj.path$workload_path
    aa <- t(wpath)%*%wpath
    obj.aa <- eigen(aa)
    P <- obj.aa$vectors
    d <- abs(obj.aa$values)
    d[d<=0] <- d[d<=0]+0.1*min(d[d>0])
    aa <- P%*%diag(d)%*%t(P)
    bb <- diag(sqrt(1/diag(aa)))
    L <- bb%*%aa%*%bb
    
    Iter <- 2
    iter <- 1
    Loss <- NULL
    NP <- NULL
    
    t.begin <- Sys.time()
    
    while(iter<=Iter){
      
      
      if(iter==1) { 
        try <- 5
        nsim.inner <- 100
        fold.max <- 1
      }
      
      if(iter==2){
        try=5
        nsim.inner <- 100
        fold.max <- 5
      } 
      
      
      fold <- 1
      
      while(fold <=fold.max){
        
        if(iter>1){
          change1 <- 1
          
          while (change1!=0) {
            change1 <- 0
            ###Perform random walk search
            obj.tmp <-  LocalPFmeasure_parallel(sch.new,sch.set,dist=dist,q.set=q.set,omega,T=T, Round=Round,model=model,nsim=1000,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
            aa <- sort(obj.tmp$V,index=TRUE)$ix
            
            ####Finish search until fails "try" times
            l <- 1
            count <- 1
            while (l <= try && count <= length(aa)) {
              kdd_tmp <- obj.tmp$KddSample[aa[count],]
              sch_tmp <- rep(sch.set,kdd_tmp) 
              q_tmp <- rep(q.set,kdd_tmp) 
              
              obj.tmp2 <-  PFmeasure_parallel(sch_tmp,dist=dist,q=q_tmp,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,w=w,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
              loss.tmp <- mean(obj.tmp2$V)
              
              
              if(loss.tmp<loss.min) {
                sch.new <- sch_tmp
                loss.min <- loss.tmp
                l <- 1
                change1 <- change1 + 1
              }else{
                
                l <- l + 1
              }
              count <- count + 1
              print(c(iter,length(sch.new),loss.tmp,loss.min))
            }
            Loss <- c(Loss,loss.min)
            NP <- c(NP,length(sch.new))
          }
        }
        
        
        
        if(iter==1) k <- 1
        
        if(iter>=2) k <- fold
        
        while (k <=fold) {
          L.d0 <- numeric(ngrid)
          idx <- sch.set%in%sch.new
          L.d0[idx] <- 1
          change <- 0
          
          if(sum(!idx)!=0){
            if(iter==1)    L.d0[!idx] <- 0.1
            
            if(iter>=2)    L.d0[!idx] <- 0.1/2^(fold-1)
            
            Lnew <- diag(sqrt(L.d0))%*%L%*%diag(sqrt(L.d0))
            if(iter==1) {
              obj.tmp <-  OptPFmeasure_parallel(Lnew,sch.new,sch.set,dist=dist,q.set=q.set,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,fix.nslots=FALSE,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K) 
            }else
              if(iter==2&fold==1)      obj.tmp <-  OptPFmeasure_parallel(Lnew,sch.new,sch.set,dist=dist,q.set=q.set,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,fix.nslots=FALSE,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K) else
                obj.tmp <-  OptPFmeasure_parallel(Lnew,sch.new,sch.set,dist=dist,q.set=q.set,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,fix.nslots=TRUE,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
            aa <- sort(obj.tmp$V,index=TRUE)$ix
            
            l <- 1
            count <- 1
            ####Finish search until fails "try" times
            while (l <= try&count<=length(aa)) {
              sch_tmp <- rep(sch.set,obj.tmp$KddSample[aa[count],])
              q_tmp <- rep(q.set,obj.tmp$KddSample[aa[count],])
              obj.tmp2 <-  PFmeasure_parallel(sch_tmp,dist=dist,q=q_tmp,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,w=w,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
              loss.tmp <- mean(obj.tmp2$V)
              
              
              if(loss.tmp<loss.min) {
                
                sch.new <- sch_tmp
                loss.min <- loss.tmp
                l <- 1
                change <- change + 1
              }else{
                l <- l + 1
              }
              count <- count + 1
              print(c(iter,length(sch.new),loss.tmp,loss.min,k,fold))
              
            }
            
            Loss <- c(Loss,loss.min)
            NP <- c(NP,length(sch.new))
          }
          if(change==0) k <- k +1 else { 
            
            ##Then update the L matrix
            idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new))
            q.new <- rep(q.set,idx_sch_new)  
            
            obj.path <-  WorkloadPath_parallel(sch.new,dist=dist,q=q.new,sch.set=sch.set,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
            wpath <- obj.path$workload_path
            aa <- t(wpath)%*%wpath
            obj.aa <- eigen(aa)
            P <- obj.aa$vectors
            d <- abs(obj.aa$values)
            d[d<=0] <- d[d<=0]+0.1*min(d[d>0])
            aa <- P%*%diag(d)%*%t(P)
            bb <- diag(sqrt(1/diag(aa)))
            L <- bb%*%aa%*%bb
            k <- fold + 2
          } 
          
        }
        if(k==(fold+1)) fold <- fold + 1 
        
        print(rbind(start,Loss,NP)) 
        
      }
      iter <- iter +1
    }
    
    time <- as.numeric(difftime(Sys.time(), t.begin, units = "secs"))
    
    
    idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new))
    q.new <- rep(q.set,idx_sch_new)  
    obj.new <-  PFmeasure_parallel(sch.new,dist=dist,q=q.new,omega,T=T, Round=Round,model=model,nsim=100000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
    print(mean(obj.new$V))
    print(round(sch.new))
    # stopCluster(cl)
    
    aN0 <- mean(obj.opt0$N)
    aOT0 <- mean(obj.opt0$OT)
    aWT0 <- mean(obj.opt0$WT)/(length(sch.opt0)*mean(1-q.opt0))
    aIT0 <- mean(obj.opt0$IT)
    aV0 <- mean(obj.opt0$V)
    
    aN <- mean(obj.new$N)
    aOT <- mean(obj.new$OT)
    aWT <- mean(obj.new$WT)/(length(sch.new)*mean(1-q.new))
    aIT <- mean(obj.new$IT)
    aV <- mean(obj.new$V)
    
    result$par <- list(
      T = T,
      slot_length = slot_length,
      mu = mu,
      cv = cov,
      K = K,
      eta = eta,
      service_support = c(0, Inf),
      auxiliary_beta_binomial_m = m,
      no_show_probability = q0,
      weights = w,
      nonpunc_dist = dist,
      omega = omega,
      model = model,
      B0 = B0,
      nonpunc_prob = nonpunc_prob
    )
    
    result$d <- c(result$d,d)
    result$time <- c(result$time,time)
    
    result$sch.set <- c(result$sch.set,list(sch.set)) 
    # result$par <- PAR
    result$sch.opt0 <- sch.opt0
    result$perf0 <- rbind(result$perf0,c(aN0,aOT0,aWT0,aIT0,aV0))
    
    result$perfnew <- rbind(result$perfnew,c(aN,aOT,aWT,aIT,aV))
    
    result$sch.new <- c(result$sch.new,list(sch.new))
    result$idx_sch_new <- c(result$idx_sch_new,list(idx_sch_new))
    
    result$Loss <- c(result$Loss,list(Loss))  ##Loss trajectory of iterations
    result$NP <- c(result$NP,list(NP)) ##Number of scheduled patients of iterations
    
    result$ci95 <- final_cost_ci95(obj.new)
    result5[[start]] <- result
  }
  result5
}


local_search_only_erlang <- function(T,slot_length,mu,m,q0,q_breaks=NULL,ci,cw,co,cw_walkin=0.075,
                                     omega=0,dist,SetSeed=FALSE,nonpunc_prob = c(1/3,1/3,1/3),K=3,start_type="random",B0,
                                     Round=FALSE,model="erlang",cl=NULL){
  result5 <- list()
  for (start in 1:5) {
    set.seed(start*1000) 
    result <- list()
    #Read parameters from files
    # PAR <- SET[[case]][[1]]
    ###Set up beta-binomial distribution
    
    w <- c(
      ci, cw, co, cw_walkin
    )
    
    cov <- 1/sqrt(K) # Erlang coefficient of variation
    eta <- K/mu
    
    aa <- ParBetaBin(mu,cov,m)
    alpha <- aa$alpha
    beta <- aa$beta
    ####
    #bb <- rbbinom(n=10000, size=m, alpha, beta)
    #mean(bb)
    #sd(bb)/mean(bb)
    
    smp <- 0:m
    prob <- dbbinom(smp, m, alpha , beta , log = FALSE)
    
    #
    grid <- seq(0, T-slot_length, by=slot_length)
    sch.set <- grid
    ngrid <- length(grid)
    
    if(length(q0)==1){
      q.set <- rep(q0, ngrid)
    } else if(length(q0)==ngrid && is.null(q_breaks)){
      q.set <- q0
    } else { # designed for 2 or 3 time varying time probability
      if (is.null(q_breaks)) {
        stop("q_breaks must be provided for period-specific probabilities")
      }
      
      if (length(q_breaks) != length(q0) - 1L) {
        stop("length(q_breaks) must equal length(q0) - 1")
      }
      
      if (any(q_breaks <= 0) ||
          any(q_breaks >= T) ||
          any(diff(q_breaks) <= 0)) {
        stop("q_breaks must be increasing and strictly between 0 and T")
      }
      
      q.ind <- findInterval(grid, q_breaks)+1
      q.set <- q0[q.ind]
    } 
    
    q.average <- mean(q.set)
    N0 <- max(1, round(T / ((1 - q.average) * mu)))
    if(start_type=="even"){
      ## Evenly distribute the patients over the session
      ideal_times <- (0:(N0 - 1)) * T / N0
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="random"){
      set.seed(start*1000)
      sch.opt0 <- sample(sch.set,N0)
    } else if(start_type=="front"){
      ## Evenly distribute the patients over the 1st-half session
      ideal_times <- sch.set[1:round(length(sch.set)/2)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="back"){
      ## Evenly distribute the patients over the 2nd-half session
      ideal_times <- sch.set[round(length(sch.set)/2+1):length(sch.set)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="dpp"){
      if(start!=0){
        #The initial L matrix
        L <- diag(ngrid)
        
        ##Find the initial schedule
        Lnew <- L
        L.decomposed <- decompose.matrix(Lnew)
        kdpp.sample <- k.dpp.sample(L.decomposed, k=floor(min(N0,length(sch.set))-2))
        sch.opt0 <- sch.set[kdpp.sample]
      }
    }
    
    sch.opt0 <- sort(sch.opt0)
    
    ## Convert the schedule to the number of appointments in each grid slot
    slot_index <- match(sch.opt0, grid)
    ## Patient-specific no-show probabilities
    q.opt0 <- q.set[slot_index]
    
    idx_sch_ini <- idx_sch_opt0 <- tabulate(
      slot_index,
      nbins = ngrid
    )
    
    ###Set up the no show probability and nonpunctuality
    # dist="three_point" 
    
    # cl <- makeCluster(50)
    # registerDoSNOW(cl)
    
    ###First compute the cost for the given 
    obj.opt0 <- obj.ini <-  PFmeasure_parallel(sch.opt0,dist=dist,q=q.opt0,omega,T=T, Round=Round,model=model,nsim=10000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
    loss.min <- mean(obj.ini$V)
    print(loss.min)
    
    ##Update the L matrix
    idx_sch_new <- idx_sch_ini
    sch.new <- sch.ini <- rep(sch.set,idx_sch_new) 
    q.new <- q.ini <- rep(q.set,idx_sch_new)  
    
    # obj.path <-  WorkloadPath_parallel(sch.new,dist=dist,q=q.new,sch.set=sch.set,omega,T=T, Round=TRUE,nsim=10000,model="discrete",p=prob,y=smp,cl=cl,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
    # wpath <- obj.path$workload_path
    # aa <- t(wpath)%*%wpath
    # obj.aa <- eigen(aa)
    # P <- obj.aa$vectors
    # d <- abs(obj.aa$values)
    # d[d<=0] <- d[d<=0]+0.1*min(d[d>0])
    # aa <- P%*%diag(d)%*%t(P)
    # bb <- diag(sqrt(1/diag(aa)))
    # L <- bb%*%aa%*%bb
    
    # Iter <- 2
    # iter <- 1
    
    # while(iter<=Iter){
    
    
    # if(iter==1) { 
    #   try <- 5
    #   nsim.inner <- 100
    #   fold.max <- 1
    # }
    # 
    # if(iter==2){
    #   try=5
    #   nsim.inner <- 100
    #   fold.max <- 5
    # } 
    
    
    # fold <- 1
    
    # while(fold <=fold.max){
    # while(fold <=fold.max){  
    #   if(iter>1){
    Loss <- NULL
    NP <- NULL
    nsim.inner <- 100
    change1 <- 1
    try <- 5
    t.begin <- Sys.time()
    
    while (change1!=0) {
      change1 <- 0
      ###Perform random walk search
      obj.tmp <-  LocalPFmeasure_parallel(sch.new,sch.set,dist=dist,q.set=q.set,omega,T=T, Round=Round,model=model,nsim=1000,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
      aa <- sort(obj.tmp$V,index=TRUE)$ix
      
      ####Finish search until fails "try" times
      l <- 1
      count <- 1
      while (l <= try && count <= length(aa)) {
        kdd_tmp <- obj.tmp$KddSample[aa[count],]
        sch_tmp <- rep(sch.set,kdd_tmp) 
        q_tmp <- rep(q.set,kdd_tmp) 
        
        obj.tmp2 <-  PFmeasure_parallel(sch_tmp,dist=dist,q=q_tmp,omega,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,w=w,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
        loss.tmp <- mean(obj.tmp2$V)
        
        
        if(loss.tmp<loss.min) {
          sch.new <- sch_tmp
          q.new <- q_tmp
          loss.min <- loss.tmp
          l <- 1
          change1 <- change1 + 1
        }else{
          
          l <- l + 1
        }
        count <- count + 1
        print(c(start,length(sch.new),loss.tmp,loss.min))
      }
      Loss <- c(Loss,loss.min)
      NP <- c(NP,length(sch.new))
    }
    #   }
    # 
    # iter <- iter +1
    # }
    
    time <- as.numeric(difftime(Sys.time(), t.begin, units = "secs"))
    
    # idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new))
    slot_id_new <- match(sch.new, sch.set)
    if (anyNA(slot_id_new)) {
      stop("The final schedule contains times outside sch.set")
    }
    idx_sch_new <- tabulate(slot_id_new,nbins = ngrid)
    q.new <- rep(q.set,idx_sch_new)  
    obj.new <-  PFmeasure_parallel(sch.new,dist=dist,q=q.new,omega,T=T, Round=Round,model=model,nsim=100000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,SetSeed=TRUE,eta=eta,K=K)
    print(mean(obj.new$V))
    print(round(sch.new))
    # stopCluster(cl)
    
    aN0 <- mean(obj.opt0$N)
    aOT0 <- mean(obj.opt0$OT)
    aWT0 <- mean(obj.opt0$WT)/(length(sch.opt0)*mean(1-q.opt0))
    aIT0 <- mean(obj.opt0$IT)
    aV0 <- mean(obj.opt0$V)
    
    aN <- mean(obj.new$N)
    aOT <- mean(obj.new$OT)
    aWT <- mean(obj.new$WT)/(length(sch.new)*mean(1-q.new))
    aIT <- mean(obj.new$IT)
    aV <- mean(obj.new$V)
    
    result$par <- list(
      T = T,
      slot_length = slot_length,
      mu = mu,
      cv = cov,
      K = K,
      eta = eta,
      service_support = c(0, Inf),
      auxiliary_beta_binomial_m = m,
      no_show_probability = q0,
      weights = w,
      nonpunc_dist = dist,
      omega = omega,
      model = model,
      B0 = B0,
      nonpunc_prob = nonpunc_prob
    )
    
    result$slot_length <- slot_length
    result$time <- c(result$time,time)
    
    result$sch.set <- c(result$sch.set,list(sch.set)) 
    # result$par <- PAR
    result$sch.opt0 <- sch.opt0
    result$perf0 <- rbind(result$perf0,c(aN0,aOT0,aWT0,aIT0,aV0))
    
    result$perfnew <- rbind(result$perfnew,c(aN,aOT,aWT,aIT,aV))
    
    result$sch.new <- c(result$sch.new,list(sch.new))
    result$idx_sch_new <- c(result$idx_sch_new,list(idx_sch_new))
    
    result$Loss <- c(result$Loss,list(Loss))  ##Loss trajectory of iterations
    result$NP <- c(result$NP,list(NP)) ##Number of scheduled patients of iterations
    
    result$ci95 <- final_cost_ci95(obj.new)
    result5[[start]] <- result
  }
  result5
}


local_search_only_punc_erlang <- function(T,slot_length,mu,m,q0,q_breaks=NULL,ci,cw,co,cw_walkin=0.075,
                                          omega=0,dist,SetSeed=FALSE,nonpunc_prob = c(1/3,1/3,1/3),K=3,start_type="random",B0,
                                          Round=FALSE,model="erlang",cl=NULL){
  result5 <- list()
  for (start in 1:5) {
    set.seed(start*1000) 
    result <- list()
    #Read parameters from files
    # PAR <- SET[[case]][[1]]
    ###Set up beta-binomial distribution
    
    w <- c(
      ci, cw, co, cw_walkin
    )
    
    cov <- 1/sqrt(K) # Erlang coefficient of variation
    eta <- K/mu
    
    aa <- ParBetaBin(mu,cov,m)
    alpha <- aa$alpha
    beta <- aa$beta
    ####
    #bb <- rbbinom(n=10000, size=m, alpha, beta)
    #mean(bb)
    #sd(bb)/mean(bb)
    
    smp <- 0:m
    prob <- dbbinom(smp, m, alpha , beta , log = FALSE)
    
    #
    grid <- seq(0, T-slot_length, by=slot_length)
    sch.set <- grid
    ngrid <- length(grid)
    
    if(length(q0)==1){
      q.set <- rep(q0, ngrid)
    } else if(length(q0)==ngrid && is.null(q_breaks)){
      q.set <- q0
    } else { # designed for 2 or 3 time varying time probability
      if (is.null(q_breaks)) {
        stop("q_breaks must be provided for period-specific probabilities")
      }
      
      if (length(q_breaks) != length(q0) - 1L) {
        stop("length(q_breaks) must equal length(q0) - 1")
      }
      
      if (any(q_breaks <= 0) ||
          any(q_breaks >= T) ||
          any(diff(q_breaks) <= 0)) {
        stop("q_breaks must be increasing and strictly between 0 and T")
      }
      
      q.ind <- findInterval(grid, q_breaks)+1
      q.set <- q0[q.ind]
    } 
    
    q.average <- mean(q.set)
    N0 <- max(1, round(T / ((1 - q.average) * mu)))
    if(start_type=="even"){
      ## Evenly distribute the patients over the session
      ideal_times <- (0:(N0 - 1)) * T / N0
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="random"){
      set.seed(start*1000)
      sch.opt0 <- sample(sch.set,N0)
    } else if(start_type=="front"){
      ## Evenly distribute the patients over the 1st-half session
      ideal_times <- sch.set[1:round(length(sch.set)/2)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="back"){
      ## Evenly distribute the patients over the 2nd-half session
      ideal_times <- sch.set[round(length(sch.set)/2+1):length(sch.set)]
      ## Round each time to the nearest appointment slot
      sch.opt0 <- slot_length * floor(ideal_times / slot_length + 0.5)
      sch.opt0 <- pmin(sch.opt0, T - slot_length)
    } else if(start_type=="dpp"){
      if(start!=0){
        #The initial L matrix
        L <- diag(ngrid)
        
        ##Find the initial schedule
        Lnew <- L
        L.decomposed <- decompose.matrix(Lnew)
        kdpp.sample <- k.dpp.sample(L.decomposed, k=floor(min(N0,length(sch.set))-2))
        sch.opt0 <- sch.set[kdpp.sample]
      }
    }
    
    sch.opt0 <- sort(sch.opt0)
    
    ## Convert the schedule to the number of appointments in each grid slot
    slot_index <- match(sch.opt0, grid)
    ## Patient-specific no-show probabilities
    q.opt0 <- q.set[slot_index]
    
    idx_sch_ini <- idx_sch_opt0 <- tabulate(
      slot_index,
      nbins = ngrid
    )
    
    ###Set up the no show probability and nonpunctuality
    # dist="three_point" 
    
    # cl <- makeCluster(50)
    # registerDoSNOW(cl)
    
    ###First compute the cost for the given 
    obj.opt0 <- obj.ini <-  PFmeasure_parallel(sch.opt0,dist="uniform",q=q.opt0,omega=0,T=T, Round=Round,model=model,nsim=10000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
    loss.min <- mean(obj.ini$V)
    print(loss.min)
    
    ##Update the L matrix
    idx_sch_new <- idx_sch_ini
    sch.new <- sch.ini <- rep(sch.set,idx_sch_new) 
    q.new <- q.ini <- rep(q.set,idx_sch_new)  
    
    # obj.path <-  WorkloadPath_parallel(sch.new,dist=dist,q=q.new,sch.set=sch.set,omega,T=T, Round=TRUE,nsim=10000,model="discrete",p=prob,y=smp,cl=cl,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
    # wpath <- obj.path$workload_path
    # aa <- t(wpath)%*%wpath
    # obj.aa <- eigen(aa)
    # P <- obj.aa$vectors
    # d <- abs(obj.aa$values)
    # d[d<=0] <- d[d<=0]+0.1*min(d[d>0])
    # aa <- P%*%diag(d)%*%t(P)
    # bb <- diag(sqrt(1/diag(aa)))
    # L <- bb%*%aa%*%bb
    
    # Iter <- 2
    # iter <- 1
    
    # while(iter<=Iter){
    
    
    # if(iter==1) { 
    #   try <- 5
    #   nsim.inner <- 100
    #   fold.max <- 1
    # }
    # 
    # if(iter==2){
    #   try=5
    #   nsim.inner <- 100
    #   fold.max <- 5
    # } 
    
    
    # fold <- 1
    
    # while(fold <=fold.max){
    # while(fold <=fold.max){  
    #   if(iter>1){
    Loss <- NULL
    NP <- NULL
    nsim.inner <- 100
    change1 <- 1
    try <- 5
    t.begin <- Sys.time()
    
    while (change1!=0) {
      change1 <- 0
      ###Perform random walk search
      obj.tmp <-  LocalPFmeasure_parallel(sch.new,sch.set,dist="uniform",q.set=q.set,omega=0,T=T, Round=Round,model=model,nsim=1000,p=prob,y=smp,cl=cl,nsim.inner = nsim.inner,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
      aa <- sort(obj.tmp$V,index=TRUE)$ix
      
      ####Finish search until fails "try" times
      l <- 1
      count <- 1
      while (l <= try && count <= length(aa)) {
        kdd_tmp <- obj.tmp$KddSample[aa[count],]
        sch_tmp <- rep(sch.set,kdd_tmp) 
        q_tmp <- rep(q.set,kdd_tmp) 
        
        obj.tmp2 <-  PFmeasure_parallel(sch_tmp,dist="uniform",q=q_tmp,omega=0,T=T, Round=Round,model=model,nsim=B0,p=prob,y=smp,cl=cl,w=w,SetSeed=SetSeed,slot_length=slot_length,nonpunc_prob=nonpunc_prob,eta=eta,K=K)
        loss.tmp <- mean(obj.tmp2$V)
        
        
        if(loss.tmp<loss.min) {
          sch.new <- sch_tmp
          q.new <- q_tmp
          loss.min <- loss.tmp
          l <- 1
          change1 <- change1 + 1
        }else{
          
          l <- l + 1
        }
        count <- count + 1
        print(c(start,length(sch.new),loss.tmp,loss.min))
      }
      Loss <- c(Loss,loss.min)
      NP <- c(NP,length(sch.new))
    }
    #   }
    # 
    # iter <- iter +1
    # }
    
    time <- as.numeric(difftime(Sys.time(), t.begin, units = "secs"))
    
    # idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new))
    slot_id_new <- match(sch.new, sch.set)
    if (anyNA(slot_id_new)) {
      stop("The final schedule contains times outside sch.set")
    }
    idx_sch_new <- tabulate(slot_id_new,nbins = ngrid)
    q.new <- rep(q.set,idx_sch_new)  
    obj.new <-  PFmeasure_parallel(sch.new,dist=dist,q=q.new,omega,T=T, Round=Round,model=model,nsim=100000,p=prob,y=smp,cl=cl,w=w,slot_length=slot_length,nonpunc_prob=nonpunc_prob,SetSeed=TRUE,eta=eta,K=K)
    print(mean(obj.new$V))
    print(round(sch.new))
    # stopCluster(cl)
    
    aN0 <- mean(obj.opt0$N)
    aOT0 <- mean(obj.opt0$OT)
    aWT0 <- mean(obj.opt0$WT)/(length(sch.opt0)*mean(1-q.opt0))
    aIT0 <- mean(obj.opt0$IT)
    aV0 <- mean(obj.opt0$V)
    
    aN <- mean(obj.new$N)
    aOT <- mean(obj.new$OT)
    aWT <- mean(obj.new$WT)/(length(sch.new)*mean(1-q.new))
    aIT <- mean(obj.new$IT)
    aV <- mean(obj.new$V)
    
    result$par <- list(
      T = T,
      slot_length = slot_length,
      mu = mu,
      cv = cov,
      K = K,
      eta = eta,
      service_support = c(0, Inf),
      auxiliary_beta_binomial_m = m,
      no_show_probability = q0,
      weights = w,
      nonpunc_dist = dist,
      omega = omega,
      model = model,
      B0 = B0,
      nonpunc_prob = nonpunc_prob
    )
    
    result$slot_length <- slot_length
    result$time <- c(result$time,time)
    
    result$sch.set <- c(result$sch.set,list(sch.set)) 
    # result$par <- PAR
    result$sch.opt0 <- sch.opt0
    result$perf0 <- rbind(result$perf0,c(aN0,aOT0,aWT0,aIT0,aV0))
    
    result$perfnew <- rbind(result$perfnew,c(aN,aOT,aWT,aIT,aV))
    
    result$sch.new <- c(result$sch.new,list(sch.new))
    result$idx_sch_new <- c(result$idx_sch_new,list(idx_sch_new))
    
    result$Loss <- c(result$Loss,list(Loss))  ##Loss trajectory of iterations
    result$NP <- c(result$NP,list(NP)) ##Number of scheduled patients of iterations
    
    result$ci95 <- final_cost_ci95(obj.new)
    result5[[start]] <- result
  }
  result5
}

# All model settings below match the supplied Erlang experiment, except for
# the requested sweeps, B0 = 1000, and the selected uniform arrival interval.
run_erlang_case <- function(mu, cw, dist, omega, cl) {
  if (!dist %in% c("uniform", "uniform_early", "uniform_symmetric")) {
    stop("Unknown uniform arrival mode")
  }
  common <- list(
    T = 480, slot_length = 5, mu = mu, m = 30,
    q0 = 0.1, q_breaks = NULL,
    ci = 1, cw = cw, co = 1, cw_walkin = 0,
    omega = omega, dist = dist,
    nonpunc_prob = c(1/2, 1/2), K = 3,
    start_type = "even", B0 = 1000,
    Round = FALSE, model = "erlang_real", cl = cl
  )
  list(
    settings = common[names(common) != "cl"],
    result_local_search_plus_dpp = do.call(local_search_plus_dpp_erlang, common),
    result_local_search_only = do.call(local_search_only_erlang, common),
    result_punctual_local_search = do.call(local_search_only_punc_erlang, common)
  )
}

# Study 1: interval -> cw -> result_all indexed by mu (5, ..., 14).
run_erlang_study_1 <- function(cl, checkpoint_file = NULL) {
  arrival_cases <- list(
    late_0_8 = list(dist = "uniform", omega = 8),
    late_0_15 = list(dist = "uniform", omega = 15)
  )
  cw_values <- c(0.1, 0.5, 0.9)
  mu_values <- 5:14
  result_all_list <- if (!is.null(checkpoint_file) && file.exists(checkpoint_file))
    readRDS(checkpoint_file) else list()
  for (arrival_name in names(arrival_cases)) {
    if (is.null(result_all_list[[arrival_name]])) result_all_list[[arrival_name]] <- list()
    for (cw_value in cw_values) {
      cw_name <- sprintf("cw_%.1f", cw_value)
      if (is.null(result_all_list[[arrival_name]][[cw_name]]))
        result_all_list[[arrival_name]][[cw_name]] <- list()
      for (mu_value in mu_values) {
        mu_name <- paste0("mu_", mu_value)
        if (!is.null(result_all_list[[arrival_name]][[cw_name]][[mu_name]])) next
        result_all_list[[arrival_name]][[cw_name]][[mu_name]] <-
          run_erlang_case(mu = mu_value, cw = cw_value,
                          dist = arrival_cases[[arrival_name]]$dist,
                          omega = arrival_cases[[arrival_name]]$omega, cl = cl)
        if (!is.null(checkpoint_file)) saveRDS(result_all_list, checkpoint_file)
      }
    }
  }
  result_all_list
}

# Study 2: interval -> mu -> result_all indexed by cw (0.1, ..., 1.0).
run_erlang_study_2 <- function(cl, checkpoint_file = NULL) {
  arrival_cases <- list(
    late_0_8 = list(dist = "uniform", omega = 8),
    late_0_15 = list(dist = "uniform", omega = 15)
  )
  mu_values <- c(8, 10, 12)
  cw_values <- seq(0.1, 1.0, by = 0.1)
  result_all_list <- if (!is.null(checkpoint_file) && file.exists(checkpoint_file))
    readRDS(checkpoint_file) else list()
  for (arrival_name in names(arrival_cases)) {
    if (is.null(result_all_list[[arrival_name]])) result_all_list[[arrival_name]] <- list()
    for (mu_value in mu_values) {
      mu_name <- paste0("mu_", mu_value)
      if (is.null(result_all_list[[arrival_name]][[mu_name]]))
        result_all_list[[arrival_name]][[mu_name]] <- list()
      for (cw_value in cw_values) {
        cw_name <- sprintf("cw_%.1f", cw_value)
        if (!is.null(result_all_list[[arrival_name]][[mu_name]][[cw_name]])) next
        result_all_list[[arrival_name]][[mu_name]][[cw_name]] <-
          run_erlang_case(mu = mu_value, cw = cw_value,
                          dist = arrival_cases[[arrival_name]]$dist,
                          omega = arrival_cases[[arrival_name]]$omega, cl = cl)
        if (!is.null(checkpoint_file)) saveRDS(result_all_list, checkpoint_file)
      }
    }
  }
  result_all_list
}

# Study 3: interval -> result_all indexed by mu (5, ..., 14); cw = 0.1.
run_erlang_study_3 <- function(cl, checkpoint_file = NULL) {
  arrival_cases <- list(
    late_0_8 = list(dist = "uniform", omega = 8),
    early_minus8_0 = list(dist = "uniform_early", omega = 8),
    symmetric_minus8_8 = list(dist = "uniform_symmetric", omega = 8)
  )
  mu_values <- 5:14
  result_all_list <- if (!is.null(checkpoint_file) && file.exists(checkpoint_file))
    readRDS(checkpoint_file) else list()
  for (arrival_name in names(arrival_cases)) {
    if (is.null(result_all_list[[arrival_name]])) result_all_list[[arrival_name]] <- list()
    for (mu_value in mu_values) {
      mu_name <- paste0("mu_", mu_value)
      if (!is.null(result_all_list[[arrival_name]][[mu_name]])) next
      result_all_list[[arrival_name]][[mu_name]] <-
        run_erlang_case(mu = mu_value, cw = 0.1,
                        dist = arrival_cases[[arrival_name]]$dist,
                        omega = arrival_cases[[arrival_name]]$omega, cl = cl)
      if (!is.null(checkpoint_file)) saveRDS(result_all_list, checkpoint_file)
    }
  }
  result_all_list
}
