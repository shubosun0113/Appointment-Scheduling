library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(RhpcBLASctl)
library(Rcpp)
library(RcppArmadillo)
library(distributions3)

Rcpp::sourceCpp("workload_discrete.cpp")

###Given mean and Coefficent of variation (sigma/mu), find the beta binomial parameters
ParBetaBin <- function(mu,cov,m){
  
  ft <- function(par){
    alpha <- par[1]
    beta <- par[2]
    (mu-m*alpha/(alpha+beta))^2+(cov^2*mu^2-m*alpha*beta*(m+alpha+beta)/(alpha+beta)^2/(alpha+beta+1))^2
  }
  
  obj <- optim(c(0.5,0.5),ft)
  parhat <- obj$par
  list(alpha=parhat[1],beta=parhat[2],fn=obj$value)
}



###Generate arrival times from schedule
GenArrival <- function(sch,dist="uniform",q,omega,T, Round=FALSE,seed=NULL,q.u=NULL,u.set=NULL,
                       slot_length, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  if(!is.null(seed)) set.seed(seed)
  n <- length(sch)
  sch.sub <- sch[which(runif(n)>q)]
  n.sub <- length(sch.sub)
  arrival <- numeric(0)
  
  if (dist == "two_point") {
    gen.rand <- function(N, par) {
      prob2 <- nonpunc_prob / sum(nonpunc_prob)
      sample(
        x = c(0, 5),
        size = N,
        replace = TRUE,
        prob = prob2
      )
    }
  }
  
  if (dist == "three_point") {
    gen.rand <- function(N, par) {
      prob3 <- nonpunc_prob / sum(nonpunc_prob)
      sample(
        x = c(-slot_length, 0, slot_length),
        size = N,
        replace = TRUE,
        prob = prob3
      )
    }
  }
  
  # The uniform variants share omega but use different offset intervals.
  if(dist %in% c("uniform", "uniform_early", "uniform_symmetric")) gen.rand <- function(N,par){
    runif(N,par[1],par[2])
  } 
  
  if(dist=="asyuniform") gen.rand <- function(N,par){
    abs(runif(N,par[1],par[2]))
  }
  
  if(dist=="normal") gen.rand <- function(N,par){
    tmp <- rnorm(N,0,omega/3)
    idx <- tmp<par[1]
    tmp[idx] <- par[1]
    idx <- tmp>par[2]
    tmp[idx] <- par[2]
    tmp
  } 
  
  if(n.sub!=0){
    # arrival <- sch.sub[1:round(n.sub/2)]+gen.rand(round(n.sub/2),c(-omega,0))
    # arrival <- c(arrival, sch.sub[(round(n.sub/2)+1):n.sub]+gen.rand(n.sub-round(n.sub/2),c(0,omega)))
    offset_bounds <- c(0, omega)
    if (dist == "uniform_early") offset_bounds <- c(-omega, 0)
    if (dist == "uniform_symmetric") offset_bounds <- c(-omega, omega)
    arrival <- sch.sub + gen.rand(n.sub, offset_bounds)
    arrival <- arrival[arrival<=T]
    idx <- which(arrival<0)
    if(length(idx)!=0) arrival[idx] <- 0
    
    # if (dist == "uniform") {
    #   
    #   # Classify by scheduled appointment time.
    #   # At T = 480, the dividing time is minute 240.
    #   before_midpoint <- sch.sub < T / 2
    #   
    #   offset <- numeric(n.sub)
    #   
    #   # Earlier appointments: late arrivals
    #   offset[before_midpoint] <- runif(
    #     sum(before_midpoint),
    #     min = 0,
    #     max = omega
    #   )
    #   
    #   # Later appointments: early arrivals
    #   offset[!before_midpoint] <- runif(
    #     sum(!before_midpoint),
    #     min = -omega,
    #     max = 0
    #   )
    #   
    # } else {
    #   
    #   # Existing behavior for the other distributions
    #   offset <- gen.rand(n.sub, c(-omega, omega))
    #   
    # }
    
    # arrival <- sch.sub + offset
    
    arrival <- arrival[arrival <= T]
    arrival[arrival < 0] <- 0
  }
  
  
  idx.walkins <- rep(0,length(arrival))
  
  if(!is.null(q.u)){
    walkins <- u.set[runif(length(q.u))<q.u]+10^(-9) ###If walk ins and scheduled patients show up at the same time, scheduled patient receives service first
    arrival <- c(arrival,walkins)  
    idx.walkins <- c(idx.walkins,rep(1,length(walkins)))
  }
  
  idx <- sort(arrival,index=TRUE)$ix
  arrival <- arrival[idx]
  idx.walkins <- idx.walkins[idx]
  
  if(Round) arrival <- round(arrival)
  
  list(arrival=arrival,idx.walkins=idx.walkins)
}


###Given the schedule, compute the work load path: parallel version
WorkloadPath_parallel <- function(sch,dist="uniform",q,sch.set,omega,T, Round=FALSE,nsim=1,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,
                                  slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  close <- FALSE
  if(is.null(cl)){
    close <- TRUE
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
    
    RhpcBLASctl::blas_set_num_threads(1)
    
    source("Main.R")
    if(SetSeed)     arrival <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = i,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)$arrival else
      arrival <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)$arrival
    
    N <- length(arrival)
    
    if(N>0){
      
      
      if(model=="discrete")  {
        #          obj=Workload_discrete(arrival,p,y,T)
        obj=wldpath_discrete(arrival,p,y,T,grid=sch.set,Nmax=2000)
        m.s <- sum(y*p)
      }
      
      if(model=="erlang")  {
        obj=wldpath_erlang(arrival, K,  eta,  T,grid=sch.set,earlystop = 1)
        m.s <- K/eta
        if(obj$flag==1){
          N.R <- 100
          X <- Erlang(K, eta)
          R  <- matrix(random(X, N*N.R),nrow=N.R)
          obj=wldpath_real(arrival,R, T, grid=sch.set)
        }
        
      }
      
      if(model=="erlang_real")  {
        N.R <- 100
        X <- Erlang(K, eta)
        R  <- matrix(random(X, N*N.R),nrow=N.R)
        obj=wldpath_real(arrival,R, T, grid=sch.set)
        m.s <- K/eta
      }
      
      Z <- c(obj$workload_path)
      
    }else{
      
      Z <- sch.set*0
    }
    
    c(Z)  
  }
  #Stop the cluster
  close(pb)
  if(close){
    stopCluster(cl)
  }
  
  list(workload_path=result,sch.set=sch.set)
  
}


###Compute performance measurements of a given schedule through simulation: parallel version
PFmeasure_parallel <- function(sch,dist="uniform",q,omega,T, Round=FALSE,nsim=1,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,w=c(1,0.1,1,0.075),SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,
                               slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  close <- FALSE
  if(is.null(cl)){
    close <- TRUE
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
    
    RhpcBLASctl::blas_set_num_threads(1)
    
    source("Main.R")
    if(SetSeed)     obj <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = i,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob) else
      obj <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
    
    arrival <- obj$arrival
    idx.walkins <- obj$idx.walkins
    
    N <- length(arrival)
    
    if(N>0){
      
      if(model=="discrete")  {
        #          obj=Workload_discrete(arrival,p,y,T)
        obj=wld_discrete(arrival,p,y,T,Nmax=2000)
        m.s <- sum(y*p)
      }
      
      
      if(model=="erlang")  {
        obj=wld_erlang(arrival, K,  eta,  T,earlystop = 1)
        m.s <- K/eta
        if(obj$flag==1){
          N.R <- 100
          X <- Erlang(K, eta)
          R  <- matrix(random(X, N*N.R),nrow=N.R)
          obj=wld_real(arrival,R, T)
        }
        
      }
      
      if(model=="erlang_real")  {
        N.R <- 100
        X <- Erlang(K, eta)
        R  <- matrix(random(X, N*N.R),nrow=N.R)
        obj=wld_real(arrival,R, T)
        m.s <- K/eta
      }
      
      
      Z <- c(obj$workload)
      
      OT <- obj$zT
      WT <- sum((Z-m.s)[idx.walkins==0])
      WTu <- sum((Z-m.s)[idx.walkins==1])
      
      IT <- OT+T-N*m.s
      
      AWT <- (WT+WTu)/N
    }else{
      
      IT <- T
      OT <- 0
      AWT= WT=WTu <- 0
      N<-0
    }
    
    c(OT,IT,WT,AWT,N,WTu)  
  }
  #Stop the cluster
  close(pb)
  if(close){
    stopCluster(cl)
  }
  OT<- result[,1]
  IT<- result[,2]
  WT<- result[,3]
  AWT<- result[,4]
  N<- result[,5]
  WTu<- result[,6]
  
  
  V <- OT*w[3]+IT*w[1]+WT*w[2]+WTu*w[4]
  
  list(aV=mean(V),aOT=mean(OT),aWT=mean(AWT),aIT=mean(IT),OT=OT,IT=IT,WT=WT,WTu=WTu,V=V,N=N)
  
}



###Find the optimal sch through k-dpp sampling: parallel version
OptPFmeasure_parallel <- function(L,sch.new,sch.set, dist="uniform",q.set,omega,T, Round=FALSE,nsim=1,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,nsim.inner=1000,fix.nslots=FALSE,w=c(1,0.1,1,0.075),SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,MaxIter=2,
                                  slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  
  ###First remove all schedued appointments by one
  if(length(sch.new)>1)
    idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new)) else  idx_sch_new <- as.numeric(sapply(sch.set, `==`, sch.new)) 
    
    idx_sch_new_keep <- idx_sch_new
    
    idx <- which(idx_sch_new!=0)
    idx_sch_new_keep[idx] <- idx_sch_new[idx]-1
    
    ##determine the number of potions to be replaced
    k <- sum(idx_sch_new) - sum(idx_sch_new_keep)
    
    L.decomposed <- decompose.matrix(L)
    L.decomposed$lambda <- abs(L.decomposed$lambda)  ##Make sure L is positive definite
    
    close <- FALSE
    if(is.null(cl)){
      close <- TRUE
      cl <- makeCluster(numCores)
      registerDoSNOW(cl)
    }
    pb <- txtProgressBar(max = nsim, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    
    result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
      #  result <- foreach(i = 1:nsim,.options.snow = opts)%dopar%{
      
      RhpcBLASctl::blas_set_num_threads(1)
      source("helper-functions.R")
      source("kdpp-sample.R")
      source("Main.R")
      if(fix.nslots)  tmp <- k.dpp.sample(L.decomposed, k=k) else{
        r <- runif(1)
        if(r<=1/3) tmp <- k.dpp.sample(L.decomposed, k=max(k-1,1))
        if(r>1/3&r<=2/3) tmp <- k.dpp.sample(L.decomposed, k=k)
        if(r>2/3) tmp <- k.dpp.sample(L.decomposed, k=min(k+1,length(sch.set)))
      }
      
      idx_add <- sch.set*0
      idx_add[tmp] <- 1
      
      idx_sch_new_keep+idx_add
    }
    
    result <- rbind(idx_sch_new,result)
    
    idx <- !duplicated(result)
    KddSample <- result[idx,]
    KddSample <- KddSample[-1,,drop=F] ###Remove the current schedule
    
    nsim <- dim(KddSample)[1]
    for(iter in 1:MaxIter){
      
      if(iter==2){
        nsim.inner <- floor(nsim.inner*10)
        nsim <- min(max(100,floor(nsim/10)),nsim)
        aa <- sort(V,index=TRUE)$ix
        KddSample <- KddSample[aa[1:nsim],,drop=F]
      }
      
      result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
        
        RhpcBLASctl::blas_set_num_threads(1)
        
        source("helper-functions.R")
        source("kdpp-sample.R")
        source("Main.R")
        kdpp.sample <- KddSample[i,]
        
        sch <- rep(sch.set,kdpp.sample)
        q <- rep(q.set,kdpp.sample)
        
        N <- OT <- WT<-WTu <- IT <- AWT <- numeric(nsim.inner)
        for (j in 1:nsim.inner) {
          
          if(SetSeed)     obj <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = j,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob) else
            obj <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
          
          arrival <- obj$arrival
          idx.walkins <- obj$idx.walkins
          
          
          N[j] <- length(arrival)
          
          if(N[j]>0){
            
            if(model=="discrete")  {
              #          obj=Workload_discrete(arrival,p,y,T)
              obj=wld_discrete(arrival,p,y,T,Nmax=2000)
              m.s <- sum(y*p)
            }
            
            if(model=="erlang")  {
              obj=wld_erlang(arrival, K,  eta,  T,earlystop = 1)
              m.s <- K/eta
              if(obj$flag==1){
                N.R <- 100
                X <- Erlang(K, eta)
                R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
                obj=wld_real(arrival,R, T)
              }
              
            }
            
            if(model=="erlang_real")  {
              N.R <- 100
              X <- Erlang(K, eta)
              R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
              obj=wld_real(arrival,R, T)
              m.s <- K/eta
            }
            
            
            Z <- c(obj$workload)
            
            OT[j] <- obj$zT
            
            WT[j] <- sum((Z-m.s)[idx.walkins==0])
            WTu[j] <- sum((Z-m.s)[idx.walkins==1])
            
            IT[j] <- OT[j]+T-N[j]*m.s
            
            AWT[j] <- (WT[j]+WTu[j])/N[j]
          }else{
            
            IT[j] <- T
            OT[j] <- 0
            AWT[j]= WT[j] =WTu[j]<- 0
            N[j]<-0
          }
          
        }
        c(mean(OT),mean(IT),mean(WT),mean(AWT),mean(N),mean(WTu))  
        
      }
      
      #Stop the cluster
      close(pb)
      if(close){
        stopCluster(cl)
      }
      
      OT<- result[,1]
      IT<- result[,2]
      WT<- result[,3]
      AWT<- result[,4]
      N<- result[,5]
      WTu<- result[,6]
      
      V <- OT*w[3]+IT*w[1]+WT*w[2]+WTu*w[4]
    }
    idx.min <- which.min(V)
    
    kdd_opt <-  KddSample[idx.min,] 
    sch_opt <- rep(sch.set,kdd_opt)
    
    V_opt <- V[idx.min]
    OT_opt <- OT[idx.min]
    WT_opt <- WT[idx.min]
    WTu_opt <- WTu[idx.min]
    AWT_opt <- AWT[idx.min]
    IT_opt <- IT[idx.min]
    N_opt <- N[idx.min]
    
    list(V_opt=V_opt,OT_opt=OT_opt,WT_opt=WT_opt,WTu_opt=WTu_opt,IT_opt=IT_opt,AWT_opt=AWT_opt,sch_opt=sch_opt,OT=OT,IT=IT,WT=WT,WTu=WTu,V=V,N=N,kdd_opt=kdd_opt,KddSample=KddSample)
}

###Find the optimal sch through RandomWalk of the new schedule: parallel version
RwPFmeasure_parallel <- function(sch.new,sch.set, dist="uniform",q.set,omega,T, Round=FALSE,nsim=1000,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,nsim.inner=1000,w=c(1,0.1,1,0.075),SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,
                                 slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  
  k <- length(sch.new)
  ngrid <- length(sch.set)
  
  close <- FALSE
  if(is.null(cl)){
    close <- TRUE
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }
  pb <- txtProgressBar(max = nsim, style = 3)
  
  
  
  idx.sch.new <- sch.set*0
  idx.sch.new[which(sch.set%in%sch.new)] <- 1
  
  Result <- NULL
  
  
  ###First find all places that can perform random walk 
  set.emp <- setdiff(sch.set,sch.new)
  
  rw.set <- NULL
  for(i in 1:k){
    if(i==1) tmp <- c(set.emp[which(set.emp<sch.new[i])],  set.emp[which((set.emp>sch.new[i])&(set.emp<sch.new[i+1]))]) else
      if(i==k) tmp <- c(set.emp[which((set.emp>sch.new[i-1])&(set.emp<sch.new[i]))],  set.emp[which(set.emp>sch.new[i])]) else   
        tmp <- c(set.emp[which((set.emp>sch.new[i-1])&(set.emp<sch.new[i]))],  set.emp[which((set.emp>sch.new[i])&(set.emp<sch.new[i+1]))])
      if(length(tmp)!=0) rw.set <- c(rw.set,i)
  }
  
  
  n.rw <- length(rw.set)
  
  ### find all random walks upto order 2
  for (j in 1:n.rw) {
    for (l in j:n.rw) {
      
      idxj <-  which(sch.set%in%sch.new[rw.set[j]])
      idxl <-  which(sch.set%in%sch.new[rw.set[l]])
      
      
      if(j==l){
        idx.sch.tmp <- idx.sch.new
        idx.sch.tmp[idxj] <- idx.sch.tmp[idxj]-1
        if(idxj+1<=ngrid){
          idx.sch.tmp[idxj+1] <- idx.sch.tmp[idxj+1]+1
          Result <- rbind(Result,idx.sch.tmp)
        }
        idx.sch.tmp <- idx.sch.new
        idx.sch.tmp[idxj] <- idx.sch.tmp[idxj]-1
        if(idxj-1>=1){
          idx.sch.tmp[idxj-1] <- idx.sch.tmp[idxj-1]+1
          Result <- rbind(Result,idx.sch.tmp)
        }
      }else{
        
        idx.sch.tmp <- idx.sch.new
        idx.sch.tmp[idxj] <- idx.sch.tmp[idxj]-1
        if(idxj+1<=ngrid){
          idx.sch.tmp[idxj+1] <- idx.sch.tmp[idxj+1]+1
          
          idx.sch.tmp2 <- idx.sch.tmp
          idx.sch.tmp2[idxl] <- idx.sch.tmp2[idxl]-1
          
          if(idxl+1<=ngrid){
            
            idx.sch.tmp2[idxl+1] <- idx.sch.tmp2[idxl+1]+1
            Result <- rbind(Result,idx.sch.tmp2)
          }
          
          idx.sch.tmp2 <- idx.sch.tmp
          idx.sch.tmp2[idxl] <- idx.sch.tmp2[idxl]-1
          
          if(idxl-1>=1){
            idx.sch.tmp2[idxl-1] <- idx.sch.tmp2[idxl-1]+1
            Result <- rbind(Result,idx.sch.tmp2)
          }
        }
        
        
        
        idx.sch.tmp <- idx.sch.new
        idx.sch.tmp[idxj] <- idx.sch.tmp[idxj]-1
        if(idxj-1>=1){
          idx.sch.tmp[idxj-1] <- idx.sch.tmp[idxj-1]+1
          
          idx.sch.tmp2 <- idx.sch.tmp
          idx.sch.tmp2[idxl] <- idx.sch.tmp2[idxl]-1
          
          if(idxl+1<=ngrid){
            
            idx.sch.tmp2[idxl+1] <- idx.sch.tmp2[idxl+1]+1
            Result <- rbind(Result,idx.sch.tmp2)
          }
          
          idx.sch.tmp2 <- idx.sch.tmp
          idx.sch.tmp2[idxl] <- idx.sch.tmp2[idxl]-1
          
          if(idxl-1>=1){
            idx.sch.tmp2[idxl-1] <- idx.sch.tmp2[idxl-1]+1
            Result <- rbind(Result,idx.sch.tmp2)
          }
        }
        
        
      } 
      
    }
  }
  
  
  ### Find all places that can perform swap
  set.emp <- setdiff(sch.set,sch.new)
  n1 <- length(sch.new)
  n2 <- length(set.emp)
  
  ### find all random walks upto order 2
  for (j in 1:n1) {
    for (l in 1:n2) {
      idx.sch.tmp <- idx.sch.new
      
      idxj <- which(sch.set==sch.new[j])
      idxl <- which(sch.set==set.emp[l])
      
      idx.sch.tmp[idxl] <- idx.sch.new[idxj]
      idx.sch.tmp[idxj] <- idx.sch.new[idxl]
      
      Result <- rbind(Result,idx.sch.tmp)
    }     
  }
  
  idx <- rowSums(Result>1)==0
  
  Result <- Result[idx,]
  Result <- unique(Result)
  
  
  ft <- function(x){
    round(sch.set[x==1])
  }
  
  tmp <- t(apply(Result, 1, ft))
  
  idx <- !duplicated(tmp)
  Result <- Result[idx,]
  
  
  if(dim(Result)[1]>nsim) {
    idx <- sample(1:dim(Result)[1],size=nsim)
    Result <- Result[idx,]
  }
  
  ####Add or Delete a location from the schedule
  
  for (i in 1:ngrid) {
    tmp <- idx.sch.new
    if(idx.sch.new[i]==0){
      tmp[i] <- 1
      Result <- rbind(Result,tmp)
    }
    
    if(idx.sch.new[i]==1){
      tmp[i] <- 0
      Result <- rbind(Result,tmp)
    }
  }
  
  Result <- as.matrix(Result)
  
  tmp <- t(apply(Result, 1, ft))
  
  idx <- !duplicated(tmp)
  Result <- Result[idx,]
  
  nsim <- dim(Result)[1]
  
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  for(iter in 1:2){
    
    
    if(iter==2){
      
      nsim.inner <- floor(nsim.inner*10)
      nsim <- min(max(floor(nsim/10),100),nsim)
      aa <- sort(V,index=TRUE)$ix
      
      Result <- Result[aa[1:nsim],,drop=F]
    }
    
    
    result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
      
      RhpcBLASctl::blas_set_num_threads(1)
      
      source("helper-functions.R")
      source("kdpp-sample.R")
      source("Main.R")
      kdpp.sample <- which(Result[i,]==1)
      sch <- sch.set[kdpp.sample]
      q <- q.set[kdpp.sample]
      N <- OT <- WT <- IT <- WTu <- AWT <- numeric(nsim.inner)
      for (j in 1:nsim.inner) {
        
        if(SetSeed)     obj <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = j,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob) else
          obj <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
        
        
        arrival <- obj$arrival
        idx.walkins <- obj$idx.walkins
        
        
        N[j] <- length(arrival)
        
        if(N[j]>0){
          
          if(model=="discrete")  {
            #          obj=Workload_discrete(arrival,p,y,T)
            obj=wld_discrete(arrival,p,y,T,Nmax=2000)
            m.s <- sum(y*p)
          }
          
          if(model=="erlang")  {
            obj=wld_erlang(arrival, K,  eta,  T,earlystop = 1)
            m.s <- K/eta
            if(obj$flag==1){
              N.R <- 100
              X <- Erlang(K, eta)
              R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
              obj=wld_real(arrival,R, T)
            }
            
          }
          
          if(model=="erlang_real")  {
            N.R <- 100
            X <- Erlang(K, eta)
            R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
            obj=wld_real(arrival,R, T)
            m.s <- K/eta
          }
          
          
          Z <- c(obj$workload)
          
          OT[j] <- obj$zT
          WT[j] <- sum((Z-m.s)[idx.walkins==0])
          WTu[j] <- sum((Z-m.s)[idx.walkins==1])
          IT[j] <- OT[j]+T-N[j]*m.s
          
          AWT[j] <- (WT[j]+WTu[j])/N[j]
        }else{
          
          IT[j] <- T
          OT[j] <- 0
          AWT[j]= WT[j]=WTu[j] <- 0
          N[j]<-0
        }
        
      }
      c(mean(OT),mean(IT),mean(WT),mean(AWT),mean(N),mean(WTu))  
      
    }
    
    #Stop the cluster
    close(pb)
    if(close){
      stopCluster(cl)
    }
    
    OT<- result[,1]
    IT<- result[,2]
    WT<- result[,3]
    WTu<- result[,6]
    AWT<- result[,4]
    N<- result[,5]
    
    V <- OT*w[3]+IT*w[1]+WT*w[2]+WTu*w[4]
  }
  idx.min <- which.min(V)
  kdd_opt <-  which(Result[idx.min,]==1) 
  sch_opt <- sch.set[kdd_opt]
  
  V_opt <- V[idx.min]
  OT_opt <- OT[idx.min]
  WT_opt <- WT[idx.min]
  WTu_opt <- WTu[idx.min]
  AWT_opt <- AWT[idx.min]
  IT_opt <- IT[idx.min]
  N_opt <- N[idx.min]
  
  list(V_opt=V_opt,OT_opt=OT_opt,WT_opt=WT_opt,WTu_opt=WTu_opt,IT_opt=IT_opt,AWT_opt=AWT_opt,sch_opt=sch_opt,OT=OT,IT=IT,WT=WT,WTu=WTu,V=V,N=N,kdd_opt=kdd_opt,KddSample=Result)
}


###Find the optimal sch through swapping of the new schedule: parallel version
SwapPFmeasure_parallel <- function(sch.new,sch.set, dist="uniform",q.set,omega,T, Round=FALSE,nsim=1000,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,nsim.inner=1000,w=c(1,0.1,1,0.075),SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,
                                   slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  
  k <- length(sch.new)
  ngrid <- length(sch.set)
  
  close <- FALSE
  if(is.null(cl)){
    close <- TRUE
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }
  pb <- txtProgressBar(max = nsim, style = 3)
  
  
  
  idx.sch.new <- sch.set*0
  idx.sch.new[which(sch.set%in%sch.new)] <- 1
  
  Result <- NULL
  
  ### Find all places that can perform swap
  set.emp <- setdiff(sch.set,sch.new)
  n1 <- length(sch.new)
  n2 <- length(set.emp)
  
  ### find all random walks upto order 2
  for (j in 1:n1) {
    for (l in 1:n2) {
      idx.sch.tmp <- idx.sch.new
      
      idxj <- which(sch.set==sch.new[j])
      idxl <- which(sch.set==set.emp[l])
      
      idx.sch.tmp[idxl] <- idx.sch.new[idxj]
      idx.sch.tmp[idxj] <- idx.sch.new[idxl]
      
      Result <- rbind(Result,idx.sch.tmp)
    }     
  }
  
  idx <- rowSums(Result>1)==0
  
  Result <- Result[idx,]
  Result <- unique(Result)
  
  
  ft <- function(x){
    round(sch.set[x==1])
  }
  
  tmp <- t(apply(Result, 1, ft))
  
  idx <- !duplicated(tmp)
  Result <- Result[idx,]
  
  
  if(dim(Result)[1]>nsim) {
    idx <- sample(1:dim(Result)[1],size=nsim)
    Result <- Result[idx,]
  }
  
  ####Add or Delete a location from the schedule
  
  for (i in 1:ngrid) {
    tmp <- idx.sch.new
    if(idx.sch.new[i]==0){
      tmp[i] <- 1
      Result <- rbind(Result,tmp)
    }
    
    if(idx.sch.new[i]==1){
      tmp[i] <- 0
      Result <- rbind(Result,tmp)
    }
  }
  
  Result <- as.matrix(Result)
  
  tmp <- t(apply(Result, 1, ft))
  
  idx <- !duplicated(tmp)
  Result <- Result[idx,]
  
  nsim <- dim(Result)[1]
  
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  
  for(iter in 1:2){
    
    
    if(iter==2){
      
      nsim.inner <- floor(nsim.inner*10)
      nsim <- min(max(floor(nsim/10),100),nsim)
      aa <- sort(V,index=TRUE)$ix
      
      Result <- Result[aa[1:nsim],,drop=F]
    }
    
    
    result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
      
      RhpcBLASctl::blas_set_num_threads(1)
      
      source("helper-functions.R")
      source("kdpp-sample.R")
      source("Main.R")
      kdpp.sample <- which(Result[i,]==1)
      sch <- sch.set[kdpp.sample]
      q <- q.set[kdpp.sample]
      N <- OT <- WT <- IT <- WTu <- AWT <- numeric(nsim.inner)
      for (j in 1:nsim.inner) {
        
        if(SetSeed)     obj <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = j,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob) else
          obj <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
        
        
        arrival <- obj$arrival
        idx.walkins <- obj$idx.walkins
        
        
        N[j] <- length(arrival)
        
        if(N[j]>0){
          
          if(model=="discrete")  {
            #          obj=Workload_discrete(arrival,p,y,T)
            obj=wld_discrete(arrival,p,y,T,Nmax=2000)
            m.s <- sum(y*p)
          }
          
          if(model=="erlang")  {
            obj=wld_erlang(arrival, K,  eta,  T,earlystop = 1)
            m.s <- K/eta
            if(obj$flag==1){
              N.R <- 100
              X <- Erlang(K, eta)
              R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
              obj=wld_real(arrival,R, T)
            }
            
          }
          
          if(model=="erlang_real")  {
            N.R <- 100
            X <- Erlang(K, eta)
            R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
            obj=wld_real(arrival,R, T)
            m.s <- K/eta
          }
          
          
          Z <- c(obj$workload)
          
          OT[j] <- obj$zT
          WT[j] <- sum((Z-m.s)[idx.walkins==0])
          WTu[j] <- sum((Z-m.s)[idx.walkins==1])
          IT[j] <- OT[j]+T-N[j]*m.s
          
          AWT[j] <- (WT[j]+WTu[j])/N[j]
        }else{
          
          IT[j] <- T
          OT[j] <- 0
          AWT[j]= WT[j]=WTu[j] <- 0
          N[j]<-0
        }
        
      }
      c(mean(OT),mean(IT),mean(WT),mean(AWT),mean(N),mean(WTu))  
      
    }
    
    #Stop the cluster
    close(pb)
    if(close){
      stopCluster(cl)
    }
    
    OT<- result[,1]
    IT<- result[,2]
    WT<- result[,3]
    WTu<- result[,6]
    AWT<- result[,4]
    N<- result[,5]
    
    V <- OT*w[3]+IT*w[1]+WT*w[2]+WTu*w[4]
  }
  idx.min <- which.min(V)
  kdd_opt <-  which(Result[idx.min,]==1) 
  sch_opt <- sch.set[kdd_opt]
  
  V_opt <- V[idx.min]
  OT_opt <- OT[idx.min]
  WT_opt <- WT[idx.min]
  WTu_opt <- WTu[idx.min]
  AWT_opt <- AWT[idx.min]
  IT_opt <- IT[idx.min]
  N_opt <- N[idx.min]
  
  list(V_opt=V_opt,OT_opt=OT_opt,WT_opt=WT_opt,WTu_opt=WTu_opt,IT_opt=IT_opt,AWT_opt=AWT_opt,sch_opt=sch_opt,OT=OT,IT=IT,WT=WT,WTu=WTu,V=V,N=N,kdd_opt=kdd_opt,KddSample=Result)
}



###Find the optimal sch through swapping of the new schedule: parallel version
LocalPFmeasure_parallel <- function(sch.new,sch.set, dist="uniform",q.set,omega,T, Round=FALSE,nsim=1000,model="discrete",p=NULL,y=NULL,numCores=10,cl=NULL,nsim.inner=1000,w=c(1,0.1,1,0.075),SetSeed=FALSE,q.u=NULL,u.set=NULL,eta=NULL,K=NULL,
                                    slot_length = 15, nonpunc_prob = c(1/3, 1/3, 1/3)){
  
  
  k <- length(sch.new)
  ngrid <- length(sch.set)
  
  close <- FALSE
  if(is.null(cl)){
    close <- TRUE
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
  }
  pb <- txtProgressBar(max = nsim, style = 3)
  
  ##Find the index set of schedule
  if(length(sch.new)>1)
    idx_sch_new <- colSums(sapply(sch.set, `==`, sch.new)) else  idx_sch_new <- as.numeric(sapply(sch.set, `==`, sch.new))   
  
  Result <- NULL
  
  ####Change one location from the schedule
  
  for (i in 1:ngrid) {
    tmp <- idx_sch_new
    tmp[i] <- idx_sch_new[i]+1
    Result <- rbind(Result,tmp)
    
    if(idx_sch_new[i]!=0){
      tmp[i] <- idx_sch_new[i]-1
      Result <- rbind(Result,tmp)
    }
  }
  
  ##Change two loactions from the schedule
  ncoms <- combn(1:ngrid,2)
  
  progress <- function(count) setTxtProgressBar(pb, count)
  opts <- list(progress = progress)
  
  result <- foreach(i = 1:dim(ncoms)[2],.options.snow = opts,.combine = rbind)%dopar%{
    
    RhpcBLASctl::blas_set_num_threads(1)
    
    idx1 <-  ncoms[1,i]
    idx2 <-  ncoms[2,i]
    
    tmp2 <- NULL
    
    tmp <- idx_sch_new
    tmp[idx1] <- idx_sch_new[idx1]+1
    tmp[idx2] <- idx_sch_new[idx2]+1
    tmp2 <- rbind(tmp2,tmp)
    
    tmp <- idx_sch_new
    tmp[idx1] <- idx_sch_new[idx1]-1
    tmp[idx2] <- idx_sch_new[idx2]+1
    tmp2 <- rbind(tmp2,tmp)
    
    tmp <- idx_sch_new
    tmp[idx1] <- idx_sch_new[idx1]+1
    tmp[idx2] <- idx_sch_new[idx2]-1
    tmp2 <- rbind(tmp2,tmp)
    
    tmp <- idx_sch_new
    tmp[idx1] <- idx_sch_new[idx1]-1
    tmp[idx2] <- idx_sch_new[idx2]-1
    tmp2 <- rbind(tmp2,tmp)
    
    idx <- rowSums(tmp2<0) != 0
    if(sum(idx)!=0)
      tmp2<- tmp2[!idx,,drop=F]
    
    tmp2
  }
  
  Result <- rbind(Result,result)
  
  Result <- unique(Result)
  
  
  
  ### adjust 3 to 4 locations at a time
  if(dim(Result)[1]<nsim){
    # k.max <- 10
    k.max <- min(10,ngrid)
    k <- 3
    while (k <= k.max) {
      
      l <- rep(list(c(-1,1,0)), k)
      opk <- as.matrix(expand.grid(l))
      
      tmp <- t(as.matrix(idx_sch_new)%*%(rep(1,dim(opk)[1])))
      
      result <- foreach(i = 1:(ngrid-k+1),.options.snow = opts,.combine = rbind)%dopar%{
        
        RhpcBLASctl::blas_set_num_threads(1)
        idx.sub <- i:(i+k-1)
        set.sub <- idx_sch_new[idx.sub]
        tmp2 <- tmp
        for (j in 1:k) {
          tmp2[,idx.sub[j]] <- tmp[,idx.sub[j]]+opk[,j] 
        }
        
        idx <- rowSums(tmp2<0) != 0
        if(sum(idx)!=0)
          tmp2<- tmp2[!idx,,drop=F]
        
        tmp2
        
      }
      
      Result <- rbind(Result,result)
      
      Result <- unique(Result)
      if(dim(Result)[1]>nsim)  k <- k.max+1 else k <- k+1
    }
  }
  
  nsim <- dim(Result)[1]
  
  for(iter in 1:2){
    
    
    if(iter==2){
      
      nsim.inner <- floor(nsim.inner*10)
      nsim <- min(max(floor(nsim/10),100),nsim)
      aa <- sort(V,index=TRUE)$ix
      
      Result <- Result[aa[1:nsim],,drop=F]
    }
    
    
    result <- foreach(i = 1:nsim,.options.snow = opts,.combine = rbind)%dopar%{
      
      RhpcBLASctl::blas_set_num_threads(1)
      
      source("helper-functions.R")
      source("kdpp-sample.R")
      source("Main.R")
      
      sch <- rep(sch.set,Result[i,]) 
      q <- rep(q.set,Result[i,]) 
      
      N <- OT <- WT <- IT <- WTu <- AWT <- numeric(nsim.inner)
      for (j in 1:nsim.inner) {
        
        if(SetSeed)     obj <- GenArrival(sch,dist=dist,q,omega,T,Round,seed = j,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob) else
          obj <- GenArrival(sch,dist=dist,q,omega,T,Round,q.u=q.u,u.set=u.set,slot_length=slot_length,nonpunc_prob=nonpunc_prob)
        
        
        arrival <- obj$arrival
        idx.walkins <- obj$idx.walkins
        
        
        N[j] <- length(arrival)
        
        if(N[j]>0){
          
          if(model=="discrete")  {
            #          obj=Workload_discrete(arrival,p,y,T)
            obj=wld_discrete(arrival,p,y,T,Nmax=2000)
            m.s <- sum(y*p)
          }
          
          if(model=="erlang")  {
            obj=wld_erlang(arrival, K,  eta,  T,earlystop = 1)
            m.s <- K/eta
            if(obj$flag==1){
              N.R <- 100
              X <- Erlang(K, eta)
              R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
              obj=wld_real(arrival,R, T)
            }
            
          }
          
          if(model=="erlang_real")  {
            N.R <- 100
            X <- Erlang(K, eta)
            R  <- matrix(random(X, N[j]*N.R),nrow=N.R)
            obj=wld_real(arrival,R, T)
            m.s <- K/eta
          }
          
          
          Z <- c(obj$workload)
          
          OT[j] <- obj$zT
          WT[j] <- sum((Z-m.s)[idx.walkins==0])
          WTu[j] <- sum((Z-m.s)[idx.walkins==1])
          IT[j] <- OT[j]+T-N[j]*m.s
          
          AWT[j] <- (WT[j]+WTu[j])/N[j]
        }else{
          
          IT[j] <- T
          OT[j] <- 0
          AWT[j]= WT[j]=WTu[j] <- 0
          N[j]<-0
        }
        
      }
      c(mean(OT),mean(IT),mean(WT),mean(AWT),mean(N),mean(WTu))  
      
    }
    
    #Stop the cluster
    close(pb)
    if(close){
      stopCluster(cl)
    }
    
    OT<- result[,1]
    IT<- result[,2]
    WT<- result[,3]
    WTu<- result[,6]
    AWT<- result[,4]
    N<- result[,5]
    
    V <- OT*w[3]+IT*w[1]+WT*w[2]+WTu*w[4]
  }
  idx.min <- which.min(V)
  kdd_opt <-  Result[idx.min,]
  sch_opt <- rep(sch.set,kdd_opt) 
  
  V_opt <- V[idx.min]
  OT_opt <- OT[idx.min]
  WT_opt <- WT[idx.min]
  WTu_opt <- WTu[idx.min]
  AWT_opt <- AWT[idx.min]
  IT_opt <- IT[idx.min]
  N_opt <- N[idx.min]
  
  list(V_opt=V_opt,OT_opt=OT_opt,WT_opt=WT_opt,WTu_opt=WTu_opt,IT_opt=IT_opt,AWT_opt=AWT_opt,sch_opt=sch_opt,OT=OT,IT=IT,WT=WT,WTu=WTu,V=V,N=N,kdd_opt=kdd_opt,KddSample=Result)
}

##Compute the work load of real time

Workload_real <- function(arrival, R){
  N <- length(arrival)
  Z <- numeric(N)
  Z[1] <- R[1]
  if(N>1){
    for (i in 2:N) {
      Z[i]  <- max(Z[i-1]-(arrival[i]-arrival[i-1]),0)+R[i]   
    }
  }
  Z
}



##Compute the workload process for discrete distribution and discrete arrival time
## y are integers starts with 0

Workload_discrete <- function(arrival,p,y,T,path=FALSE,grid=NULL){
  N <- length(arrival)
  n.y <- length(y)
  if(N>1)
    z <- seq(0,max(diff(arrival))+10*max(y),by=1) else z <- seq(0,10*max(y),by=1)
  n.z <- length(z)
  
  g <- G <- matrix(0,N,n.z)
  
  g[1,1:n.y] <- p
  G[1,] <- 1
  G[1,1:n.y] <- cumsum(p)
  ind <- outer(z,y,FUN="-")>0
  
  if(N>1){
    for(i in 2:N){
      gap <- arrival[i]-arrival[i-1]
      g[i,1:n.y] <- G[i-1,gap+1]*p
      
      ft <- function(x){
        val <- numeric(length(x))
        idx <- which((x>=0)&(x<=n.z-1)) 
        
        if(length(idx)!=0){
          val[idx] <- g[i-1,x[idx]+1]
        }
        val 
      }
      
      
      tmp <- outer(z,y,FUN="-")+gap 
      
      idx <- (tmp>=0)&(tmp<=n.z-1) 
      tmp <- tmp*idx
      tmp[idx] <- ft(tmp[idx])
      tmp <- tmp*ind
      g[i,] <- g[i,]+c(tmp%*%p)
      G[i,] <- cumsum(g[i,])
    }
  }
  
  if(N>1) workload <- rowSums(1-G) else workload <- sum(1-G)
  
  gap <- T-arrival[N]
  if(gap<=n.z-1) zT <- sum(1-G[N,(1+gap):n.z]) else zT=0
  
  workload_path <- NULL
  if(path){
    ngrid <- length(grid)
    workload_path <- numeric(ngrid)
    
    for (i in 1:ngrid) {
      if(grid[i]%in%arrival){
        workload_path[i] <- workload[which(arrival==grid[i])]
      }else{
        
        if(sum(arrival<grid[i])==0) workload_path[i]=0 else{
          
          idx <- max(which(arrival<grid[i]))
          
          gap <- grid[i]-arrival[idx]
          if(gap<=n.z-1) workload_path[i] <- sum(1-G[idx,(1+gap):n.z]) else workload_path[i]=0  
        }
        
        
      }
    }
  }
  
  list(g=g,G=G,workload=workload,z=z,zT=zT,workload_path=workload_path)
}


##Compute the workload process for exponential distribution

Workload_exp <- function(arrival,eta){
  N <- length(arrival)
  alpha <-beta <- matrix(0,nrow=N,ncol=N)
  beta[1,1] <- 1
  if(N == 1){return(beta)}
  for(i in 2:N){
    for(j in 1:i){
      if(j==1){
        l <- 2:i
        alpha[i,j] <- 1-sum(beta[i-1,l-1]* exp(-lfactorial(l-1)+(l-1)*log(eta*(arrival[i]-arrival[i-1]))-eta*(arrival[i]-arrival[i-1])))
      }else{
        alpha[i,j] <- beta[i-1,j-1]*exp(-eta*(arrival[i]-arrival[i-1]))
      }
    }
    
    for(j in 1:i){
      l <- j:i
      beta[i,j] <- sum(alpha[i,l]* exp(-lfactorial(l-j)+(l-j)*log(eta*(arrival[i]-arrival[i-1]))))
      
    }
  }
  list(beta=beta,alpha=alpha,workload=rowSums(beta)/eta)
}

##Compute the workload process for Erlang(K,eta) distribution

Workload_erlang <- function(arrival,K,eta,T){
  N <- length(arrival)
  alpha <-beta <- matrix(0,nrow=N,ncol=N*K)
  beta[1,1:K] <- 1
  if(N>1){
    for(i in 2:N){
      for(j in 1:(i*K)){
        if(j<=K){
          l <- 1:((i-1)*K)
          alpha[i,j] <- -sum((-1)^(K-j)*beta[i-1,l]* exp(-lfactorial(l-1)+(K+l-j)*log(eta*(arrival[i]-arrival[i-1]))-lfactorial(K-j)-log(K+l-j)-eta*(arrival[i]-arrival[i-1])))
        }else{
          l <- 1:K
          alpha[i,j] <- beta[i-1,j-K]* sum((-1)^(K-l)*exp(-log(j-l)-lfactorial(l-1)-lfactorial(K-l)+lfactorial(j-1)-lfactorial(j-K-1)-eta*(arrival[i]-arrival[i-1])))
        }
      }
      
      for(j in 1:(i*K)){
        if(j<=K){
          l <- j:(i*K)
          beta[i,j] <- 1+sum(alpha[i,l]* exp(-lfactorial(l-j)+(l-j)*log(eta*(arrival[i]-arrival[i-1]))))
        }else{
          l <- j:(i*K)
          beta[i,j] <- sum(alpha[i,l]* exp(-lfactorial(l-j)+(l-j)*log(eta*(arrival[i]-arrival[i-1]))))
        }
      }
      
      ###quality control (remove second or latter peaks)
      idx <- min(which(beta[i,]< 0.001|c(-0.1,diff(beta[i,]))>1))
      
      idx1 <- which(abs(beta[i,idx:(i*K)])>0.05)
      
      if(length(idx1)!=0){
        print(i)
        idx2 <- which.min( abs(beta[i,1:(idx+min(idx1)-1)]))
        beta[i,1+idx2:(i*K-1)] <- 0
      }
    }
  }
  
  zT <- 0
  for (j in 1:(N*K)) {
    l <- 1:j
    zT <- zT+beta[N,j]*sum(exp((j-l)*log(eta*(T-arrival[N]))-eta*(T-arrival[N])-lfactorial(j-l)))/eta
  }
  
  list(beta=beta,alpha=alpha,workload=rowSums(beta)/eta,zT=zT)
}



Workload_erlang_tmp <- function(arrival,K,eta){
  N <- length(arrival)
  lgamma <-  alpha <-beta <- matrix(0,nrow=N,ncol=N*K)
  
  for (i in 2:N) {
    j <- 1:(i*K)
    #  lgamma[i,j] <- max(-j*log(eta*(arrival[i]-arrival[i-1])),-log(1000000000000))
    lgamma[i,j] <- -K*log(eta*(arrival[i]-arrival[i-1]))
  }
  lgamma <- lgamma*0+1
  
  beta[1,1:K] <- exp(-lgamma[1,1:K])
  if(N == 1){return(beta)}
  for(i in 2:N){
    for(j in 1:(i*K)){
      
      if(j<=K){
        l <- 1:((i-1)*K)
        alpha[i,j] <- -sum((-1)^(K-j)*beta[i-1,l]* exp(-lfactorial(l-1)+(K+l-j)*log(eta*(arrival[i]-arrival[i-1]))-lfactorial(K-j)-log(K+l-j)-eta*(arrival[i]-arrival[i-1])+lgamma[i-1,l]-lgamma[i,j]))
      }else{
        l <- 1:K
        alpha[i,j] <- beta[i-1,j-K]* sum((-1)^(K-l)*exp(-log(j-l)-lfactorial(l-1)-lfactorial(K-l)+lfactorial(j-1)-lfactorial(j-K-1)-eta*(arrival[i]-arrival[i-1])+lgamma[i-1,j-K]-lgamma[i,j]))
      }
    }
    
    for(j in 1:(i*K)){
      if(j<=K){
        l <- j:(i*K)
        beta[i,j] <- exp(-lgamma[i,j])+sum(alpha[i,l]* exp(-lfactorial(l-j)+(l-j)*log(eta*(arrival[i]-arrival[i-1]))+lgamma[i,l]-lgamma[i,j]))
      }else{
        l <- j:(i*K)
        beta[i,j] <- sum(alpha[i,l]* exp(-lfactorial(l-j)+(l-j)*log(eta*(arrival[i]-arrival[i-1]))+lgamma[i,l]-lgamma[i,j]))
      }
    }
  }
  
  beta <- beta*exp(lgamma)
  
  list(beta=beta,alpha=alpha,workload=rowSums(beta)/eta)
}


##Compute the workload process for Hyperexponential distribution

Workload_hyperexp <- function(arrival,K,p,eta){
  n <- length(arrival)
  
  alpha <- beta <- array(0,dim=c(n,n,K))
  beta[1,1,] <- p
  
  for(i in 2:n){
    for(k in 1:K){
      for(j in 1:i){
        if(j==1){
          l <- 2:i
          alpha[i,j,k] <- p[k]-p[k]*sum(beta[i-1,l-1,k]*exp(-lfactorial(l-1)-eta[k]*(arrival[i]-arrival[i-1])+(l-1)*log(eta[k]*(arrival[i]-arrival[i-1]))))
          for(l in 1:K){
            if(l != k){
              for(h in 2:i){
                m <- 0:(h-2)
                alpha[i,j,k]  <-   alpha[i,j,k] +p[k]*beta[i-1,h-1,l]*sum(exp((h-2)*log(eta[l])+log(eta[k])+(h-2-m)*log(arrival[i]-arrival[i-1])-eta[l]*(arrival[i]-arrival[i-1])-lfactorial(h-2-m)-(m+1)*log(abs(eta[l]-eta[k])))*(sign(eta[l]-eta[k])^(m+1)))
                alpha[i,j,k]  <-   alpha[i,j,k] -beta[i-1,h-1,k]*p[l]*exp(log(eta[l])+(h-2)*log(eta[k])-eta[k]*(arrival[i]-arrival[i-1])-(h-1)*log(abs(eta[k]-eta[l])))*(sign(eta[k]-eta[l])^(h-1))
              }
            }
          }
          
        }else if(j==i){
          alpha[i,j,k] <- beta[i-1,j-1,k]*p[k]*exp(-eta[k]*(arrival[i]-arrival[i-1]))
        }else{
          alpha[i,j,k] <- beta[i-1,j-1,k]*p[k]*exp(-eta[k]*(arrival[i]-arrival[i-1]))
          for(l in 1:K){
            if(l != k){
              h <- (j+1):i
              alpha[i,j,k]  <- alpha[i,j,k] -p[l]*sum(beta[i-1,h-1,k]*exp(log(eta[l])+(h-1-j)*log(eta[k])-eta[k]*(arrival[i]-arrival[i-1])-(h-j)*log(abs(eta[k]-eta[l])))*sign(eta[k]-eta[l])^(h-j))
            }
          }
        }
      }
    }
    
    for(k in 1:K){
      for(j in 1:i){
        l <- j:i
        
        beta[i,j,k] <- sum(alpha[i,l,k]*exp((l-j)*log(eta[k]*(arrival[i]-arrival[i-1]))-lfactorial(l-j)))
      }
    }   
  }
  workload <- c((apply(beta, c(1,3), sum)%*%(1/eta)))
  
  list(beta=beta,alpha=alpha,workload=workload)
  
}
