#include "RcppArmadillo.h"
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double lfact(int x){
double a=0;
  if(x>0){
  for(int i=0;i<x;i++){
    a=a+log(i+1);
  }
  }
  return a;
}
  
  
// [[Rcpp::export]]
Rcpp::List wld_discrete(arma::vec arrival,arma::vec p,arma::vec y, int T, int Nmax)
{
 
 int N=arrival.n_elem;
 int ny=y.n_elem;
 arma::mat G;
 G.ones(N,Nmax);

 G(0,0)= p(0);
 
 for(int i=1;i<ny;i++){
    G(0,i)=0;
   for(int j=0;j<=i;j++){
     G(0,i) =G(0,i)+ p(j);
   }
 }
 
 for(int i=1;i<N;i++){
  
  double tol=1;
   int z=1;
   int gap=arrival(i)-arrival(i-1);
     
   G(i,0)= G(i-1,gap)*p(0);
   
   while((tol>0.00001)&(z<Nmax)){
     
     G(i,z)=0;
     for(int j=0;j<(z+1);j++){
       if(j<ny)       G(i,z)=G(i,z)+G(i-1,z-j+gap)*p(j);
     }
     
     tol=1-G(i,z);
     z=z+1;
   }

 } 
 
  //find the workload
  arma::vec workload(N);
 
  for(int i=0;i<N;i++){
 workload(i)=sum(1-G.row(i));   
  }
  
  //find the workload at time T
  int gap = T-arrival(N-1);
  double zT=0;
  if(gap<=Nmax-1) {
    for(int i=gap;i<Nmax;i++){
      zT=zT+1-G.row(N-1)(i);
    }
  } 
  
    Rcpp::List out;
  out["G"] = G;
  out["workload"] = workload;
  out["zT"] = zT;
  return out;
}


// [[Rcpp::export]]
Rcpp::List wldpath_discrete(arma::vec arrival,arma::vec p,arma::vec y, int T, arma::vec grid, int Nmax)
{
  
  int N=arrival.n_elem;
  int ny=y.n_elem;
  arma::mat G;
  G.ones(N,Nmax);
  
  G(0,0)= p(0);
  
  for(int i=1;i<ny;i++){
    G(0,i)=0;
    for(int j=0;j<=i;j++){
      G(0,i) =G(0,i)+ p(j);
    }
  }
  
  for(int i=1;i<N;i++){
    
    double tol=1;
    int z=1;
    int gap=arrival(i)-arrival(i-1);
    
    G(i,0)= G(i-1,gap)*p(0);
    
    while((tol>0.00001)&(z<Nmax)){
      
      G(i,z)=0;
      for(int j=0;j<(z+1);j++){
        if(j<ny)       G(i,z)=G(i,z)+G(i-1,z-j+gap)*p(j);
      }
      
      tol=1-G(i,z);
      z=z+1;
    }
    
  } 
  
  //find the workload
  arma::vec workload(N);
  
  for(int i=0;i<N;i++){
    workload(i)=sum(1-G.row(i));   
  }
  
  //find the workload at time T
  int gap = T-arrival(N-1);
  double zT=0;
  if(gap<=Nmax-1) {
    for(int i=gap;i<Nmax;i++){
      zT=zT+1-G.row(N-1)(i);
    }
  } 
  
  //Find the workload path
    int ngrid=grid.n_elem;
    arma::vec workload_path(ngrid);
    
    int current=0;
    
    for(int i=0;i<ngrid;i++){
      
      if(grid(i)<arrival(0)) workload_path(i)=0; 
      if(current< (N-1)) {
        while((grid(i)>=arrival(current))&(grid(i)>= arrival(current+1))) {
          current=current +1; 
        if(current>=(N-1)) break;
        }
      }
      
       if(current< (N-1)){
          
        if((grid(i)>=arrival(current))&(grid(i)< arrival(current+1))) {
          gap = grid(i)-arrival(current);
          if(gap<=Nmax-1) {
            for(int j=gap;j<Nmax;j++){
              workload_path(i)=workload_path(i)+1-G.row(current)(j);
            }
          }
        }
        
      }
      
      if(grid(i)>=arrival(N-1)){
        gap = grid(i)-arrival(N-1);
        if(gap<=Nmax-1) {
          for(int j=gap;j<Nmax;j++){
            workload_path(i)=workload_path(i)+1-G.row(N-1)(j);
          }
        }
      } 
      
    }
  
  Rcpp::List out;
  out["G"] = G;
  out["workload"] = workload;
  out["zT"] = zT;
  out["workload_path"] = workload_path;
  return out;
}


// [[Rcpp::export]]
Rcpp::List wld_erlang(arma::vec arrival, int K, double eta, int T,int earlystop)
{
  
  int N=arrival.n_elem;
  arma::mat alpha;
  alpha.zeros(N,N*K);
  arma::mat beta;
  beta.zeros(N,N*K);
  double small=pow(10,-50);
  int flag=0;
  
  for(int j=0;j<K;j++){
      beta(0,j)=1;
  }
  
  if(N>2){
    
    for(int i=2;i<=N;i++){
      
      for(int j=1;j<=K;j++){
        for(int l=1;l<=((i-1)*K);l++){
          if(beta(i-2,l-1)!=0)
          alpha(i-1,j-1)=alpha(i-1,j-1)-beta(i-2,l-1)*exp(-lfact(l-1)+(K+l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small))-lfact(K-j)-log(K+l-j)-eta*(arrival(i-1)-arrival(i-2)))*pow(-1,K-j);
        }
        
      }
      
      for(int j=K+1;j<=(i*K);j++){
        
        for(int l=1;l<=K;l++){
          if(beta(i-2,j-1-K)!=0)
          alpha(i-1,j-1)=alpha(i-1,j-1)+beta(i-2,j-1-K)*exp(-eta*(arrival(i-1)-arrival(i-2))-log(j-l)-lfact(l-1)-lfact(K-l)+lfact(j-1)-lfact(j-K-1))*pow(-1,K-l);
        }
        
      }
      
      
      for(int j=1;j<=K;j++){
        beta(i-1,j-1)=1;
        for(int l=j;l<=(i*K);l++){
          if(alpha(i-1,l-1)!=0)
          beta(i-1,j-1)=beta(i-1,j-1)+alpha(i-1,l-1)*exp(-lfact(l-j)+(l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small)));
        }
        if(beta(i-1,j-1)<-0.01) flag=1;
      }
      
      for(int j=K+1;j<=(i*K);j++){
        for(int l=j;l<=(i*K);l++){
          if(alpha(i-1,l-1)!=0)
          beta(i-1,j-1)=beta(i-1,j-1)+alpha(i-1,l-1)*exp(-lfact(l-j)+(l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small)));
        }
        if(beta(i-1,j-1)<-0.01) flag=1;
      }
     
     if((earlystop==1)&(flag==1)) break;
    }
    
    
  }
  
  //find the workload
  arma::vec workload(N);
  double zT=0;
  
  if((earlystop==0)|(flag==0)){
  
  for(int i=1;i<=N;i++){
    workload(i-1)=sum(beta.row(i-1))/eta;
  }
  
  //find the workload at time T
  for(int j=1;j<=(N*K);j++){
    for(int l=1;l<=j;l++){
      if(beta(N-1,j-1)!=0)
      zT=zT+beta(N-1,j-1)*exp((j-l-1)*log(eta)+(j-l)*log(T-arrival(N-1)+small)-eta*(T-arrival(N-1))-lfact(j-l));
    }
    }
  }
  
  Rcpp::List out;
  out["alpha"] = alpha;
  out["flag"] = flag;
  out["beta"] = beta;
  out["workload"] = workload;
  out["zT"] = zT;
  return out;
}



// [[Rcpp::export]]
Rcpp::List wld_real(arma::vec arrival, arma::mat R, int T)
{
  
  int N=arrival.n_elem;
  int nsim=R.n_rows;
  arma::mat Z(nsim,N);

  for(int i=0;i<nsim;i++){
    Z(i,0)=R(i,0);  
    for(int j=1;j<N;j++){
      Z(i,j)=R(i,j);
      double tmp=Z(i,j-1)-(arrival(j)-arrival(j-1));
      if(tmp>0) Z(i,j)=Z(i,j)+tmp;
    }
  }
 
 //find the workload
 arma::vec workload(N);
  double zT=0;
  
 for(int j=0;j<N;j++)
   workload(j)=mean(Z.col(j));

 for(int i=0;i<nsim;i++){
   double tmp=(Z(i,N-1)-(T-arrival(N-1)))/nsim;
   if(tmp>0) zT=zT+tmp;
 }

  Rcpp::List out;
  out["Z"] = Z;
  out["workload"] = workload;
  out["zT"] = zT;
  return out;
}





// [[Rcpp::export]]
Rcpp::List wldpath_erlang(arma::vec arrival, int K, double eta, int T,arma::vec grid,int earlystop)
{
  
  int N=arrival.n_elem;
  arma::mat alpha;
  alpha.zeros(N,N*K);
  arma::mat beta;
  beta.zeros(N,N*K);
  double small=pow(10,-50);
  int flag=0;
  
  for(int j=0;j<K;j++){
    beta(0,j)=1;
  }
  
  if(N>2){
    
    for(int i=2;i<=N;i++){
      
      for(int j=1;j<=K;j++){
        
        for(int l=1;l<=((i-1)*K);l++){
          if(beta(i-2,l-1)!=0)
          alpha(i-1,j-1)=alpha(i-1,j-1)-beta(i-2,l-1)*exp(-lfact(l-1)+(K+l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small))-lfact(K-j)-log(K+l-j)-eta*(arrival(i-1)-arrival(i-2)))*pow(-1,K-j);
        }
        
      }
      
      for(int j=K+1;j<=(i*K);j++){
        
        for(int l=1;l<=K;l++){
          if(beta(i-2,j-1-K)!=0)
                      alpha(i-1,j-1)=alpha(i-1,j-1)+beta(i-2,j-1-K)*exp(-eta*(arrival(i-1)-arrival(i-2))-log(j-l)-lfact(l-1)-lfact(K-l)+lfact(j-1)-lfact(j-K-1))*pow(-1,K-l);
        }
        
      }
      
      
      for(int j=1;j<=K;j++){
        beta(i-1,j-1)=1;
        for(int l=j;l<=(i*K);l++){
          if(alpha(i-1,l-1)!=0)
          beta(i-1,j-1)=beta(i-1,j-1)+alpha(i-1,l-1)*exp(-lfact(l-j)+(l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small)));
        }
        if(beta(i-1,j-1)<-0.01) flag=1;
      }
      
      for(int j=K+1;j<=(i*K);j++){
        for(int l=j;l<=(i*K);l++){
          if(alpha(i-1,l-1)!=0)
          beta(i-1,j-1)=beta(i-1,j-1)+alpha(i-1,l-1)*exp(-lfact(l-j)+(l-j)*log(eta*(arrival(i-1)-arrival(i-2)+small)));
        }
        if(beta(i-1,j-1)<-0.01) flag=1;
      }
      
      if((earlystop==1)&(flag==1)) break;
    }
    
    
  }
  
  //find the workload
  arma::vec workload(N);
  double zT=0;
  //Find the workload path
  int ngrid=grid.n_elem;
  arma::vec workload_path(ngrid);
  
  if((earlystop==0)|(flag==0)){
    
    double gap;
    
    for(int i=1;i<=N;i++){
      workload(i-1)=sum(beta.row(i-1))/eta;
    }
    
    //find the workload at time T
    for(int j=1;j<=(N*K);j++){
      for(int l=1;l<=j;l++){
        if(beta(N-1,j-1)!=0)
        zT=zT+beta(N-1,j-1)*exp((j-l-1)*log(eta)+(j-l)*log(T-arrival(N-1)+small)-eta*(T-arrival(N-1))-lfact(j-l));
      }
    }
  
  int current=0;
  
  for(int i=0;i<ngrid;i++){
    
    if(grid(i)<arrival(0)) workload_path(i)=0; 
    if(current< (N-1)) {
      while((grid(i)>=arrival(current))&(grid(i)>= arrival(current+1))) {
        current=current +1; 
        if(current>=(N-1)) break;
      }
    }
    
    if(current< (N-1)){
      
      if((grid(i)>=arrival(current))&(grid(i)< arrival(current+1))) {
         gap = grid(i)-arrival(current);
        
        //find the workload at time T
        for(int j=1;j<=((current+1)*K);j++){
          for(int l=1;l<=j;l++){
if(beta(current,j-1)!=0)
  workload_path(i)=workload_path(i)+beta(current,j-1)*exp((j-l-1)*log(eta)+(j-l)*log(gap+small)-eta*gap-lfact(j-l));
          }
        }
        
      }
      
    }
    
    if(grid(i)>=arrival(N-1)){
      gap = grid(i)-arrival(N-1);
      //find the workload at time T
      for(int j=1;j<=(N*K);j++){
        for(int l=1;l<=j;l++){
          if(beta(N-1,j-1)!=0)
          workload_path(i)=workload_path(i)+beta(N-1,j-1)*exp((j-l-1)*log(eta)+(j-l)*log(gap+small)-eta*gap-lfact(j-l));
        }
      }
      
    } 
    
  }
  }
  Rcpp::List out;
  out["beta"] = beta;
  out["workload"] = workload;
  out["zT"] = zT;
  out["flag"] = flag;
  out["workload_path"] = workload_path;
  return out;
}



// [[Rcpp::export]]
Rcpp::List wldpath_real(arma::vec arrival, arma::mat R, int T,arma::vec grid)
{
  
  int N=arrival.n_elem;
  int nsim=R.n_rows;
  arma::mat Z(nsim,N);
  
  for(int i=0;i<nsim;i++){
    Z(i,0)=R(i,0);  
    for(int j=1;j<N;j++){
      Z(i,j)=R(i,j);
      double tmp=Z(i,j-1)-(arrival(j)-arrival(j-1));
      if(tmp>0) Z(i,j)=Z(i,j)+tmp;
    }
  }
  
  //find the workload
  arma::vec workload(N);
  double zT=0;
  
  for(int j=0;j<N;j++)
    workload(j)=mean(Z.col(j));
  
  for(int i=0;i<nsim;i++){
    double tmp=(Z(i,N-1)-(T-arrival(N-1)))/nsim;
    if(tmp>0) zT=zT+tmp;
  }
  
  //Find the workload path
  int ngrid=grid.n_elem;
  arma::vec workload_path(ngrid);
  
    
    int current=0;
    double gap;
    for(int i=0;i<ngrid;i++){
      
      if(grid(i)<arrival(0)) workload_path(i)=0; 
      if(current< (N-1)) {
        while((grid(i)>=arrival(current))&(grid(i)>= arrival(current+1))) {
          current=current +1; 
          if(current>=(N-1)) break;
        }
      }
      
      if(current< (N-1)){
        
        if((grid(i)>=arrival(current))&(grid(i)< arrival(current+1))) {
          gap = grid(i)-arrival(current);
          
          //find the workload at time T
          for(int j=0;j<nsim;j++){
            double tmp=Z(j,current)-gap;
            if(tmp>0) workload_path(i)=workload_path(i)+tmp/nsim;
          }
          
        }
        
      }
      
      if(grid(i)>=arrival(N-1)){
        gap = grid(i)-arrival(N-1);
        //find the workload at time T
        for(int j=0;j<nsim;j++){
          double tmp=Z(j,N-1)-gap;
          if(tmp>0) workload_path(i)=workload_path(i)+tmp/nsim;
        }
        
      } 
      
    }
    
  Rcpp::List out;
  out["Z"] = Z;
  out["workload"] = workload;
  out["zT"] = zT;
  out["workload_path"] = workload_path;
  return out;
}

