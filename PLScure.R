list.of.packages <- c("MASS","plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(MASS)
library(plyr)

PLScureSIM<-function(seed=NA,n,alpha,beta,gamma,scen=1){
  if(is.numeric(seed)){set.seed(seed)}
  if(!scen%in%c(1,2,3)){print("You have to input a value (1, 2, or 3) into argument `scen`");break}
  
  alpha_norm<-sqrt(sum(alpha^2))
  gamma_norm<-sqrt(sum(gamma^2))
  balpha<-alpha/alpha_norm
  bgamma<-gamma/gamma_norm
  
  data<-NULL
  
  for(i in c(1:n)){
    X1<-2*rbinom(n = 1,size = 1,prob = 0.5)-1
    X2<-rnorm(n = 1,mean = 0,sd = 1)
    X3<-rnorm(n = 1,mean = 0,sd = 1)
    X4<-rnorm(n = 1,mean = 0,sd = 1)
    XX<-c(X1,X2,X3,X4)
    ZZ<-c(X3,X4)
    
    SI1<-XX%*%balpha
    SI2<-ZZ%*%bgamma
    
    #incidence link function
    if(scen==1){eta<-exp(SI1)/(1+exp(SI1))}
    if(scen==2){eta<-0.4+0.6*pnorm(2*SI1-1)}
    if(scen==3){eta<-0.5*pnorm(2*SI1+3)+0.5*pnorm(2*SI1-3)}
    
    uncured<-rbinom(n = 1,size = 1,prob = eta)
    
    Ci<-rexp(n = 1,rate = 0.2)
    if(uncured==0){
      cure   <-1
      Ti     <-Inf
    }else if(uncured==1){
      cure   <-0
      St     <-runif(n = 1,min = 0,max = 1)
      if(scen==1){Ti<-((-log(St))/exp(c(X1,X2)%*%beta+SI2)/0.8)^(1/1.2)}
      if(scen==2){Ti<-((-log(St))/exp(c(X1,X2)%*%beta+log(1+SI2^2))/0.8)^(1/1.2)}
      if(scen==3){Ti<-((-log(St))/exp(c(X1,X2)%*%beta+sin(3/2*SI2))/0.8)^(1/1.2)}
    }
    
    Yi<-min(Ti,Ci)
    cen<-ifelse(Ti<=Ci,1,0)
    
    data<-rbind(data,c(i,Yi,cen,XX,cure))
  }
  data<-data.frame(data)
  names(data)<-c("id","Yi","cen","X1","X2","X3","X4","cure")
  return(data)
}

Logit <-function(x){ifelse(exp(as.matrix(x))==Inf,1,as.numeric(exp(as.matrix(x))/(1+exp(as.matrix(x)))))}

BP_G<-function(para,m10,x.scaled){
  psi   <-c(para[1],para[1]+cumsum(exp(para[2:(m10+1)])))
  ii    <-c(0:m10)
  cii   <-choose(m10,ii)
  At    <-t(sapply(c(1:length(x.scaled)),function(o){cii*x.scaled[o]^ii*(1-x.scaled[o])^(m10-ii)}))
  return(as.numeric(At%*%psi))
}

BP_H<-function(para,m20,z.scaled){
  phi   <-para
  jj    <-c(0:m20)
  cjj   <-choose(m20,jj)
  Bt    <-t(sapply(c(1:length(z.scaled)),function(o){cjj*z.scaled[o]^jj*(1-z.scaled[o])^(m20-jj)}))
  return(as.numeric(Bt%*%phi))
}

obs.like<-function(psi, phi, cumhaz, alpha, beta, gamma, X, W, Z, p, q, r, m10, m20, Xmin, Xmax, Zmin, Zmax, cen, Yi, maxT){
  betaW     <-as.numeric(W%*%beta) 
  SI1       <-as.numeric(X%*%alpha)
  SI1.scaled<-(SI1-Xmin)/(Xmax-Xmin)
  
  if(r==1){SI2.scaled<-(as.numeric(Z)-Zmin)/(Zmax-Zmin)}else{
    SI2       <-as.numeric(Z%*%gamma)
    SI2.scaled<-(SI2-Zmin)/(Zmax-Zmin)
  }
  
  GSI <-Logit(BP_G(para = psi,m10 = m10,x.scaled = SI1.scaled))
  HSI <-BP_H(para = phi,m20 = m20,z.scaled = SI2.scaled)
  
  lambda<-c(cumhaz[1],diff(cumhaz))
  
  return(ifelse(cen==1,(log(GSI)+log(lambda)+betaW+HSI-cumhaz*exp(betaW+HSI)), #log(1-GSI+GSI*exp(-cumhaz*exp(betaW+HSI)))))
                ifelse(Yi>maxT,log(1-GSI),log(1-GSI+GSI*exp(-cumhaz*exp(betaW+HSI))))))
}

ini.psi<-function(par, SI1.scaled, m10){sum((Logit(BP_G(para = par[1:(m10+1)],m10 = m10,x.scaled = SI1.scaled))-seq(0,1,length.out=10))^2)}

Q1_psi<-function(par, Eu, m10, SI1.scaled){
  GI<-BP_G(para = par,m10 = m10,x.scaled = SI1.scaled)
  -sum(Eu*GI-log(1+exp(GI))) #check to be correct
}

Q3_beta_phi<-function(par, q, m20, W, n, Eu, Nd, cumhaz.d, SI2.scaled, Rset){
  beta <-par[1:q]; phi <-par[(q+1):(q+1+m20)]
  betaW<-apply(matrix(rep(beta,n) ,nrow=n,byrow=TRUE)*W,MARGIN = 1,sum)
  HSI  <-BP_H(para = phi,m20 = m20,z.scaled = SI2.scaled)
  MM<-0
  for(k in Nd){MM<-MM+betaW[k]+HSI[k]-log(sum((Eu*exp(betaW+HSI))[Rset[[k]]]))}
  return(MM)
}

EM<-function(psi.d, phi.d, cumhaz.d, alpha.d, beta.d, gamma.d, p, q, r, X, W, Z, m10, m20, Xmin, Xmax, Zmin, Zmax, 
             n, Rset, cen, Yi, maxT, NPupdate, tol){
  conv<-999
  c1.d<-0
  c2.d<-0
  times<-0
  repeat{
    betaW     <-as.numeric(W%*%beta.d) 
    SI1       <-as.numeric(X%*%alpha.d)
    SI1.scaled<-(SI1-Xmin)/(Xmax-Xmin)
    
    if(r==1){SI2.scaled<-(as.numeric(Z)-Zmin)/(Zmax-Zmin)}else{
      SI2       <-as.numeric(Z%*%gamma.d)
      SI2.scaled<-(SI2-Zmin)/(Zmax-Zmin)
    }
    
    #E-step
    GSI <-Logit(BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled))
    HSI <-BP_H(para = phi.d,m20 = m20,z.scaled = SI2.scaled)
    Su  <-ifelse(Yi>maxT,0,exp(-cumhaz.d*exp(betaW+HSI)))
    Eu  <-cen+(1-cen)*ifelse(GSI==1,1,(GSI*Su)/(1-GSI+GSI*Su))
    
    cumhaz.dummy<-cumhaz.d; phi.dummy<-phi.d; psi.dummy<-psi.d; alpha.dummy<-alpha.d; beta.dummy <-beta.d; gamma.dummy<-gamma.d
    
    XER<-try(optim(par = psi.d, fn = Q1_psi,gr = NULL, method="BFGS",m10=m10,SI1.scaled=SI1.scaled,Eu=Eu, 
                   hessian=FALSE, control=list(reltol=10^{-8},maxit=1000)),silent = TRUE)
    if(inherits(XER, "try-error")){NULL}else{psi.d<-XER$par}
    
    #E-step
    GSI <-Logit(BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled))
    Eu  <-cen+(1-cen)*ifelse(GSI==1,1,(GSI*Su)/(1-GSI+GSI*Su))
    
    if(NPupdate==FALSE){
      tb<-0
      repeat{
        alpha.dummy1<-alpha.d
        c1.dummy1   <-c1.d
        SI1         <-as.numeric(X%*%alpha.d)
        SI1.scaled  <-(SI1-Xmin)/(Xmax-Xmin)
        
        psi.star<-c(psi.d[1],psi.d[1]+cumsum(exp(psi.d[2:(m10+1)])))
        GSI<-Logit(BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled))
        GI <-BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled)
        
        GSI.prime<-sapply(c(1:n), function(i){m10*sum(psi.star*sapply(c(0:m10),function(j){
          ifelse(choose(m10-1,j-1)==0,0,choose(m10-1,j-1)*(SI1.scaled[i])^(j-1)*(1-SI1.scaled[i])^(m10-j))-
            ifelse(choose(m10-1,j)==0,0,choose(m10-1,j)*(SI1.scaled[i])^j*(1-SI1.scaled[i])^(m10-1-j))}))}) #checked to be correct
        
        GSI.double.prime<-sapply(c(1:n), function(i){m10*(m10-1)*sum(psi.star*sapply(c(0:m10),function(j){
          ifelse(choose(m10-2,j-2)==0,0,choose(m10-2,j-2)*(SI1.scaled[i])^(j-2)*(1-SI1.scaled[i])^(m10-j))-
            2*ifelse(choose(m10-2,j-1)==0,0,choose(m10-2,j-1)*(SI1.scaled[i])^(j-1)*(1-SI1.scaled[i])^(m10-2-(j-1)))+
            ifelse(choose(m10-2,j)==0,0,choose(m10-2,j)*(SI1.scaled[i])^j*(1-SI1.scaled[i])^(m10-2-j))}))}) #checked to be correct
        
        dQ1a<-as.numeric(sapply(c(1:p), function(o){(Xmax-Xmin)^{-1}*sum(GSI.prime*X[,o]*(Eu-GSI))}))
        
        dQ1a2<-matrix(0,nrow=p+1, ncol=p+1)
        for(i in c(1:p)){
          for(j in c(1:i)){
            dQ1a2[i,j]<-dQ1a2[j,i]<-(Xmax-Xmin)^{-2}*sum(X[,i]*X[,j]*(-GSI.prime^2*exp(GI)*(1-GSI)^2+(Eu-GSI)*GSI.double.prime))
          }}
        
        dQ1a2<-dQ1a2+cbind(rbind(diag(x = 2*c1.d,nrow = p),0),0)
        dQ1a2[(p+1),1:p]<-2*alpha.d; dQ1a2[1:p,(p+1)]<-2*alpha.d
        
        if((sum(is.finite(dQ1a2))!=length(dQ1a2))|(sum(is.finite(dQ1a))!=length(dQ1a))){break}
        alpha_c1_inc<-as.numeric(ginv(dQ1a2)%*%(c((dQ1a+2*c1.d*alpha.d), (alpha.d%*%alpha.d-1))))
        
        K<-0
        repeat{
          alpha_c1_test<-as.numeric(c(alpha.d, c1.d)-2^{-K}*alpha_c1_inc)
          alpha_test<-alpha_c1_test[1:p]
          c1_test   <-alpha_c1_test[p+1]
          
          if(max(abs(as.numeric(X%*%alpha_test)))<=Xmax|K>100){
            alpha.d <-alpha_test
            c1.d    <-c1_test
            break
          }
          K<-K+1
        }
        tb<-tb+1
        if((max(abs(c(alpha.d-alpha.dummy1, c1.d-c1.dummy1)))<10^{-4})|tb==100){break}
      }
      if(tb==100){conv<-1}
    }
    
    Bk<-matrix(0,nrow = n,ncol = (m20+1))
    for(i in c(1:n)){
      for(j in c(1:(m20+1))){
        Bk[i,j]<-choose(m20,(j-1))*(SI2.scaled[i])^(j-1)*(1-SI2.scaled[i])^(m20-j+1)
      }}
    
    beta.dummy1<-beta.d
    phi.dummy1 <-phi.d
    
    Nd<-which(cen==1)
    
    dLp_beta<-sapply(c(1:q),function(k){sum(sapply(Nd, function(o1){
      W[o1,k]-(sum((Eu*exp(betaW+HSI))[Rset[[o1]]]))^{-1}*(sum((Eu*W[,k]*exp(betaW+HSI))[Rset[[o1]]]))
    }) )})
    
    dLp_phi<-sapply(c(1:(m20+1)),function(k){sum(sapply(Nd, function(o1){
      Bk[o1,k]-(sum((Eu*exp(betaW+HSI))[Rset[[o1]]]))^{-1}*(sum((Eu*Bk[,k]*exp(betaW+HSI))[Rset[[o1]]]))
    }) )})
    
    dLp2<-matrix(0,nrow=q+(m20+1), ncol=q+(m20+1))
    
    vk<-(Eu*exp(betaW+HSI))
    for(i in c(1:(q+m20+1))){
      for(j in c(1:i)){
        if(i<=q){Ai<-W[,i]}else{Ai<-Bk[,(i-q)]}
        if(j<=q){Aj<-W[,j]}else{Aj<-Bk[,(j-q)]}
        MM<-0
        for(k in Nd){
          Vk  <-sum(vk[Rset[[k]]])
          Vki <-sum((Ai*vk)[Rset[[k]]])
          Vkj <-sum((Aj*vk)[Rset[[k]]])
          Vkij<-sum((Ai*Aj*vk)[Rset[[k]]])
          MM<-MM+Vk^{-1}*(Vkij-Vki*Vk^{-1}*Vkj)
        }
        dLp2[j,i]<-dLp2[i,j]<- -MM
      }}
    
    current.Q3<-Q3_beta_phi(par = c(beta.d,phi.d), q=q, W=W, m20=m20, n=n, Eu=Eu, Nd=Nd, 
                            cumhaz.d=cumhaz.d, SI2.scaled=SI2.scaled, Rset=Rset)
    
    if(NPupdate==FALSE){
      if((sum(is.finite(dLp2))!=length(dLp2))|(sum(is.finite(c(dLp_beta,dLp_phi)))!=length(c(dLp_beta,dLp_phi)))){conv<-1}else{
        beta_phi_inc<-as.numeric(ginv(dLp2)%*%(c(dLp_beta,dLp_phi)))
        K<-0
        repeat{
          beta_phi_test<-as.numeric(c(beta.d, phi.d)-2^{-K}*beta_phi_inc)
          beta_test<-beta_phi_test[1:q]
          phi_test <-beta_phi_test[(q+1):(q+1+m20)]-BP_H(para=beta_phi_test[(q+1):(q+1+m20)],m20 = m20,z.scaled = 0.5)
          
          new.Q3<-Q3_beta_phi(par = c(beta_test,phi_test), q=q, W=W, m20=m20, n=n, Eu=Eu, Nd=Nd,  
                              cumhaz.d=cumhaz.d, SI2.scaled=SI2.scaled, Rset=Rset)
          
          if((is.finite(new.Q3)&(new.Q3>=current.Q3))|K>100){
            beta.d <-beta_test
            phi.d  <-phi_test
            break
          }
          K<-K+1
        }
      }}else{
        phi_inc<-as.numeric(ginv(dLp2[(q+1):(q+m20+1),(q+1):(q+m20+1)])%*%(c(dLp_phi)))
        K<-0
        repeat{
          phi_test0<-as.numeric(phi.d-2^{-K}*phi_inc)
          phi_test <-phi_test0-BP_H(para=phi_test0,m20 = m20,z.scaled = 0.5)
          
          new.Q3<-Q3_beta_phi(par = c(beta.d,phi_test), q=q, W=W, m20=m20, n=n, Eu=Eu, Nd=Nd, 
                              cumhaz.d=cumhaz.d, SI2.scaled=SI2.scaled, Rset=Rset)
          
          if((is.finite(new.Q3)&(new.Q3>=current.Q3))|K>100){
            phi.d  <-phi_test
            break
          }
          K<-K+1
        }
      }
    
    #E-step
    GSI <-Logit(BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled))
    HSI <-BP_H(para = phi.d,m20 = m20,z.scaled = SI2.scaled)
    Su  <-ifelse(Yi>maxT,0,exp(-cumhaz.d*exp(betaW+HSI)))
    Eu  <-cen+(1-cen)*ifelse(GSI==1,1,(GSI*Su)/(1-GSI+GSI*Su))
    
    #update baseline cumulative hazard
    cum.den <-as.numeric(lapply(X = Rset, FUN = function(o){sum((Eu*exp(betaW+HSI))[o])}))
    cumhaz.d<-cumsum(ifelse(cen==0,0,1/cum.den))
    
    if(r!=1){
      if(NPupdate==FALSE){
        td<-0
        repeat{
          gamma.dummy1<-gamma.d
          c2.dummy1   <-c2.d
          SI2         <-as.numeric(Z%*%gamma.d)
          SI2.scaled  <-(SI2-Zmin)/(Zmax-Zmin)
          
          HSI<-BP_H(para = phi.d,m20 = m20,z.scaled = SI2.scaled)
          
          HSI.prime<-sapply(c(1:n), function(i){m20*sum(phi.d*sapply(c(0:m20),function(j){
            ifelse(choose(m20-1,j-1)==0,0,choose(m20-1,j-1)*(SI2.scaled[i])^(j-1)*(1-(SI2.scaled[i]))^(m20-j))-
              ifelse(choose(m20-1,j)==0,0,choose(m20-1,j)*(SI2.scaled[i])^j*(1-(SI2.scaled[i]))^(m20-1-j))}))})
          
          HSI.double.prime<-sapply(c(1:n), function(i){m20*(m20-1)*sum(phi.d*sapply(c(0:m20),function(j){
            ifelse(choose(m20-2,j-2)==0,0,choose(m20-2,j-2)*((SI2.scaled[i]))^(j-2)*(1-(SI2.scaled[i]))^(m20-j))-
              2*ifelse(choose(m20-2,j-1)==0,0,choose(m20-2,j-1)*((SI2.scaled[i]))^(j-1)*(1-(SI2.scaled[i]))^(m20-1-j))+
              ifelse(choose(m20-2,j)==0,0,choose(m20-2,j)*((SI2.scaled[i]))^j*(1-(SI2.scaled[i]))^(m20-2-j))}))})
          
          dLp_gamma<-sapply(c(1:r), function(k){sum((Zmax-Zmin)^{-1}*HSI.prime*Z[,k]*(cen-Eu*cumhaz.d*exp(betaW+HSI)))})
          
          dLp2_gamma<-matrix(0,nrow=r+1, ncol=r+1)
          for(i in c(1:r)){
            for(j in c(1:i)){
              dLp2_gamma[i,j]<-dLp2_gamma[j,i]<-(Zmax-Zmin)^{-2}*
                sum(Z[,i]*Z[,j]*(cen*HSI.double.prime-Eu*cumhaz.d*exp(betaW+HSI)*(HSI.double.prime+HSI.prime^2)))
            }}
          
          dLp2_gamma<-dLp2_gamma+cbind(rbind(diag(x = 2*c2.d,nrow = r),0),0)
          dLp2_gamma[(r+1),1:r]<-2*gamma.d; dLp2_gamma[1:r,(r+1)]<-2*gamma.d
          
          if((sum(is.finite(dLp2_gamma))!=length(dLp2_gamma))|(sum(is.finite(dLp_gamma))!=length(dLp_gamma))){break}
          
          gamma_c2_inc<-as.numeric(ginv(dLp2_gamma)%*%(c((dLp_gamma+2*c2.d*gamma.d), (gamma.d%*%gamma.d-1))))
          
          K<-0
          repeat{
            gamma_c2_test<-as.numeric(c(gamma.d, c2.d)-2^{-K}*gamma_c2_inc)
            gamma_test<-gamma_c2_test[1:r]
            c2_test   <-gamma_c2_test[r+1]
            
            if(max(abs(as.numeric(Z%*%gamma_test)))<=Zmax|K>100){
              gamma.d <-gamma_test
              c2.d    <-c2_test
              break
            }
            K<-K+1
          }
          
          td<-td+1
          if((max(abs(c(gamma.d-gamma.dummy1, c2.d-c2.dummy1)))<10^{-4})|td==100){break}
        }
        
        if(gamma.d[1]<0){
          gamma.d <- (-gamma.d)
          phi.d <- rev(phi.d)
        }
        if(td==100){conv<-1}
      }
    }
    
    dist<-max(abs(c(psi.d-psi.dummy, phi.d-phi.dummy, alpha.d-alpha.dummy, beta.d-beta.dummy, gamma.d-gamma.dummy)))
    
    times<-times+1
    like<-obs.like(psi = psi.d,phi = phi.d,cumhaz = cumhaz.d,alpha = alpha.d,beta = beta.d,gamma = gamma.d,
                   X = X,W = W,Z = Z,m10 = m10,m20 = m20, p=p,q=q,r=r,
                   Xmin = Xmin, Xmax = Xmax, Zmin= Zmin, Zmax = Zmax, cen = cen, Yi=Yi, maxT=maxT)
    
    if(conv==1){break}
    if(dist<tol){conv<-0;break}
    if((NPupdate==FALSE)&(times>100)){conv<-0;break}
    if((NPupdate==TRUE)&(times>100)){conv<-0;break}
  }
  inc.prob<-cbind(SI1,GSI)
  H.value <-cbind(SI2,HSI)
  return(list(like=like,times=times,cumhaz.h=cumhaz.d, psi.h=psi.d, phi.h=phi.d,
              alpha.h=alpha.d, beta.h=beta.d, gamma.h=gamma.d,conv=conv,Pi=inc.prob,HgamZ=H.value))
}

PLScureEST<-function(X,W,Z,Yi,cen,K1=3,K2=5,M2=1:5,tolerance=10^{-4},attempt=10,SE_est=TRUE,TRACE=FALSE){
  M1<-K1*floor(length(cen)^{1/4})
  
  X<-data.frame(scale(X))
  W<-data.frame(scale(W))
  Z<-data.frame(scale(Z))
  
  ascd<-order(Yi)
  X<-as.matrix(as.matrix(X)[ascd,])
  W<-as.matrix(as.matrix(W)[ascd,])
  Z<-as.matrix(as.matrix(Z)[ascd,])
  Yi<-Yi[ascd]
  cen<-cen[ascd]
  
  n <-length(cen)
  
  p<-ncol(X)
  q<-ncol(W)
  r<-ncol(Z)
  
  if(p==1){Xmax<- max(X); Xmin<- min(X)}else{Xmax<- quantile(sqrt(apply(X^2,MARGIN = 1,sum)),probs = 1); Xmin<- -Xmax} 
  if(r==1){Zmax<- max(Z); Zmin<- min(Z)}else{Zmax<- quantile(sqrt(apply(Z^2,MARGIN = 1,sum)),probs = 1); Zmin<- -Zmax}
  
  maxT <-max(cen*Yi)
  
  Rset <- vector("list", length=n)
  for(i in 1:length(Rset)){Rset[[i]] <- c(1:length(cen))[Yi>=Yi[i]]}
  
  ini<-0
  store.full<-c()
  Pi.store<-HgamZ.store<-vector(mode = 'list', length = length(M1)*length(M2))
  for(m10 in M1){
    for(m20 in M2){
      l<-S<-0; discard.old<-discard<-0
      repeat{
        if((TRACE==TRUE)&(discard.old==discard)){print(paste0("m10: ", m10,", m20: ", m20, ", Attempt: ", S+1))}
        discard.old<-discard
        set.seed(m10*100+m20*10+l)
        N_a<-rnorm(n = p,0,1)
        N_g<-rnorm(n = r,0,1); N_g[1]<-abs(N_g[1])
        true.alpha<-N_a/sqrt(sum(N_a^2))
        true.beta <-rep(0,q)
        true.gamma<-N_g/sqrt(sum(N_g^2))
        
        SI.ini<-seq(0,1,length.out=12)[2:11]
        psi0 <-optim(par = rep(0,m10+1),fn = ini.psi, m10=m10, method = "BFGS",SI1.scaled=SI.ini,control = list(reltol=10^{-4},maxit=1000))$par
        
        EMstore<-EM(psi.d = psi0, phi.d = rep(0,m20+1), cumhaz.d = cumsum(ifelse(cen==1,1/sum(cen),0)), alpha.d = true.alpha,
                    beta.d = true.beta, gamma.d = true.gamma, p = p, q = q, r = r, X = X, W = W, Z = Z, m10 = m10, m20 = m20, n=n,
                    Rset=Rset, Xmin = Xmin, Xmax = Xmax, Zmin = Zmin, Zmax = Zmax, cen = cen,Yi = Yi, maxT = maxT, NPupdate=FALSE, tol=tolerance)
        
        aa<-c(m10,m20,S,sum(EMstore$like),
              EMstore$times,EMstore$cumhaz.h,EMstore$psi.h,rep(NA,max(M1)-m10),EMstore$phi.h,rep(NA,max(M2)-m20),
              EMstore$alpha.h,EMstore$beta.h,EMstore$gamma.h)
        
        if(EMstore$conv==1){discard<-discard+1}
        if(EMstore$conv==0){
           ini<-ini+1
          Pi.store[[ini]]<-EMstore$Pi
          HgamZ.store[[ini]]<-EMstore$HgamZ
          store.full<-rbind(store.full,c(aa,discard))
          S<-S+1
        }
        l<-l+1
        if(S==attempt){break}
      }
    }
  }
  
  store.case<-as.data.frame(store.full)
  names(store.case)<-c("m10","m20","l","like","iteration",paste0("Lam",c(1:n)),
                       paste0("psi",c(1:(max(M1)+1))), paste0("phi",c(1:(max(M2)+1))),
                       paste0("alpha",c(1:p)),paste0("beta",c(1:q)),paste0("gamma",c(1:r)),"discard")
  store.case$AICs<- -2*store.case$like+2*(store.case$m10+store.case$m20+p+q+r)
  
  selected<-which(store.case$AICs==min(store.case$AICs))[1]
  store.star<-as.matrix(store.case[selected,])
  m10       <-store.star[1]
  m20       <-store.star[2]
  like.star <-store.star[4]
  times.star<-store.star[5]
  cumhaz.d  <-store.star[6:(n+5)]
  psi.d     <-store.star[(n+5+1):(n+5+m10+1)]
  phi.d     <-store.star[(n+5+max(M1)+2):(n+5+max(M1)+m20+2)]
  alpha.d   <-store.star[(n+5+max(M1)+max(M2)+3):(n+5+max(M1)+max(M2)+p+2)]
  beta.d    <-store.star[(n+5+max(M1)+max(M2)+p+3):(n+5+max(M1)+max(M2)+p+q+2)]
  gamma.d   <-store.star[(n+5+max(M1)+max(M2)+p+q+3):(n+5+max(M1)+max(M2)+p+q+r+2)]
  Pi.d      <-as.data.frame(Pi.store[[selected]])
  HgamZ.d   <-as.data.frame(HgamZ.store[[selected]])
  names(Pi.d)   <-c("alphaX","P(U=1)")
  names(HgamZ.d)<-c("gammaZ","H(gammaZ)")
  
  point.Euclidean<-c(alpha.d,beta.d,gamma.d)
  
  grid<-1000
  XX  <-seq(Xmin,Xmax, length.out = grid)
  ZZ  <-seq(Zmin,Zmax, length.out = grid)
  SI1.scaled<-(XX-Xmin)/(Xmax-Xmin)
  SI2.scaled<-(ZZ-Zmin)/(Zmax-Zmin)
  
  G.matrix<-numeric(grid)
  for(j in c(1:grid)){G.matrix[j]<-Logit(BP_G(para = psi.d,m10 = m10,x.scaled = SI1.scaled[j]))}
  plot(x = XX, y = G.matrix,xlim=c(Xmin,Xmax),ylim=c(-0.05,1),main="Uncured probabiltiy", type="l", ylab="", xlab="")
  title(ylab = expression(paste(hat(pi),"(",alpha^"T","X)")), line = 2.5, cex.lab = 1)
  title(xlab = expression(paste(alpha^"T","X")), line = 2.5, cex.lab = 1)
  points(x = X%*%alpha.d, y = rep(-0.05,length(X%*%alpha.d)),pch="|",cex=0.6)
  
  H.matrix<-numeric(grid)
  for(j in c(1:grid)){H.matrix[j]<-BP_H(para = phi.d, m20 = m20, z.scaled = SI2.scaled[j])}
  plot(x = ZZ, y = H.matrix,type="l",xlim=c(Zmin,Zmax),ylim=c(min(H.matrix)-0.05,max(H.matrix)), ylab="", xlab="",main="H-function")
  title(ylab = expression(paste(hat("H"),"(",gamma^"T","Z)")), line = 2.5, cex.lab = 1)
  title(xlab = expression(paste(gamma^"T","Z")), line = 2.5, cex.lab = 1)
  points(x = Z%*%gamma.d, y = rep(min(H.matrix)-0.05,length(Z%*%gamma.d)),pch="|",cex=0.6)
  
  if(SE_est==FALSE){
    return(list(alpha=as.numeric(alpha.d), beta=as.numeric(beta.d), gamma=as.numeric(gamma.d), 
                psi=c(psi.d[1],psi.d[1]+cumsum(exp(psi.d[2:(m10+1)]))), phi=phi.d,
                m10.star=m10, m20.star=m20, Xmax=Xmax, Zmax=Zmax,
                likelihood=like.star, AIC=-2*like.star+2*(length(point.Euclidean)+m10+m20+2),
                Pi=Pi.d, HgamZ=HgamZ.d, AIC.by.m10.m20=ddply(store.case,m10~m20,summarize,AIC=min(AICs)) ))
    
  }else if(SE_est==TRUE){
    grad.old<-replicate(obs.like(psi = psi.d, phi = phi.d, gamma = gamma.d,cumhaz = cumhaz.d, 
                                 alpha = alpha.d, beta = beta.d,X=X, W=W, Z=Z, m10=m10, m20=m20,p=p,q=q,r=r,
                                 Xmin=Xmin, Xmax=Xmax, Zmin=Zmin, Zmax=Zmax, cen=cen,Yi=Yi, maxT=maxT),n=(p-1)+q+(r-1))
    
    aI         <-which.max(abs(alpha.d))
    grad.new   <-matrix(0,nrow = n,ncol = (p-1)+q+(r-1))
    K2n        <-K2*n^{-1/2}
    hn         <-c(pmin(K2n, abs(alpha.d[-aI])), rep(K2n,q), pmin(K2n, abs(gamma.d[-1])))
    eu.sign    <-c(-sign(alpha.d[-aI]),rep(1,q),-sign(gamma.d[-1]))
    
    for(I in c(1:(p-1+q+r-1))){
      Eu.par    <-c(alpha.d[-aI],beta.d, gamma.d[-1])
      Eu.par[I] <-Eu.par[I]+hn[I]*eu.sign[I]
      perb.alpha<-append(Eu.par[1:(p-1)], sign(alpha.d[aI])*sqrt(1-sum(Eu.par[1:(p-1)]^2)),after = aI-1)
      
      perb.beta <-Eu.par[p:(p-1+q)]
      if(r==1){perb.gamma<-1}else{perb.gamma<-c(sqrt(1-sum(Eu.par[(p+q):(p+q+r-2)]^2)),Eu.par[(p+q):(p+q+r-2)])}
      
      EMstore<-EM(psi.d = psi.d, phi.d = phi.d, cumhaz.d = cumhaz.d, alpha.d = perb.alpha, 
                  beta.d = perb.beta, gamma.d = perb.gamma, p=p, q=q, r=r, X = X,W = W,Z = Z,m10 = m10,m20 = m20, n=n,
                  Rset=Rset, Xmin = Xmin,Xmax = Xmax,Zmin = Zmin,Zmax = Zmax,cen = cen,Yi = Yi,maxT = maxT, NPupdate=TRUE, tol=tolerance)
      
      grad.new[,I]<-EMstore$like
    }
    
    dli <-(grad.new-grad.old)/matrix(rep(hn*eu.sign,n),byrow=TRUE,nrow=n)
    
    dli2<-matrix(0,(p-1)+q+(r-1),(p-1)+q+(r-1))
    for(j in c(1:n)){dli2<-dli2+as.vector(dli[j,])%o%as.vector(dli[j,])}
    
    var.est    <-ginv(dli2)
    se.est     <-sqrt(diag(var.est))
    
    if(p==1){
      SE.vector<-c(NA,se.est[1:q])
    }else{
      a_ratio<-as.matrix(alpha.d[-aI]/alpha.d[aI])
      sd.alpha.p<-sqrt(t(a_ratio)%*%var.est[c(1:(p-1)),c(1:(p-1))]%*%a_ratio)
      SE.vector<-c(append(se.est[1:(p-1)],sd.alpha.p,aI-1),se.est[p:(p-1+q)])
    }
    if(r==1){
      SE.vector<-c(SE.vector,NA)
    }else{
      sd.gamma.r<-sqrt(t(as.vector(gamma.d[-1]/gamma.d[1]))%*%var.est[c((p+q):(p+q+r-2)),c((p+q):(p+q+r-2))]%*%t(t(as.vector(gamma.d[-1]/gamma.d[1]))))
      SE.vector<-c(SE.vector,sd.gamma.r, se.est[c((p+q):(p+q+r-2))])
    }
    CI.left   <-point.Euclidean-1.96*SE.vector
    CI.right  <-point.Euclidean+1.96*SE.vector
    
    results<-round(data.frame(estimate=point.Euclidean, SE=SE.vector, CI.left=CI.left, CI.right=CI.right),3)
    results<-cbind(data.frame(coef=c(paste0("alpha",c(1:p)),paste0("beta",c(1:q)),paste0("gamma",c(1:r)))),results)
    
    return(list(summary=results, alpha=as.numeric(alpha.d), alpha.se=as.numeric(SE.vector[1:p]), 
                beta=as.numeric(beta.d), beta.se=as.numeric(SE.vector[(p+1):(p+q)]),
                gamma=as.numeric(gamma.d), gamma.se=as.numeric(SE.vector[(p+q+1):(p+q+r)]),
                psi=c(psi.d[1],psi.d[1]+cumsum(exp(psi.d[2:(m10+1)]))), phi=phi.d,
                m10.star=m10, m20.star=m20, Xmax=Xmax, Zmax=Zmax, 
                likelihood=like.star, AIC=-2*like.star+2*(length(point.Euclidean)+m10+m20+2),
                Pi=Pi.d, HgamZ=HgamZ.d, AIC.by.m10.m20=ddply(store.case,m10~m20,summarize,AIC=min(AICs)) ))
  }
}
