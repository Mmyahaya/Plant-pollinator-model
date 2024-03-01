library(deSolve)
library(bipartite)


M=10# No. of plant
N=15 # No. of animal
{# Intrinsic growth rate
  rP=matrix(runif(M,0,0.1), nr=M)
  rA=matrix(runif(N,0,0.1), nr=N)
  # Density dependent
  CP=matrix(runif(M,0,0.1), nr=M)
  CA=matrix(runif(N,0,0.1), nr=N)
  # floral resource production rate
  a=matrix(rbeta(M,2,5),M,1)
  #floral resource decay rate
  w=matrix(runif(M,0,.1), nr=M)
  #conversion rate
  sigma_P=matrix(runif(M*N,0,.1), nr=M)
  sigma_A=matrix(runif(M*N,0,.1), nr=M)
  #Initial foraging effort
  bet0<-matrix(1/M,M,N)
  # Initial floral resource abundance
  Fi0<-matrix(runif(M,0,1), nr=M)
  #Initial density
  XP0<-matrix(runif(M,0,.1), nr=M)
  XA0<-matrix(runif(N,0,.1), nr=N)
  #Animal adaptation rate
  G=matrix(runif(N,1,2), nr=N)
  #Initial densities
  X0=rbind(XP0,XA0)}

# Mutualistic model function
lotka<-function(t,y,parameters){
  with(as.list(c(y,parameters)),{
    Xx<-y
    XP<-Xx[1:M]
    XA<-Xx[(M+1):(M+N)]
    Fi<-Xx[(M+N+1):(2*M+N)]
    bet<-matrix(Xx[(2*M+N+1):(M*N+2*M+N)],M,N)
    dXP<-rP-CP*XP+(((sigma_P*bet)%*%XA)*(Fi/XP))+0.1*rnorm(M,0,1) #diffusion added
    dXA<-rA-CA*XA+(t(sigma_A*bet)%*%Fi)+0.1*rnorm(N,0,1) # diffusion added
    dX<-Xx[1:(M+N)]*rbind(dXP,dXA)
    dFi<-a*XP-w*Fi-rowSums(bet*(Fi%*%t(XA)))
    dbet<-(bet%*%diag(c(G)))*sweep(sigma_A*Fi,2,colSums(bet*sigma_A*Fi),"-")
    return(list(c(dX,dFi,dbet))) 
  })
}


# Simulation of the model equation
yini=c(c(X0),c(Fi0),c(bet0))
times=seq(0,10000,0.1)
parameters=list(rP=rP,rA=rA,CP=CP,CA=CA,a=a,w=w,sigma_A=sigma_A,
                sigma_P=sigma_P, G=G)
solution<-rk4(y=yini, times=times, func=lotka, parms=parameters)


#Extraction of each state variables
{X<-as.matrix(solution[,2:(M+N+1)])
  Fi<-as.matrix(solution[,(M+N+2):(2*M+N+1)])
  ForEffMatA<-as.matrix(solution[,(2*M+N+2):(M*N+2*M+N+1)])
  XAF<-as.matrix(X[length(times),(M+1):(M+N)],nr=N)
  FiF<-as.matrix(Fi[length(times),],nr=M)
  bet<-matrix(ForEffMatA[length(times),],M,N)
  V<-bet*(FiF%*%t(XAF))
  V[V<0]<-0}

#Computation of network structural properties 
H2fun(V, H2_integer = FALSE)[1] #Specialisation
computeModules(V)@likelihood  # Modularity
nested(V,method = "WNODA") # Nestedness
networklevel(V, index = "weighted connectance")

#Plot of plant density
matplot(X[,1:M], type = "l",lwd=2,lty = "solid" , pch=1,col = 1:M, 
        main=paste(" Plant : P=",M,"A=", N), ylab = "Density", xlab = "Time")


#Plot of animal density
matplot(X[,(M+1):(M+N)], type = "l",lwd=2,lty = "solid" , pch = 1, col=1:N, 
        main=paste(" Animal : P=",M,"A=", N), ylab = "Density", xlab = "Time")


#Plot of animal foraging effort
matplot(ForEffMatA, type = "l", lwd=2,lty ="solid" ,pch = 1, col = rep(1:N, each=M), 
        main=paste(" Animal foraging effort : P=",M,"A=", N) ,
        ylab = "Foraging Effort ", xlab = "Time")

#Plot of floral resource abundance
matplot(Fi, type = "l",lty = "solid" ,lwd=2, pch = 1, col=1:M,
        main=paste("Floral resource: P=",M, "A=",N), 
        ylab = "Floral resource", xlab = "Time")



