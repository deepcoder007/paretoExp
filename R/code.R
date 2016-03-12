# The solution of assignment from Statistical Simulation Lab
# Code for the bivariate case 
# The original parameters of the distribution
require(POT)
require(rootSolve)

alpha_0_orig = 2
alpha_1_orig = 2
alpha_2_orig = 2

mu_0_orig = 0
mu_1_orig = 1
mu_2_orig = 1

sigma_0_orig = 1
sigma_1_orig = 1
sigma_2_orig = 1

N = 100     # The number of samples to generate

# Generate the input data first 
# package for the generalized pareto distribution

Y<-data.frame(rgpd(N,loc=mu_0_orig,scale=sigma_0_orig,shape=alpha_0_orig),
              rgpd(N,loc=mu_1_orig,scale=sigma_1_orig,shape=alpha_1_orig),
              rgpd(N,loc=mu_2_orig,scale=sigma_2_orig,shape=alpha_2_orig))

colnames(Y)<-c("Y0","Y1","Y2")

X1 <- 0
X2 <- 0
for(i in 1:N)
{
  X1[i]<-min( Y$Y0[i]*sigma_1_orig+mu_1_orig , Y$Y1[i])
  X2[i]<-min( Y$Y0[i]*sigma_2_orig+mu_2_orig , Y$Y2[i] )
}

X<-data.frame(X1,X2)
#  X is the data frame of the sample
print("Sample Data Generated")


#------------------------------------------------
#   Start the EM Algo
#------------------------------------------------

# The theta_0 vector with theta_0 = (alpha_0,alpha_1,alpha_2,mu_1,mu_2,sigma_1,sigma_2)
theta<-data.frame(2,3,4,2,3,2,3)
colnames(theta)<-c("alpha_0","alpha_1","alpha_2","mu_1","mu_2","sigma_1","sigma_2")

# The w vector is w = (w_0,w_1,w_2) 
w<-data.frame(0,0,0)
colnames(w)<-c("w_0","w_1","w_2")

# the last calculated value of theta will be theta[i]
for(i in 1:10)
{
  print(c("Iteration ",i," Started"))
  theta_curr<-0   # the new vector to be initialized
  z_1 <-  (X$X1-theta[i,4])/theta[i,6]
  z_2 <-  (X$X2-theta[i,5])/theta[i,7]
  
  diff_z <- (z_1-z_2)  # The diff between z
  w_0 <- 0
  w_1 <- 0
  w_2 <- 0
  for(j in 1:N )
  {
    if( diff_z[j] ==0 ) {
      w_0 <- w_0 + 1
    } 
    if( diff_z[j] > 0 ) {
      w_1 <- w_1 + 1
    }
    if( diff_z[j] < 0 ) {
      w_2 <- w_2 + 1
    }
  }
  w<-rbind(w,c(w_0,w_1,w_2))  # The ordinalities of the z1(=,>,<)z2
  print(c("The (w_0,w_1,w_2) are :",w_0,w_1,w_2))
  
  z <- 0
  for ( j in 1:N )  # for all the X_i
  {
    
    z[j] <- max( ((X[j,1]-theta[i,4])/theta[i,6]),  
                 ((X[j,2]-theta[i,5])/theta[i,7]) )
  }
  
  theta_curr[1] <- N/( sum(log(1+z))+
                      (theta[i,2]*w_1)/(theta[i,1]*(theta[i,1]+theta[i,2]))+
                      (theta[i,3]*w_2)/(theta[i,1]*(theta[i,1]+theta[i,3])) )
  
  theta_curr[2] <- N/( sum(log(1+z_1))+
                         (theta[i,1]*w_1)/(theta[i,2]*(theta[i,1]+theta[i,2]))+
                         (w_0)/(theta[i,2]) )
  
  theta_curr[3] <- N/( sum(log(1+z_2))+
                         (theta[i,1]*w_2)/(theta[i,3]*(theta[i,1]+theta[i,2]))+
                         (w_0)/(theta[i,3]) )
  
  theta_curr[4] <- min(X$X1)
  theta_curr[5] <- min(X$X2)
  
  print(c("Current value of theta is : ",theta_curr))
  
  # solving a numerical equation to find the value of theta_curr[6] and theta_curr[7]
  
  func_sig_1 <- function(x) {
    sm<-0
    for(j in 1:N) {
      sm <- sm + 1/(x+X$X1[j]-theta_curr[4])
    }
    ans <- sm*x
    ans <- ans - (N*(theta_curr[1]+theta_curr[2]))/(theta_curr[1]+theta_curr[2]+1)
    return(ans)
  }
  
  func_sig_2 <- function(x) {
    sm<-0
    for(j in 1:N) {
      sm <- sm + 1/(x+X$X2[j]-theta_curr[5])
    }
    ans <- sm*x
    ans <- ans - (N*(theta_curr[1]+theta_curr[3]))/(theta_curr[1]+theta_curr[3]+1)
    return(ans)
  }
  
  print(c("Sigma1 : ",uniroot.all(func_sig_1,c(0,1000))))
  print(c("Sigma2 : ",uniroot.all(func_sig_2,c(0,1000))))
  break
  theta_curr[6] <- uniroot.all(func_sig_1,c(-10,10)) 
  theta_curr[7] <- uniroot.all(func_sig_2,c(-10,10))

  theta<-rbind(theta,theta_curr)
  print(theta_curr)
}



# MSE calculation
err<-(theta-theta_curr)
err2<-(err**2)
mse_alpha_0 <- mean(err2[,1])
mse_alpha_1 <- mean(err2[,2])
mse_alpha_2 <- mean(err2[,3])
mse_mu_1 <- mean(err2[,4])
mse_mu_2 <- mean(err2[,5])
mse_sigma_1 <- mean(err2[,6])
mse_sigma_2 <- mean(err2[,7])


