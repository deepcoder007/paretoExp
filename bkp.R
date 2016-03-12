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

theta_orig = c(alpha_0_orig,alpha_1_orig,alpha_2_orig,mu_1_orig,mu_2_orig,sigma_1_orig,sigma_2_orig)

theta_orig_frame <- data.frame(alpha_0_orig,alpha_1_orig,alpha_2_orig,mu_1_orig,mu_2_orig,sigma_1_orig,sigma_2_orig)
colnames(theta_orig_frame)<-c("alpha_0","alpha_1","alpha_2","mu_1","mu_2","sigma_1","sigma_2")


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

# loglikelihood of observed data set
Q_theta <- function( theta_t, X) {
  # at time (t)
  z_1_t <-  (X$X1-theta_t$mu_1)/theta_t$sigma_1
  z_2_t <-  (X$X2-theta_t$mu_2)/theta_t$sigma_2
  
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
  
  n <- dim(X)[1]
  z <- 0
  for ( j in 1:n )  # for all the X_i
  {
    
    z[j] <- max( ((X[j,1]-theta[i,4])/theta_curr[6]),  
                 ((X[j,2]-theta[i,5])/theta_curr[7]) )
  }
  
  e <- N*log(theta_t$alpha_0*theta_t$alpha_1*theta_t$alpha_2)
  for(i in 1:n) 
  {
    e <- e - theta_t$alpha_0*( log(1+z[i]) +  (theta_t$alpha_2*w_2)/(theta_t$alpha_0*(theta_t$alpha_0+theta_t$alpha_2)) + (theta_t$alpha_1*w_1)/(theta_t$alpha_0*(theta_t$alpha_0+theta_t$alpha_1)) )
    e <- e - theta_t$alpha_1*( log(1+z_1_t[i]) + (theta_t$alpha_0*w_1)/(theta_t$alpha_1*(theta_t$alpha_0 + theta_t$alpha_1)) + (w_0)/(theta_t$alpha_1) )
    e <- e - theta_t$alpha_2*( log(1+z_2_t[i]) + (theta_t$alpha_0*w_2)/(theta_t$alpha_2*(theta_t$alpha_0 + theta_t$alpha_2)) + (w_0)/(theta_t$alpha_2) )
  }
  return(e);
}

#------------------------------------------------
#   Start the EM Algo
#------------------------------------------------

# The theta_0 vector with theta_0 = (alpha_0,alpha_1,alpha_2,mu_1,mu_2,sigma_1,sigma_2)
theta<-data.frame(4,6,6,4,4,4,4)
theta<-rbind(theta,c(2,3,3,2,2,2,2))
colnames(theta)<-c("alpha_0","alpha_1","alpha_2","mu_1","mu_2","sigma_1","sigma_2")


# The w vector is w = (w_0,w_1,w_2) 
w<-data.frame(0,0,0)
colnames(w)<-c("w_0","w_1","w_2")

i<-2
err<-100000000
err_struct<-err
# the last calculated value of theta will be theta[i]
while( err > 0.0000001 )
{
  # ------------  E-STEP  --------------------------
  
  print(c("Iteration ",i-1," Started"))
  theta_curr<-0   # the new vector to be initialized
  
  #--------------- M-Step (1) ------------------------------
  
  # trivial to calculate [mu]
  theta_curr[4] <- min(X$X1)
  theta_curr[5] <- min(X$X2)
  
  # Calculate the [sigma]
  
  sm2 <- (N*(theta[i,1]+theta[i,2]))/(theta[i,1]+theta[i,2]+1)
  sm <- 0
  for(j in 1:N) {
    sm <- sm + 1/(theta[i,6]+X$X1[j]-theta[i,4])
  }
  if( abs(sm)<0.001 )
  {
    print(c("Small value 1: ",sm))
  }
  theta_curr[6] <- ( sm2 /sm )
  
  sm <- 0
  sm2 <- (N*(theta[i,1]+theta[i,3]))/(theta[i,1]+theta[i,3]+1)
  for(j in 1:N) {
    sm <- sm + 1/(theta[i,7]+X$X2[j]-theta[i,5])
  }
  
  if( abs(sm)<0.1 )
  {
    print(c("Small value 2 :",sm))
  }
  theta_curr[7] <- ( sm2 / sm )
  
  
  #-------------------------------------------
  #    E Step
  #-------------------------------------------
  
  z_1 <-  (X$X1-theta[i,4])/theta_curr[6]
  z_2 <-  (X$X2-theta[i,5])/theta_curr[7]
  
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
  
  z <- 0
  for ( j in 1:N )  # for all the X_i
  {
    
    z[j] <- max( ((X[j,1]-theta[i,4])/theta_curr[6]),  
                 ((X[j,2]-theta[i,5])/theta_curr[7]) )
  }
  
  #----------------------------------------------------
  #      M Step (2)
  #-------------------------------------------------------
  
  
  # calculate the [alpha]
  theta_curr[1] <- N/( sum(log(1+z))+
                         (theta[i,2]*w_1)/(theta[i,1]*(theta[i,1]+theta[i,2]))+
                         (theta[i,3]*w_2)/(theta[i,1]*(theta[i,1]+theta[i,3])) )
  
  theta_curr[2] <- N/( sum(log(1+z_1))+
                         (theta[i,1]*w_1)/(theta[i,2]*(theta[i,1]+theta[i,2]))+
                         (w_0)/(theta[i,2]) )
  
  theta_curr[3] <- N/( sum(log(1+z_2))+
                         (theta[i,1]*w_2)/(theta[i,3]*(theta[i,1]+theta[i,2]))+
                         (w_0)/(theta[i,3]) )
  
  
  # solving a numerical equation to find the value of theta_curr[6] and theta_curr[7]
  
  # Use fixed point iteration to calculate sigma
  
  
  theta<-rbind(theta,theta_curr)
  #  print(c("Theta is : ",theta_curr))
  err <-abs((Q_theta(theta[i,],X)-Q_theta(theta[i-1,],X))/(Q_theta(theta[i-1,],X)))
  err_struct <-c(err_struct,err)
  print(c("[Q_theta_t+1 - Q_theta_t] :",err))
  i <- i+1
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