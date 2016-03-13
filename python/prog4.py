# Python program for c

# This code is having stopping criteria as Q(theta)

# Param vector will have 7 elements as given in the seq:
# param = [alpha0,alpha1,alpha2,mu1,mu2,sigma1,sigma2]
#

# NOTE: config is of the form [ iterCnt, N , alpha0 , alpha1, alpha2 , mu1 , mu2 , sigma1 , sigma2 ]

import math
import random
import csv
import sys


# get samples from the pareto distribution
def getPa(N,alpha,mu,sigma):
    out = []
    for i in range(N):
        u = math.exp(-math.log(random.random())/alpha)
        out.append( mu + sigma*(u-1) )
    return out

# Generate N bivariate datapoints
def generate(N,param):   
    dist1 = getPa(N,param[0],0.0,1.0)
    dist2 = getPa(N,param[1],param[3],param[5])
    dist3 = getPa(N,param[2],param[4],param[6])
    out = []
    for i in range(N):
        x = []
        y = param[5]*dist1[i]+param[3]
        x.append( min(y,dist2[i]) )
        y = param[6]*dist1[i]+param[4]
        x.append( min(y,dist3[i]) )
        out.append(tuple(x))
    return out

# returns a param of (alpha0_(k+1),alpha1_(k+1),alpha2_(k+1),mu1_k,mu2_k,sigma1_k,sigma2_k)
def EM1Pas(N,data,param):
    sum0,sum1,sum2 = 0.0,0.0,0.0
    for i in range(N):
        z1 = ( data[i][0] - param[3] )/param[5]
        z2 = ( data[i][1] - param[4] )/param[6]
        sum1 += math.log( 1+z1 )
        sum2 += math.log( 1+z2 )
        sum0 += math.log( 1+max(z1,z2) )

    neg = N*param[0]/(param[0]+param[1]+param[2])
    n1  = N*param[2]/(param[0]+param[1]+param[2])
    n2  = N*param[1]/(param[0]+param[1]+param[2])
    a0n = N/(sum0 + n2*param[2]/((param[0]+param[2])*param[0]) + n1*param[1]/((param[0]+param[1])*param[0]) )
    a1n = N/(sum1 + param[0]*n1/((param[0]+param[1])*param[1]) + neg/param[1])
    a2n = N/(sum2 + param[0]*n2/((param[0]+param[2])*param[2]) + neg/param[2])
    return [a0n,a1n,a2n]+param[3:]


def iterPa(N,data,param):
    s1 = param[5]
    s2 = param[6]

    while True :
        s1b = s1
        suma = 0
        for x in data :
            suma += 1/(s1b+x[0]-param[3])
        s1 = N*(param[0]+param[1])/((param[0]+param[1]+1)*suma)
        if abs(s1b-s1)<0.0001 :
            break

    while True:
        s2b = s2
        suma = 0
        for x in data:
            suma += 1/(s2b+x[1]-param[4])
        s2 = N*(param[0]+param[2])/((param[0]+param[2]+1)*suma)
        if abs(s2b-s2)<0.0001 :
            break
    return (param[0:5]+[s1,s2])

# likelihood of these N data points , given these value of parameters
def Q(N,param,data):
    sum0,sum1,sum2 = 0.0,0.0,0.0
    for i in range(N):
        z1 = ( data[i][0] - param[3] )/param[5]
        z2 = ( data[i][1] - param[4] )/param[6]
        sum1 += math.log( 1+z1 )
        sum2 += math.log( 1+z2 )
        sum0 += math.log( 1+max(z1,z2) )

    neg = N*param[0]/(param[0]+param[1]+param[2])
    n1  = N*param[2]/(param[0]+param[1]+param[2])
    n2  = N*param[1]/(param[0]+param[1]+param[2])

    out = 0.0
    out += N*math.log(param[0]*param[1]*param[2])
    out -= ( param[0]*(sum0 + n2*param[2]/((param[0]+param[2])*param[0]) + n1*param[1]/((param[0]+param[1])*param[0])) )
    out -= ( param[1]*(sum1 + param[0]*n1/((param[0]+param[1])*param[1]) + neg/param[1] ) )
    out -= ( param[2]*(sum2 + param[0]*n2/((param[0]+param[2])*param[2]) + neg/param[2] ) )
    return out



if __name__=='__main__':
    args = sys.argv
    configs = open(args[1])   # read the configuration from this file
    counter = 0              # to keep track of the line number
    
    f=open('output_tex_direct_final.txt','wb',0) # The latex output

    # For each line in the configuration file
    for config in configs:  
        counter += 1
        config = config.split(',')
        # config is of the form [ iterCnt, N , alpha0 , alpha1, alpha2 , mu1 , mu2 , sigma1 , sigma2 ]

        iterCnt = int(config[0])
        N = int(config[1])
        AE=[0 for x in range(7)]       # The approximated values
        MSE=[0 for x in range(7)]      # The mean-squared error values
     
        param_real = map(lambda x : float(x) , config[2:] )
        param_init = [1.5,1.5,1.5,1000.0,1000.0,1.5,1.5]
        EMiter = 0.0
        f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
        f.write('\\subsection{{ For parameter : N = {} , $\\alpha_0$ = {} , $\\alpha_1$= {} , $\\alpha_2$ = {} , $\\mu_1$ = {} , $\\mu_2$ = {} , $\\sigma_1$ = {} , $\\sigma_2$ = {} }}\n'.format(N,param_real[0],param_real[1],param_real[2],param_real[3],param_real[4],param_real[5],param_real[6])  )

      # writer = csv.writer(f)
      # writer.writerow(["alpha0","alpha1","alpha2","mu1","mu2","sigma1","sigma2"])
     
        # For each iteration of the EM Algorithm
        for iter in range(iterCnt):
            param = param_init[:]
            data = generate(N,param_real)
            param[3] = min(map(lambda x: x[0] , data ))
            param[4] = min(map(lambda x: x[1] , data ))
     
            # One iteration of the EM Algorithm inside the loop
            while True:
                EMiter += 1.0
                param_old = param[:]
                param = EM1Pas(N,data,param)
                param = iterPa(N,data,param)
                if( abs((Q(N,param,data)-Q(N,param_old,data))/Q(N,param_old,data)) < 0.000001 ):
                    break
     
            for i in range(7):
                AE[i]+=param[i]
                MSE[i]+= (param[i]-param_real[i])**2
     
     
        for i in range(7):
            AE[i]  /= iterCnt
            MSE[i] /= iterCnt

        f.write('No. of iterations  = '+str(EMiter/iterCnt)+'\n\n')
        f.write('\\begin{center}\n')
        f.write('\\begin{tabular}{|| c c c c c c c c ||}\n')
        f.write('\\hline\n')
        f.write('Value & $\\alpha_0 $ & $\\alpha_1$ & $\\alpha_2$ & $\\mu_1$ & $\\mu_2$ & $\\sigma_1$ & $\\sigma_2$ \\\\ \n')
        f.write('\\hline\\hline\n')
        f.write('AE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f}  & {:.5f} \\\\ \n'.format(AE[0],AE[1],AE[2],AE[3],AE[4],AE[5],AE[6]) )
        f.write('MSE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f}  & {:.5f} \\\\ \n'.format(MSE[0],MSE[1],MSE[2],MSE[3],MSE[4],MSE[5],MSE[6]) )
        f.write('\\hline\n')

        f.write('\\end{tabular}\n')
        f.write('\\end{center}\n')
        f.write('\n\n')
        f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

f.close()

     
     
     
     
