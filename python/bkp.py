# Python program for c

# Param vector will have 7 elements as given in the seq:
# param = [alpha0,alpha1,alpha2,mu1,mu2,sigma1,sigma2]
#
import math
import random
import csv

iterCnt = 50  # no. of times em

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

AE=[0 for x in range(7)]
MSE=[0 for x in range(7)]

param_real = [1.0,2.0,2.0,1.0,1.0,1.0,1.0]
param_init = [1.5,1.5,1.5,1000.0,1000.0,1.5,1.5]
EMiter = 0.0

f=open('output.csv','wb')
writer = csv.writer(f)
writer.writerow(["alpha0","alpha1","alpha2","mu1","mu2","sigma1","sigma2"])

# For each iteration of the EM Algorithm
for iter in range(iterCnt):
    N = 250 
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
        if abs(param[0]-param_old[0])<  0.000001 and abs(param[1]-param_old[1])< 0.000001 and abs(param[2]-param_old[2])< 0.000001 :
            break

    for i in range(7):
        AE[i]+=param[i]
        MSE[i]+= (param[i]-param_real[i])**2
    print param
    writer.writerow(param)

f.close()


for i in range(7):
    AE[i]  /= iterCnt
    MSE[i] /= iterCnt

print 'Average Iteration count =  ' +str(EMiter/iterCnt)
print '\n\n'
print 'MSE values: '
print 'Alpha0    : ' + str(MSE[0])
print 'Alpha1    : ' + str(MSE[1])
print 'Alpha2    : ' + str(MSE[2])
print 'Mu1       : ' + str(MSE[3])
print 'Mu2       : ' + str(MSE[4])
print 'Sigma1    : ' + str(MSE[5])
print 'Sigma2    : ' + str(MSE[6])
print '\n\n'

print 'AE values: '
print 'Alpha0    : ' + str(AE[0])
print 'Alpha1    : ' + str(AE[1])
print 'Alpha2    : ' + str(AE[2])
print 'Mu1       : ' + str(AE[3])
print 'Mu2       : ' + str(AE[4])
print 'Sigma1    : ' + str(AE[5])
print 'Sigma2    : ' + str(AE[6])



