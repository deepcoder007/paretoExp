# This is a bootstrap programme to convert csv to latex table files

import os
import re

pars = []
f = open('rpars.txt')
pars=f.readlines()


csvs = []
for f in os.listdir('output'):
    if re.search('.txt',f):
        csvs += [f]


# the output of all the lines
f=file('output3.txt','w')

# for each line in the file
for file in csvs :
    idx = 0
    print file
    try:
        idx = int(file.strip()[7:9])
    except :
        idx = int(file.strip()[7:8])
    print idx

    f.write('\n\n\n'+file+'\n\n')
    infile = open('output/'+file)
    data = infile.readlines()
    infile.close()
    curr =  pars[idx-1].strip().split(',')
    nn = int(curr[1])
    alpha0 = float(curr[2])
    alpha1 = float(curr[3])
    alpha2 = float(curr[4])
    mu1 = float(curr[5])
    mu2 = float(curr[6])
    sigma1 = float(curr[7])
    sigma2 = float(curr[8])
    f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    f.write('\\subsection{{ For parameter: N = {}, $\\alpha_0$ = {} , $\\alpha_1$ = {} , $\\alpha_2$ = {} , $\\mu_1$ = {} , $\\mu_2$ = {} , $\\sigma_1$ = {} , $\\sigma_2$ = {} }}\n '.format(nn,alpha0,alpha1,alpha2,mu1,mu2,sigma1,sigma2) )

    f.write('\\begin{center}\n')
    f.write('\\begin{tabular}{||c c c c c c c c || }\n')
    f.write('\\hline\n')
    f.write('Value & $\\alpha_0 $ & $\\alpha_1$ & $\\alpha_2$ & $\\mu_1$ & $\\mu_2$ & $\\sigma_1$ & $\\sigma_2$ \\\\ \n')
    f.write('\\hline\\hline\n')
#    print (data[2].strip().split(','))
    dt1 = map(float, (data[2].strip().split(','))[1:] )  # for AE
    dt2 = map(float, (data[3].strip().split(','))[1:] )  # for MSE
    f.write('AE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\ \n'.format(dt1[0],dt1[1],dt1[2],dt1[3],dt1[4],dt1[5],dt1[6]))
    f.write('MSE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\ \n'.format(dt2[0],dt2[1],dt2[2],dt2[3],dt2[4],dt2[5],dt2[6]))
    f.write('\\hline\n')

    f.write('\\end{tabular}\n')
    f.write('\\end{center}\n')
    f.write('\n\n\n')
    f.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        

f.close()


