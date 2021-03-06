# This is a bootstrap programme to convert csv to latex table files

import os
import re

csvs = []
for f in os.listdir('output'):
    if re.search('.txt',f):
        csvs += [f]

print csvs

# the output of all the lines
f=file('output2.txt','w')

# for each line in the file
for file in csvs :
    f.write('\n\n\n'+file+'\n\n')
    infile = open('output/'+file)
    data = infile.readlines()
    infile.close()
#    data = data[1:]       # remove the header
    f.write('\\begin{center}\n')
    f.write('\\begin{tabular}{||c c c c c c c c || }\n')
    f.write('\\hline\n')
    f.write('Value & $\\alpha_0 $ & $\\alpha_1$ & $\\alpha_2$ & $\\mu_1$ & $\\mu_2$ & $\\sigma_1$ & $\\sigma_2$ \\\\ \n')
    f.write('\\hline\\hline\n')
    print file
    print data[2]
    print data[2].strip().split(',')
#    print (data[2].strip().split(','))
    dt1 = map(float, (data[2].strip().split(','))[1:] )  # for AE
    dt2 = map(float, (data[3].strip().split(','))[1:] )  # for MSE
    f.write('AE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\ \n'.format(dt1[0],dt1[1],dt1[2],dt1[3],dt1[4],dt1[5],dt1[6]))
    f.write('MSE & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\ \n'.format(dt2[0],dt2[1],dt2[2],dt2[3],dt2[4],dt2[5],dt2[6]))
    f.write('\\hline\n')

    f.write('\\end{tabular}\n')
    f.write('\\end{center}\n')
    f.write('\n\n\n')
        

f.close()


