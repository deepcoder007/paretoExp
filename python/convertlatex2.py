# This is a bootstrap programme to convert csv to latex table files

import os
import re

csvs = []
for f in os.listdir('output'):
    if re.search('.csv',f):
        csvs += [f]

print csvs

# the output of all the lines
f=file('output.txt','w')

# for each line in the file
for file in csvs :
    infile = open('output/'+file)
    data = infile.readlines()
    infile.close()
    data = data[1:]       # remove the header
    f.write('\\begin{center}\n')
    f.write('\\begin{tabular}{|| c c c c c c c || }\n')
    f.write('\\hline\n')
    f.write('$\\alpha_0 $ & $\\alpha_1$ & $\\alpha_2$ & $\\mu_1$ & $\\mu_2$ & $\\sigma_1$ & $\\sigma_2$ \\\\ \n')
    f.write('\\hline\\hline\n')
    for x in data :
        dt = map(float, x.strip().split(','))
        f.write('{:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} & {:.5f} \\\\ \n'.format(dt[0],dt[1],dt[2],dt[3],dt[4],dt[5],dt[6]))
 #       f.write(str(dt[0])+' & '+str(dt[1]) +' & '+str(dt[2])+' & '+str(dt[3])+' & '+str(dt[4])+' & '+str(dt[5])+' & '+str(dt[6])+'\\\\ \n')
        f.write('\\hline\n')
    f.write('\\end{tabular}\n')
    f.write('\\end{center}\n')
    f.write('\n\n\n')
        

f.close()


