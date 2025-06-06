import sys
import numpy as np
import scipy.stats as stats
# total arguments
name=sys.argv[1]
op=sys.argv[2]
r1=int(sys.argv[3])-1
r2=int(sys.argv[4])-1
chr=int(sys.argv[5])
p1=int(sys.argv[6])
p2=int(sys.argv[7])
with open(f'{name}', 'r') as f:
    data = [line.strip().split() for line in f.readlines()]
    # Filter out rows containing 0
    #filtered_data = [row for row in data if float(row[0]) > 0.025 and float(row[1]) > 0.025]
    data1 = [float(row[r1]) for row in data]
    data2 = [float(row[r2]) for row in data]

statistic, pvalue = stats.ttest_rel(data1, data2)

with open(f'{op}', 'a') as f:
    f.write(f'{chr}\t{p1}\t{p2}\t{pvalue}\n')
