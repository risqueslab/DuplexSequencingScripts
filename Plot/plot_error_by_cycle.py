#!/usr/bin/python

import matplotlib.pyplot as plt
from sys import argv

# Create a plot of mutations per read cycle from a text file created by 
# GATK ErrorRateByCycle tool. 
#
# Usage: >> python plot_error_by_cycle.py error_per_cycle.txt
# Last modified by Dana Nachmanson 9/2/16

script, error_file = argv
prefix = error_file.split('.')[0]
file = open(error_file, 'r+')

cycle = []
qual = []
count = 0 

for line in file:
	if count > 3 :
		if line.strip():
			cycle.append(int(line.split()[1]))
			qual.append(int(line.split()[2]))
		else:
			break
	count += 1

plt.bar(cycle, qual, alpha=0.5)

plt.ylabel('Mutations')

plt.xlabel('Read Cycle #')

plt.title(prefix + ' Mutations by Read Cycle')

plt.savefig(prefix + '.Muts_by_Cycle.png', bbox_inches='tight')
