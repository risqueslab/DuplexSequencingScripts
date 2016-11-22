#!/usr/bin/python

import matplotlib.pyplot as plt
from sys import argv

script, pileup_file = argv
prefix = pileup_file.split('.')[0]
pileup = open(pileup_file, 'r+')

pos = []
depth = []

for line in pileup:
	pos.append(int(line.split()[1]))
	depth.append(int(line.split()[3]))

plt.bar(pos, depth, alpha=0.5)

plt.ylabel('# of DCS Reads')

plt.xlabel('Genomic Coordinate')

plt.title(prefix + ' DCS Reads over genomic coordinates')

plt.savefig(prefix + '.DCS_Depth_by_pos.png', bbox_inches='tight')
